/*
 * NewMPI.h
 *
 *  Created on: Jan 28, 2016
 *      Author: marchi
 */

#ifndef LIBTRAJ_NEWMPI_H_
#define LIBTRAJ_NEWMPI_H_

#include <iostream>
#include <vector>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <complex>


namespace MPI{
template <bool MMPI>
class Intracom;
}
#ifdef HAVE_MPI
#include <mpi.h>
#else
namespace MPI{
static void * COMM_WORLD;
static void * INT;
static void * FLOAT;
static void * SUM;
static void * DOUBLE;
static void * DOUBLE_COMPLEX;
static void * IN_PLACE;
static void * CHAR;
using Intracomm=Intracom<false>;
inline static void Init(...){};
inline static void Finalize(){};
inline static double Wtime(){
	timeval tim;
	double temp=tim.tv_sec+(tim.tv_usec/1000000.0);
	return temp;
	}
}
#endif

using std::vector;
using std::string;
using std::cout;
using std::endl;

namespace MPI{
const size_t zero{0};

template <>
class Intracom<false>{
public:
	Intracom(void *){};
	size_t Get_size(){return zero;}
	size_t Get_rank(){return zero;}
	void Barrier(){};
	void Reduce(...){};
	void Bcast(...){};
	void Gather(...){};
	void Allgather(...){};
	void Send(...){};
	void Recv(...){};
	Intracom(){};
	virtual ~Intracom(){};
};


}

namespace Parallel{
typedef std::complex<double> Complex;
const bool DEBUG{false};
class NewMPI {
	bool is_parallel{false};
	MPI::Intracomm comm{0};

public:
	NewMPI(int & argc,char ** & argv){
		MPI::Init(argc,argv);
		comm=MPI::COMM_WORLD;
		if(comm.Get_size()){
			if(comm.Get_rank()){
				std::cout.setstate(std::ios_base::badbit);
				fclose(stdout);
				fclose(stderr);
			}
			std::cout << "\n" << std::endl;
			std::cout << " ------ Parallel run with " << comm.Get_size() << " CPUS " << std::endl;
			is_parallel=true;
		}
	};
	NewMPI(){
		MPI::Init();
		comm=MPI::COMM_WORLD;
		if(comm.Get_size()){

			if(comm.Get_rank()){
				std::cout.setstate(std::ios_base::badbit);
				fclose(stdout);
				fclose(stderr);
			}
			std::cout << "\n" << std::endl;
			std::cout << " ------ Parallel run with " << comm.Get_size() << " CPUS " << std::endl;
			is_parallel=true;
		}
	};
	size_t Get_Rank(){
		return comm.Get_rank();
	}
	size_t Get_Size(){
		return comm.Get_size();
	}
	bool AmI_Parallel(){return is_parallel;};

	template <class T>
	void Broadcast(T *, const int);

	template <class T>
	void ReduceSum(T *, const int);

	template <class T>
	void ReduceSum(vector<T> &);

	template <class T>
	void Gather(vector<T> & , vector<T> &);

	template <typename T>
	void AllGather(int, T *, T * );

	void Barrier(){comm.Barrier();}
	double Time(){return MPI::Wtime();};

	template <typename T>
	void Send(int dest, int dim, int tag, T * sbuffer);

	template <typename T>
	void Recv(int source,int dim, int tag, T * rbuffer);

	MPI::Intracomm * Communicator(){ return &comm;}
	void Finalize(){MPI::Finalize();};
	virtual ~NewMPI(){MPI::Finalize();};
};

} /* namespace Parallel */


namespace Parallel {

template<>
inline void NewMPI::Gather<double>(vector<double> & bufferIn, vector<double> & bufferOut){
	int nIn=static_cast<int> (bufferIn.size());
	int nbufferIn,nbufferOut;
	nbufferIn=nIn;
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(&nIn,NULL,1,MPI::INT,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,&nIn,1,MPI::INT,MPI::SUM,0);
	comm.Bcast(&nIn,1,MPI::INT,0);
	nbufferOut=nIn;
	bufferOut.clear();
	bufferOut=vector<double>(nbufferOut,0.0);
	comm.Gather(&bufferIn[0],nbufferIn,MPI::DOUBLE,&bufferOut[0],nbufferIn,MPI::DOUBLE,0);
}

template<>
inline void NewMPI::Broadcast<double>(double * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI::DOUBLE,0);
};

template<>
inline void NewMPI::Broadcast<int>(int * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI::INT,0);
	};
template<>
inline void NewMPI::Broadcast<float>(float * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI::FLOAT,0);
};

template<>
inline void NewMPI::Broadcast<char>(char * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI::CHAR,0);
	};

template<>
inline void NewMPI::ReduceSum<float>(float * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::FLOAT,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::FLOAT,MPI::SUM,0);
	comm.Barrier();
	};
template<>
inline void NewMPI::ReduceSum<double>(double * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::DOUBLE,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::DOUBLE,MPI::SUM,0);
	comm.Barrier();
};
template<>
inline void NewMPI::ReduceSum<Complex>(Complex * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::DOUBLE_COMPLEX,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::DOUBLE_COMPLEX,MPI::SUM,0);
	comm.Barrier();
};
template<>
inline void NewMPI::ReduceSum<double>(vector<double> & buffer0){
	const int nbuffer=buffer0.size();
	double * buffer;
	buffer=new double[nbuffer];
	for(int n=0;n<nbuffer;n++)
		buffer[n]=buffer0[n];
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::DOUBLE,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::DOUBLE,MPI::SUM,0);
	for(int n=0;n<nbuffer;n++)
		buffer0[n]=buffer[n];
	comm.Barrier();
	};
template<>
inline void NewMPI::ReduceSum<int>(int * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI::INT,MPI::SUM,0);
	else comm.Reduce(MPI::IN_PLACE,buffer,nbuffer,MPI::INT,MPI::SUM,0);
	};
template<>
inline void NewMPI::AllGather<int>(int n, int * sendbuff, int * recvbuf){
	comm.Allgather(sendbuff,n, MPI::INT, recvbuf, n, MPI::INT);
}
template<>
inline void NewMPI::AllGather<double>(int n, double * sendbuff, double * recvbuf){
	comm.Allgather(sendbuff,n, MPI::DOUBLE, recvbuf, n, MPI::DOUBLE);
}
template<>
inline void NewMPI::Send<float>(int dest, int nbuffer, int tag, float * x){
	comm.Send(x,nbuffer,MPI::FLOAT,dest,tag);
}
template<>
inline void NewMPI::Send<double>(int dest, int nbuffer, int tag, double * x){
	comm.Send(x,nbuffer,MPI::DOUBLE,dest,tag);
}
template<>
inline void NewMPI::Recv<float>(int source, int nbuffer, int tag, float * x){
	comm.Recv(x,nbuffer,MPI::FLOAT,source,tag);
}
template<>
inline void NewMPI::Recv<double>(int source, int nbuffer, int tag, double * x){
	comm.Recv(x,nbuffer,MPI::DOUBLE,source,tag);
}

}
#endif /* LIBTRAJ_NEWMPI_H_ */
