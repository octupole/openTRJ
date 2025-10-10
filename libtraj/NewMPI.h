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
const size_t one{1};

template <>
class Intracom<false>{
public:
	Intracom(void *){};
	size_t Get_size(){return one;}
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
	MPI_Comm comm;

public:
	NewMPI(int & argc,char ** & argv){
		MPI_Init(&argc,&argv);
		comm=MPI_COMM_WORLD;
		int size, rank;
		MPI_Comm_size(comm, &size);
		MPI_Comm_rank(comm, &rank);
		if(size > 1){
			if(rank){
				std::cout.setstate(std::ios_base::badbit);
				fclose(stdout);
				fclose(stderr);
			}
			std::cout << "\n" << std::endl;
			std::cout << " ------ Parallel run with " << size << " CPUS " << std::endl;
			is_parallel=true;
		}
	};
	NewMPI(){
		MPI_Init(NULL, NULL);
		comm=MPI_COMM_WORLD;
		int size, rank;
		MPI_Comm_size(comm, &size);
		MPI_Comm_rank(comm, &rank);
		if(size > 1){

			if(rank){
				std::cout.setstate(std::ios_base::badbit);
				fclose(stdout);
				fclose(stderr);
			}
			std::cout << "\n" << std::endl;
			std::cout << " ------ Parallel run with " << size << " CPUS " << std::endl;
			is_parallel=true;
		}
	};
	size_t Get_Rank(){
		int rank;
		MPI_Comm_rank(comm, &rank);
		return rank;
	}
	size_t Get_Size(){
		int size;
		MPI_Comm_size(comm, &size);
		return size;
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

	void Barrier(){MPI_Barrier(comm);}
	double Time(){return MPI_Wtime();};

	template <typename T>
	void Send(int dest, int dim, int tag, T * sbuffer);

	template <typename T>
	void Recv(int source,int dim, int tag, T * rbuffer);

	MPI_Comm * Communicator(){ return &comm;}
	void Finalize(){MPI_Finalize();};
	virtual ~NewMPI(){MPI_Finalize();};
};

} /* namespace Parallel */


namespace Parallel {

template<>
inline void NewMPI::Gather<double>(vector<double> & bufferIn, vector<double> & bufferOut){
	int nIn=static_cast<int> (bufferIn.size());
	int nbufferIn,nbufferOut;
	int rank;
	MPI_Comm_rank(comm, &rank);
	nbufferIn=nIn;
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(&nIn,NULL,1,MPI_INT,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,&nIn,1,MPI_INT,MPI_SUM,0,comm);
	MPI_Bcast(&nIn,1,MPI_INT,0,comm);
	nbufferOut=nIn;
	bufferOut.clear();
	bufferOut=vector<double>(nbufferOut,0.0);
	MPI_Gather(&bufferIn[0],nbufferIn,MPI_DOUBLE,&bufferOut[0],nbufferIn,MPI_DOUBLE,0,comm);
}

template<>
inline void NewMPI::Broadcast<double>(double * buffer,const int nbuffer){
	MPI_Bcast(buffer,nbuffer,MPI_DOUBLE,0,comm);
};

template<>
inline void NewMPI::Broadcast<int>(int * buffer,const int nbuffer){
	MPI_Bcast(buffer,nbuffer,MPI_INT,0,comm);
	};
template<>
inline void NewMPI::Broadcast<float>(float * buffer,const int nbuffer){
	MPI_Bcast(buffer,nbuffer,MPI_FLOAT,0,comm);
};

template<>
inline void NewMPI::Broadcast<char>(char * buffer,const int nbuffer){
	MPI_Bcast(buffer,nbuffer,MPI_CHAR,0,comm);
	};

template<>
inline void NewMPI::ReduceSum<float>(float * buffer,const int nbuffer){
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(buffer,NULL,nbuffer,MPI_FLOAT,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_FLOAT,MPI_SUM,0,comm);
	MPI_Barrier(comm);
	};
template<>
inline void NewMPI::ReduceSum<double>(double * buffer,const int nbuffer){
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(buffer,NULL,nbuffer,MPI_DOUBLE,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_DOUBLE,MPI_SUM,0,comm);
	MPI_Barrier(comm);
};
template<>
inline void NewMPI::ReduceSum<Complex>(Complex * buffer,const int nbuffer){
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(buffer,NULL,nbuffer,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
	MPI_Barrier(comm);
};
template<>
inline void NewMPI::ReduceSum<double>(vector<double> & buffer0){
	int rank;
	MPI_Comm_rank(comm, &rank);
	const int nbuffer=buffer0.size();
	double * buffer;
	buffer=new double[nbuffer];
	for(int n=0;n<nbuffer;n++)
		buffer[n]=buffer0[n];
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(buffer,NULL,nbuffer,MPI_DOUBLE,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_DOUBLE,MPI_SUM,0,comm);
	for(int n=0;n<nbuffer;n++)
		buffer0[n]=buffer[n];
	MPI_Barrier(comm);
	};
template<>
inline void NewMPI::ReduceSum<size_t>(size_t * buffer,const int nbuffer){
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(buffer,NULL,nbuffer,MPI_INT,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_INT,MPI_SUM,0,comm);
	};
template<>
inline void NewMPI::ReduceSum<int>(int * buffer,const int nbuffer){
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Barrier(comm);
	if(rank) MPI_Reduce(buffer,NULL,nbuffer,MPI_INT,MPI_SUM,0,comm);
	else MPI_Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_INT,MPI_SUM,0,comm);
	};
template<>
inline void NewMPI::AllGather<int>(int n, int * sendbuff, int * recvbuf){
	MPI_Allgather(sendbuff,n, MPI_INT, recvbuf, n, MPI_INT,comm);
}
template<>
inline void NewMPI::AllGather<double>(int n, double * sendbuff, double * recvbuf){
	MPI_Allgather(sendbuff,n, MPI_DOUBLE, recvbuf, n, MPI_DOUBLE,comm);
}
template<>
inline void NewMPI::Send<float>(int dest, int nbuffer, int tag, float * x){
	MPI_Send(x,nbuffer,MPI_FLOAT,dest,tag,comm);
}
template<>
inline void NewMPI::Send<double>(int dest, int nbuffer, int tag, double * x){
	MPI_Send(x,nbuffer,MPI_DOUBLE,dest,tag,comm);
}
template<>
inline void NewMPI::Recv<float>(int source, int nbuffer, int tag, float * x){
	MPI_Status status;
	MPI_Recv(x,nbuffer,MPI_FLOAT,source,tag,comm,&status);
}
template<>
inline void NewMPI::Recv<double>(int source, int nbuffer, int tag, double * x){
	MPI_Status status;
	MPI_Recv(x,nbuffer,MPI_DOUBLE,source,tag,comm,&status);
}

}
#endif /* LIBTRAJ_NEWMPI_H_ */
