/*
 * VoronoiMicelles.h
 *
 *  Created on: Nov 19, 2017
 *      Author: marchi
 */

#ifndef VORONOIMICELLES_H_
#define VORONOIMICELLES_H_

#include <Voronoi.h>
namespace Voro{
class VoronoiMicelles: public Voro::Voronoi {
protected:
	virtual void WriteIt(std::ofstream &);
	void ReadIt(std::ifstream &){};

	void __compShell();
	void __searchNeighs(int,int);
	void __computeAggregate();
	VoronoiMicelles(){};

public:

	VoronoiMicelles(ifstream & f){};
	VoronoiMicelles(Topol & x,bool y, Parallel::NewMPI * curr): VoronoiMicelles(x,y){}
	VoronoiMicelles(Topol &,bool);
	virtual void WriteLastJSON(std::ofstream & fout){};

	virtual void getData();
	virtual ~VoronoiMicelles();
};
}
#endif /* VORONOIMICELLES_H_ */
