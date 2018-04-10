/*
 * VoronoiBinary.h
 *
 *  Created on: Dec 2, 2017
 *      Author: marchi
 */

#ifndef VORONOI_VORONOIBINARY_H_
#define VORONOI_VORONOIBINARY_H_

#include "VoronoiMicelles.h"
#include "VoronoiMicellesJSON.h"

namespace Voro {

template <typename T>
class VoronoiBinary: public T {
	void WriteIt(std::ofstream &);
	void ReadIt(std::ifstream &);
public:
	VoronoiBinary(ifstream &);
	VoronoiBinary(ofstream &,Topol &,bool, Parallel::NewMPI *);

	virtual ~VoronoiBinary();
};

} /* namespace Voro */

#endif /* VORONOI_VORONOIBINARY_H_ */
