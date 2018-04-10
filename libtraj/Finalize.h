/*
 * Finalize.h
 *
 *  Created on: Mar 4, 2016
 *      Author: marchi
 */

#ifndef LIBTRAJ_FINALIZE_H_
#define LIBTRAJ_FINALIZE_H_
#include "NewMPI.h"
namespace Finale {

class Finalize {
public:
	static Parallel::NewMPI * CurrMPI;
	Finalize(Parallel::NewMPI * y){CurrMPI=y;};
	static void Final(){if(CurrMPI) {
		CurrMPI->Barrier();
		CurrMPI->Finalize();
	}
	exit(0);}
	virtual ~Finalize();
};

} /* namespace Finale */

#endif /* LIBTRAJ_FINALIZE_H_ */
