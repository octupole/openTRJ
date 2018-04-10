/*
 * CenterMassWigner.h
 *
 *  Created on: Jun 1, 2015
 *      Author: marchi
 */

#ifndef SRC_CENTERMASSWIGNER_H_
#define SRC_CENTERMASSWIGNER_H_

#include "CenterMassBW3.h"

template <typename T>
class CenterMassWigner: public CenterMassBW3<T> {
public:
	CenterMassWigner();
	virtual void DiffQ(ofstream &,const int , const int , const int);
	virtual void Diffusion(ofstream & fout){};
	virtual ~CenterMassWigner();
	virtual void GofR(ofstream & fout){};
};

#endif /* SRC_CENTERMASSWIGNER_H_ */
