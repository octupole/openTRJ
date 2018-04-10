/*
 * CMPick.h
 *
 *  Created on: Jun 3, 2015
 *      Author: marchi
 */

#ifndef SRC_CMPICK_H_
#define SRC_CMPICK_H_
#include "CenterMass.h"
template <typename T>
class CenterMassBW3;
template <typename T>
class CenterMassWigner;
template <typename T>
class CenterMassMicelles;

// Functor to figure out which diffusion parameters should be computed
template <typename T>
class CMPick{
public:
	CenterMass_t operator ()(CenterMass<T> *);
	virtual ~CMPick();
};

#endif /* SRC_CMPICK_H_ */
