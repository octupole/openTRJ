/*
 * CMPick.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: marchi
 */

#include "CMPick.h"
#include "CenterMassWigner.h"
#include "CenterMassMicelles.h"

template <typename T>
CenterMass_t CMPick<T>::operator ()(CenterMass<T> * p){
	size_t t=0;
	if(dynamic_cast<CenterMassWigner<T> *> (p)){
		return diffq;
	} else if(dynamic_cast<CenterMassMicelles<T>*> (p)){
		return diffMicelles;
	} else if(dynamic_cast<CenterMassBW3<T>*> (p)){
		return diffk;
	} else{
		return diff0;
	}
}

template <typename T>
CMPick<T>::~CMPick() {
	// TODO Auto-generated destructor stub
}
template class CMPick<float>;
template class CMPick<double>;
