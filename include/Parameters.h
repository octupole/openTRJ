/*
 * Parameters.h
 *
 *  Created on: Mar 12, 2012
 *      Author: marchi
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <string>
#include <cmath>

namespace Parameters{
	const double pi=3.14159265358979323846,
			twopi=2.0*pi,
			avogad=6.0225e23,boltz=1.38054e-23,
			gascon=8.3143,
    		planck=6.6256e-34,
    		elechg=1.602e-19,
    		epso=8.854e-12,
    		boxl=2.0,
    		unitm=1.0/(avogad*1000.0),
    		unitl=1.0e-9,
    		unitt=1.0e-15,
    		unite=unitm*(unitl/unitt)*(unitl/unitt),
    		unitc=4.0*pi*epso*unitl*unite/(elechg*elechg),
    		unitp=(unite/unitl*unitl*unitl)/1.0e6,
    		efact=unite*avogad,
    		efact_nm=100*unite*avogad,
    		hartree=4.35981 * 1.0e-18,
    		unitepot=unite/sqrt(unitc)/hartree,
    		lbohr=0.52917706,
    		unitefield=unitepot*lbohr,
    		unitfield=elechg/unitl/(4.0*pi*epso),
			unit_mm=avogad*(elechg*elechg/(epso*unitl)/1000.0),
    		kT300=gascon*300.0/1000.0;
}




#endif /* PARAMETERS_H_ */
