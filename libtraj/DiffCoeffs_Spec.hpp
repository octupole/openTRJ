/*
 * DiffCoeffs_Spec.hpp
 *
 *  Created on: Jul 12, 2011
 *      Author: marchi
 */

#ifndef DIFFCOEFFS_SPEC_HPP_
#define DIFFCOEFFS_SPEC_HPP_


template <> double DiffCoeffs<2>::coef[1]={1.0};
template <> double DiffCoeffs<4>::coef[2]={4.0/3.0,-1.0/6.0};
template <> double DiffCoeffs<6>::coef[3]={3.0/2.0,-3.0/10.0,1.0/30.0};
template <> double DiffCoeffs<8>::coef[4]={8.0/5.0,-2.0/5.0,8.0/105.0,-2.0/280.0};


#endif /* DIFFCOEFFS_SPEC_HPP_ */
