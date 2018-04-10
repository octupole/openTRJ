/*
 * myEnums.hpp
 *
 *  Created on: Feb 24, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_MYENUMS_HPP_
#define TRJLIB_MYENUMS_HPP_

namespace Enums{
enum Padding {zero, avgDensity, myDensity, Periodic};
enum myAtoms {Reference, Selection, fftPadding};
enum Compute {SQ, SAXS, SANS};
enum myWriteOptions { noJSON, JSON};
}



#endif /* TRJLIB_MYENUMS_HPP_ */
