/*
 * Ftypedefs.h
 *
 *  Created on: May 12, 2012
 *      Author: marchi
 */

#ifndef LINKEDCELLS_SRC_FTYPEDEFS_H_
#define LINKEDCELLS_SRC_FTYPEDEFS_H_
#include <cstdlib>

#define DIM 3
#define XX 0
#define YY 1
#define ZZ 2
enum {Slt, Sol, Ions};
const double unit_nm=0.1;
const double unit_amu=1.66053892e-27;
const double unit_rho=unit_amu/1.0e-24;


namespace Typedefs{
typedef int     	atom_id;	/* To indicate an atoms id         */




typedef float           real;

typedef real        	rvec[DIM];

typedef double       	dvec[DIM];

typedef real	    	matrix[DIM][DIM];


typedef struct {
  int nr;			/* The number of blocks			*/
  atom_id *index;		/* Array of indices (dim: nr+1) 	*/
  int nalloc_index;             /* The allocation size for index        */
} t_blockx;
}
#endif /* nodef HAVE_TYPEDEFS_H */

