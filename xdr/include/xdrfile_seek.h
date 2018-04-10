/*
 * my.h
 *
 *  Created on: Feb 13, 2016
 *      Author: marchi
 */

#ifndef XDRFILE_SEEK_H_
#define XDRFILE_SEEK_H_
#ifdef __cplusplus
extern "C"
{
#endif

#include "xdrfile.h"
#define XDR_INT_SIZE 4
#define XTC_MAGIC 1995

static const int header_size = 16;
#ifndef XDR_H_
#define XDR_H_
typedef struct XDR XDR;
#endif

  static off_t xtc_get_next_frame_start(FILE *fp, XDR *xdrs, int natoms);
  static int xtc_get_current_frame_number(FILE *fp, XDR *xdrs, int natoms, int * bOK);
  static int xtc_at_header_start(FILE *fp, XDR *xdrs,
                               int natoms, int * timestep, float * time);
  static float xtc_get_current_frame_time(FILE *fp, XDR *xdrs, int natoms, int * bOK);
  int xtc_get_next_frame_number(FILE *fp, XDR *xdrs, int natoms);
  int xdr_xtc_seek_frame(int frame, FILE *fp, XDR *xdrs, int natoms);
  int xdr_xtc_get_last_frame_number(FILE *fp, XDR *xdrs, int natoms, int * bOK);
  float xdr_xtc_get_last_frame_time(FILE *fp, XDR *xdrs, int natoms, int * bOK);

#ifdef __cplusplus
}
#endif
#endif /* XDRFILE_SEEK_H_ */

