/*
 * my.c
 *
 *  Created on: Feb 13, 2016
 *      Author: marchi
 */
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xdrfile_seek.h"


int xtc_get_next_frame_number(FILE *fp, XDR *xdrs, int natoms)
{
    off_t off;
    int       step;
    float     time;
    int       ret;

    if ((off = ftell(fp)) < 0)
    {
        return -1;
    }

    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            if (fseek(fp, off, SEEK_SET))
            {
                return -1;
            }

            return step;
        }
        else if (ret == -1)
        {
            if (fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
        }
    }
    return -1;
}
static int xtc_at_header_start(FILE *fp, XDR *xdrs,
                               int natoms, int * timestep, float * time)
{
    int       i_inp[3];
    float     f_inp[10];
    int       i;
    off_t off;


    if ((off = ftell(fp)) < 0)
    {
        return -1;
    }
    /* read magic natoms and timestep */
    for (i = 0; i < 3; i++)
    {
        if (!xdr_int(xdrs, &(i_inp[i])))
        {
            fseek(fp, off+XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* quick return */
    if (i_inp[0] != XTC_MAGIC)
    {
        if (fseek(fp, off+XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        return 0;
    }
    /* read time and box */
    for (i = 0; i < 10; i++)
    {
        if (!xdr_float(xdrs, &(f_inp[i])))
        {
            fseek(fp, off+XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* Make a rigourous check to see if we are in the beggining of a header
       Hopefully there are no ambiguous cases */
    /* This check makes use of the fact that the box matrix has 3 zeroes on the upper
       right triangle and that the first element must be nonzero unless the entire matrix is zero
     */
    if (i_inp[1] == natoms &&
        ((f_inp[1] != 0 && f_inp[6] == 0) ||
         (f_inp[1] == 0 && f_inp[5] == 0 && f_inp[9] == 0)))
    {
        if (fseek(fp, off+XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        *time     = f_inp[0];
        *timestep = i_inp[2];
        return 1;
    }
    if (fseek(fp, off+XDR_INT_SIZE, SEEK_SET))
    {
        return -1;
    }
    return 0;
}




int xdr_xtc_seek_frame(int frame, FILE *fp, XDR *xdrs, int natoms)
{
	off_t low = 0;
	off_t high, pos;


	/* round to 4 bytes */
	int        fr;
	off_t  offset;
	if (fseek(fp, 0, SEEK_END)) return -1;

	if ((high = ftell(fp)) < 0)return -1;

	/* round to 4 bytes  */
	high  /= XDR_INT_SIZE;
	high  *= XDR_INT_SIZE;
	offset = ((high/2)/XDR_INT_SIZE)*XDR_INT_SIZE;

	if (fseek(fp, offset, SEEK_SET)) return -1;

	while (1){
		fr = xtc_get_next_frame_number(fp, xdrs, natoms);
		if (fr < 0)return -1;

		if (fr != frame && llabs(low-high) > header_size){
			if (fr < frame) low = offset;
			else high = offset;

			/* round to 4 bytes */
			offset = (((high+low)/2)/4)*4;

			if (fseek(fp, offset, SEEK_SET)) return -1;
		}
		else break;
	}
	if (offset <= header_size) offset = low;

	if (fseek(fp, offset, SEEK_SET)) return -1;

	if ((pos = xtc_get_next_frame_start(fp, xdrs, natoms)) < 0){
		/* we probably hit an end of file */
		return -1;
	}

	if (fseek(fp, pos, SEEK_SET)) return -1;

	return 0;
}
static off_t xtc_get_next_frame_start(FILE *fp, XDR *xdrs, int natoms)
{
    off_t res;
    int       ret;
    int       step;
    float     time;
    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            if ((res = ftell(fp)) >= 0)
            {
                return res - XDR_INT_SIZE;
            }
            else
            {
                return res;
            }
        }
        else if (ret == -1)
        {
            return -1;
        }
    }
    return -1;
}
int xdr_xtc_get_last_frame_number(FILE *fp, XDR *xdrs, int natoms, int * bOK)
{
    int        frame;
    off_t  off;
    int        res;
    *bOK = 1;

    if ((off = ftell(fp)) < 0)
    {
        *bOK = 0;
        return -1;
    }


    if (fseek(fp, -3*XDR_INT_SIZE, SEEK_END))
    {
        *bOK = 0;
        return -1;
    }

    frame = xtc_get_current_frame_number(fp, xdrs, natoms, bOK);
    if (!bOK)
    {
        return -1;
    }


    if (fseek(fp, off, SEEK_SET))
    {
        *bOK = 0;
        return -1;
    }

    return frame;
}
static int xtc_get_current_frame_number(FILE *fp, XDR *xdrs, int natoms, int * bOK)
{
    off_t off;
    int       ret;
    int       step;
    float     time;
    *bOK = 0;

    if ((off = ftell(fp)) < 0)
    {
        return -1;
    }


    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = 1;
            if (fseek(fp, off, SEEK_SET))
            {
                *bOK = 0;
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;

        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (fseek(fp, -2*XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
    return -1;
}
static float xtc_get_current_frame_time(FILE *fp, XDR *xdrs, int natoms, int * bOK)
{
    off_t off;
    int       step;
    float     time;
    int       ret;
    *bOK = 0;

    if ((off = ftell(fp)) < 0)
    {
        return -1;
    }

    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = 1;
            if (fseek(fp, off, SEEK_SET))
            {
                *bOK = 0;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (fseek(fp, -2*XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
    return -1;
}
float xdr_xtc_get_last_frame_time(FILE *fp, XDR *xdrs, int natoms, int * bOK)
{
    float      time;
    off_t  off;
    int        res;
    *bOK = 1;
    off  = ftell(fp);
    if (off < 0)
    {
        *bOK = 0;
        return -1;
    }

    if ( (res = fseek(fp, -3*XDR_INT_SIZE, SEEK_END)) != 0)
    {
        *bOK = 0;
        return -1;
    }

    time = xtc_get_current_frame_time(fp, xdrs, natoms, bOK);
    if (!(*bOK))
    {
        return -1;
    }

    if ( (res = fseek(fp, off, SEEK_SET)) != 0)
    {
        *bOK = 0;
        return -1;
    }
    return time;
}

