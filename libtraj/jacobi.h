/*
 * jacobi.h
 *
 *  Created on: Feb 13, 2016
 *      Author: marchi
 */

#ifndef LIBTRAJ_JACOBI_H_
#define LIBTRAJ_JACOBI_H_
#include <cstdlib>
#include <string>
#include <iostream>

using namespace std;
static inline void do_rotate(double **a, int i, int j, int k, int l, double tau, double s)
{
    double g, h;
    g       = a[i][j];
    h       = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

void jacobi(double **a, int n, double d[], double **v, int *nrot){
    int    j, i;
    int    iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;
    b=new double[n];
    z=new double[n];
    for (ip = 0; ip < n; ip++)
    {
    	for (iq = 0; iq < n; iq++)
    	{
    		v[ip][iq] = 0.0;
    	}
    	v[ip][ip] = 1.0;
    }
    for (ip = 0; ip < n; ip++)
    {
    	b[ip] = d[ip] = a[ip][ip];
    	z[ip] = 0.0;
    }
    *nrot = 0;
    try{
    	for (i = 1; i <= 50; i++)
    	{
    		sm = 0.0;
    		for (ip = 0; ip < n-1; ip++)
    		{
    			for (iq = ip+1; iq < n; iq++)
    			{
    				sm += fabs(a[ip][iq]);
    			}
    		}
    		if (sm == 0.0)
    		{
    			delete [] z;
    			delete [] b;
    			return;
    		}
    		if (i < 4)
    		{
    			tresh = 0.2*sm/(n*n);
    		}
    		else
    		{
    			tresh = 0.0;
    		}
    		for (ip = 0; ip < n-1; ip++)
    		{
    			for (iq = ip+1; iq < n; iq++)
    			{
    				g = 100.0*fabs(a[ip][iq]);
    				if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
    						&& fabs(d[iq])+g == fabs(d[iq]))
    				{
    					a[ip][iq] = 0.0;
    				}
    				else if (fabs(a[ip][iq]) > tresh)
    				{
    					h = d[iq]-d[ip];
    					if (fabs(h)+g == fabs(h))
    					{
    						t = (a[ip][iq])/h;
    					}
    					else
    					{
    						theta = 0.5*h/(a[ip][iq]);
    						t     = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
    						if (theta < 0.0)
    						{
    							t = -t;
    						}
    					}
    					c         = 1.0/sqrt(1+t*t);
    					s         = t*c;
    					tau       = s/(1.0+c);
    					h         = t*a[ip][iq];
    					z[ip]    -= h;
    					z[iq]    += h;
    					d[ip]    -= h;
    					d[iq]    += h;
    					a[ip][iq] = 0.0;
    					for (j = 0; j < ip; j++)
    					{
    						do_rotate(a, j, ip, j, iq, tau, s);
    					}
    					for (j = ip+1; j < iq; j++)
    					{
    						do_rotate(a, ip, j, j, iq, tau, s);
    					}
    					for (j = iq+1; j < n; j++)
    					{
    						do_rotate(a, ip, j, iq, j, tau, s);
    					}
    					for (j = 0; j < n; j++)
    					{
    						do_rotate(v, j, ip, j, iq, tau, s);
    					}
    					++(*nrot);
    				}
    			}
    		}
    		for (ip = 0; ip < n; ip++)
    		{
    			b[ip] +=  z[ip];
    			d[ip]  =  b[ip];
    			z[ip]  =  0.0;
    		}
    	}
    	throw string("Error: Too many iterations in routine JACOBI\n");
    }catch(const string & s){
    	cout << s << endl;
    	exit(1);
    }

}

#endif /* LIBTRAJ_JACOBI_H_ */
