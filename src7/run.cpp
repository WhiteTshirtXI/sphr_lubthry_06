/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <math.h>
#include <string>
#include "../include/fd.h"
#include "../include/rw.h"

using namespace std;

int main(){
	int i, j, m, n;
	char fn[1024];
	double params[50] = {0.0};

	// read parameters
	sprintf(fn, "./params.in" );
	readInput(fn, params);

	double Ca    = params[0];
	double Bo    = params[1];
	double Ma    = params[2];
	double tstop = params[3];
	double r1    = params[4];
	double t1    = params[5];
	double dr    = params[6];
	double dt    = params[7];
	double dtrec = params[8];
	
	int J  = r1/dr;
	int N  = t1/dt + 1;
	int M = N*dt/dtrec;
	
	int J1 = J+1;
	int M1 = M+1;
	int M1J1 = M1*J1;
	int M1J  = M1*J;
	
	// allocate memory and initialize to zero
	double *R  = (double*) calloc(M1J1, sizeof(double));
	double *T  = (double*) calloc(M1J1, sizeof(double));
	double *H  = (double*) calloc(M1J1, sizeof(double));
	double *G  = (double*) calloc(M1J1, sizeof(double));
	double *F  = (double*) calloc(M1J1, sizeof(double));
	double *P  = (double*) calloc(M1J1, sizeof(double));
	double *Q  = (double*) calloc(M1J1, sizeof(double));
	double *Vs = (double*) calloc(M1J1, sizeof(double));
	
	// time evolution of h, g, and f
	fdevol(J, N, M, dr, dt, params, H, G, F);
	
	// back-calculate p, q, and vs
	fdaux (J, N, M, dr, dt, params, H, G, F, P, Q, Vs);

	// space and time domains
	fdgrid(J, N, M, dr, dt, R, T);

	// write to file
	sprintf(fn, "r" ); write(J, M, dr, dt, params, R , fn);
	sprintf(fn, "t" ); write(J, M, dr, dt, params, T , fn);
	sprintf(fn, "h" ); write(J, M, dr, dt, params, H , fn);
	sprintf(fn, "g" ); write(J, M, dr, dt, params, G , fn);
	sprintf(fn, "f" ); write(J, M, dr, dt, params, F , fn);
	sprintf(fn, "p" ); write(J, M, dr, dt, params, P , fn);
	sprintf(fn, "q" ); write(J, M, dr, dt, params, Q , fn);
	sprintf(fn, "vs"); write(J, M, dr, dt, params, Vs, fn);

	// free memory
	free(R );
	free(T );
	free(H );
	free(G );
	free(F );
	free(P );
	free(Q );
	free(Vs);
	
	return(0);
}
