/*
 * comm.h
 *
 *  Created on: Jan 07, 2014
 *
 */

#ifndef COMM_H_
#define COMM_H_

#include<stdio.h>
#include<stdlib.h>

/*- constant & func. -*/
#define maxeqn 9
#define ru 8.3141
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

inline int I0; // Number of cells in x direction for each thread, including ghost cells
inline int J0; // Number of cells in y direction, including ghost cells
inline int neqv; // Solution variables without chemical terms
inline double dxc; // Non-dimensional length of each cell in xi direction
inline double dyc; // Non-dimensional length of each cell in et direction

inline struct strct_configInt
{
	int timeOrder; // Time order of the Runge-Kutta method (?)
	int ni; // Number of cells in x direction for each thread
	int nj; // Number of cells in y direction
	int Ng; // Number of ghost cells
	int iStep0; // The initial time step of simulation
	int nStep; // Total time step of simulation
	int nRamp; // Time steps for the program march from CFL0 to CFL1 (or dt0 to dt1)
	int Samples; // Intervals for the program to output results
}config1;

inline struct strct_configDouble
{
	double molWeight; // Average mole weight of the gas
	double gam0; // Heat capacity ratio (C_V/C_p)
	double Pr0; // Prandtl number (C_p*mu/lambda)
	double t0; // The initial time of simulation
	double Lx; // Length of the rectangle in x direction
	double Ly; // Length of the rectangle in y direction
	double dt0; // The minimum time step
	double dt1; // The maximum time step
	double CFL0; // The minimum CFL number
	double CFL1; // The maximum CFL number
	double muRef, suthC1, suthC2; // The constants for Sutherland's law
	double MaRef; // Reference Mach number, which only used for the non-dimension of velocity and have no relation to the actual Mach number
	double temRef; // Reference temperature which only used for non-perfect gas to recover the dimension, while for perfect gas flow not used
	double preRef; // Reference density which only used for non-perfect gas to recover the dimension, while for perfect gas flow not used
	double p1, T1, u1, v1; // Initial condition
}config2;

/* Origin flow variables */
inline struct strct_U
{
	double ** q; // (rho, u, v, e)
	double * pre; // Pressure
	double * tem; // Temperature
	double * mu; // Coefficient of viscosity
	double * kt; // Coefficient of heat conductivity
	double * cv; // Constant pressure heat capacity
	double * rgas; // Gas constant (R/molWeight)
	double * gam; // Heat capacity ratio (C_V/C_p)
}U;

/*- viscous related variables -*/
inline struct strct_Uv
{
	// origin variables derivatives
	double* u_xi, * u_et, * v_xi, * v_et, * T_xi, * T_et;

	// geometry derivatives coefficients
	double* fu1, * fu2, * fu3, * fv1, * fv2, * fv3, * fuv, * fe1, * fe2,
		* gu1, * gu2, * gu3, * gv1, * gv2, * gv3, * guv, * ge1, * ge2;
}Uv;

/*- flux related variables -*/
inline thread_local struct strct_flux
{
	double* xix, * xiy, * etx, * ety, * yas, * rho,
		* du, * dv, * dt,
		* u, * v,
		*p, // Pressure
		*e, // Energy
		*t, // Temperature
		*gam, // Heat capacity ratio (C_V/C_p)
		* mu, * kt, ** flux;
}U1d;

/* Geometry variables (coordinations and derivatives) */
inline struct strct_metric
{
	double* x, * y; // Physical coordination of the cell (x, y)
	double * xi, * et; // Non-dimensional coordination of the cell (xi, et)
	double * x_xi, * x_et, * y_xi, * y_et; // The derivatives of the physical coordination (x, y) over the non-dimensional coordination (xi, et)
	double * yaks; // Jacobian (det(partial(x, y)/partial(xi, et)))
}mesh;

inline double** qo; // (rho, u, v, e) value from the last step
inline double** rhs; // The increment value of (rho, u, v, e)

void readjob();
void setjob();
void jobbody();
void endjob();

#endif /* COMM_H_ */