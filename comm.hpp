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
#define mpicell 3000 // Number of exchange cells for MPI
#define ru 8.3141
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#ifndef SCARF_DECLSPEC
#define SCARF_DECLSPEC extern
#endif

SCARF_DECLSPEC int neqn, nproc;
SCARF_DECLSPEC thread_local int MyID;
SCARF_DECLSPEC int I0; // Number of cells in x direction for each thread, including ghost cells
SCARF_DECLSPEC int J0; // Number of cells in y direction, including ghost cells
SCARF_DECLSPEC int neqv; // Solution variables without chemical terms
SCARF_DECLSPEC double dxc; // Non-dimensional length of each cell in xi direction
SCARF_DECLSPEC double dyc; // Non-dimensional length of each cell in et direction

SCARF_DECLSPEC struct strct_configInt
{
	int newrun, nonDi, useDt,
		ifilm, gasModel, visModel, reacModel, transModel,
	    x_sh, nmix, nspec, nreac, thermo_base;
	int timeOrder; // Time order of the Runge-Kutta method (?)
	int ni; // Number of cells in x direction for each thread
	int nj; // Number of cells in y direction
	int Ng; // Number of ghost cells
	int iStep0; // The initial time step of simulation
	int nStep; // Total time step of simulation
	int nRamp; // Time steps for the program march from CFL0 to CFL1 (or dt0 to dt1)
	int Samples; // Intervals for the program to output results
} config1;

SCARF_DECLSPEC struct strct_configDouble
{
	double x0, Sc0, Re0, condRef, diffRef;
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
	double p1, T1, u1, v1, p2, T2, u2, v2; // Initial condition
} config2;

/*- reference state -*/
SCARF_DECLSPEC double LRef, uRef, cvRef, rhoRef, tRef, temRef, preRef,
       rgasRef, muRef, condRef, diffRef, Upsilon;

/*- MPI related -*/
SCARF_DECLSPEC struct strct_gcelltype
{
	double rho, u, v, e, p, t, mu, kt, ga, qs[6], di[6];
}mpiSend_ql[mpicell], mpiSend_qr[mpicell], mpiRecv_ql[mpicell], mpiRecv_qr[mpicell];

/*- Origin flow variables -*/
SCARF_DECLSPEC struct strct_U
{
	double **qs, **di;
	double ** q; // (rho, u, v, e)
	double * pre; // Pressure
	double * tem; // Temperature
	double * mu; // Coefficient of viscosity
	double * kt; // Coefficient of heat conductivity
	double * cv; // Constant pressure heat capacity
	double * rgas; // Gas constant (R/molWeight)
	double * gam; // Heat capacity ratio (C_V/C_p)
} U, Ug;

/*- viscous related variables -*/
SCARF_DECLSPEC struct strct_Uv
{
	// origin variables derivatives
	double *u_xi, *u_et, *v_xi, *v_et, *T_xi, *T_et, **qs_xi, **qs_et;

	// geometry derivatives coefficients
	double *fu1, *fu2, *fu3, *fv1, *fv2, *fv3, *fuv, *fe1, *fe2;
	double *gu1, *gu2, *gu3, *gv1, *gv2, *gv3, *guv, *ge1, *ge2;
} Uv;

/*- flux related variables -*/
SCARF_DECLSPEC thread_local struct strct_flux
{
	double *xix, *xiy, *etx, *ety, *yas, *rho,
		*du, *dv, *dt, **dqs,
		*u, *v,
		*p, // Pressure
		*e,
		*t, // Temperature
		*gam, // Heat capacity ratio (C_V/C_p)
		*mu, *kt, **Ds, **qs, **flux;
} U1d;

/*- initial condition variables -*/
SCARF_DECLSPEC struct strct_ic
{
	double u, v, t, p, ys[6];
} inc[2];

/* Geometry variables (coordinations and derivatives) */
SCARF_DECLSPEC struct strct_metric
{
	double* x, * y; // Physical coordination of the cell (x, y)
	double * xi, * et; // Non-dimensional coordination of the cell (xi, et)
	double * x_xi, * x_et, * y_xi, * y_et; // The derivatives of the physical coordination (x, y) over the non-dimensional coordination (xi, et)
	double * yaks; // Jacobian (det(partial(x, y)/partial(xi, et)))
} mesh;

/*- others -*/
SCARF_DECLSPEC thread_local double **qso, ***dsdq;
SCARF_DECLSPEC thread_local double** qo; // (rho, u, v, e) value from the last step
SCARF_DECLSPEC thread_local double** rhs; // The increment value of (rho, u, v, e)

SCARF_DECLSPEC thread_local double **dmat1, **sour, **dtdq;

#endif /* COMM_H_ */