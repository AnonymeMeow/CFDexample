/*
 * flux.c
 *
 *  Created on:    Jan 14, 2014
 *  last modified: Sep 10, 2014
 */

#include<string.h>
#include<math.h>
#include"comm.h"
#include"chemdata.h"

/*--------------------------------------------------------------
 * Calculate Fluxes of perfect-gas flow
 * -------------------------------------------------------------*/
void flux(double **rhs)
{
	void adGhost();
	void fluxF(double **rhs);
	void fluxG(double **rhs);
	void vfluxF(double **rhs);
	void vfluxG(double **rhs);

	adGhost();
	fluxF(rhs);
	fluxG(rhs);
	if(config1.visModel != 0)
	{
		vfluxF(rhs);
		vfluxG(rhs);
	}

}

/*--------------------------------------------------------------
 * Calculate Fluxes of real-gas flow
 * -------------------------------------------------------------*/
void fluxchem(double **rhs)
{
	void adGhost();
	void fluxchemF(double **rhs);
	void fluxchemG(double **rhs);
	void vfluxchemF(double **rhs);
	void vfluxchemG(double **rhs);

	adGhost();

	fluxchemF(rhs);
	fluxchemG(rhs);
	if(config1.visModel != 0)
	{
		vfluxchemF(rhs);
		vfluxchemG(rhs);
	}

}

/*---------------------------------------------------
 * Calculate the inviscid flux in x direction
 * ------------------------------------------------*/
void fluxF(double **rhs)
{
	int    i, ir, ii, il, j, jj, jr, ic, s, k, ik, iv;
	double ph[maxeqn], phi_N[maxeqn], dsm[5], dsp[5], Fplus[6][maxeqn], UU[maxeqn],
		   Fminus[6][maxeqn], dFplus[5][maxeqn], dFminus[5][maxeqn], LF[maxeqn],
		   qave[maxeqn], f06[maxeqn], le[maxeqn][maxeqn], re[maxeqn][maxeqn],
		   phip, phim, sum1, sum2, sum3, c, alf, maxLamda, pave, tave, gave,
		   te, temp, xix, xiy;

	double lf[6] = {0., -1./12., 7./12., 7./12., -1./12., 0.};

	void boundX();
	double phin(double fa, double fb, double fc, double fd);
	void allocateFlux(int nlen, struct strct_flux *f);
	void freeFlux(int nlen, struct strct_flux *f);
	void getEigenvector(double qave[], double p, double t, double ga, double kx, 
		                double ky, double (*le)[maxeqn], double (*re)[maxeqn]);

	il = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;

    boundX();

    allocateFlux(I0, &U1d);

	for(j=config1.Ng; j<jr; j++)
	{
		for(i=0; i<I0; i++)
		{
		/*---- convert to 1D-array ----*/
			ic = i*J0 + j;

			U1d.yas[i] =  mesh.yaks[ic];
			U1d.xix[i] =  mesh.y_et[ic]/U1d.yas[i];
			U1d.xiy[i] = -mesh.x_et[ic]/U1d.yas[i];
			U1d.rho[i] =  Ug.q[ic][0];
			U1d.u[i]   =  Ug.q[ic][1];
			U1d.v[i]   =  Ug.q[ic][2];
			U1d.e[i]   =  Ug.q[ic][3];
			U1d.p[i]   =  Ug.pre[ic];
			U1d.t[i]   =  Ug.tem[ic];
			U1d.gam[i] =  Ug.gam[ic];
		}

		for(i=il; i<ir; i++) // loop for all the i faces
		{
		/*---- 1. obtain the averaged value ----*/
			xix     = 0.5*(U1d.xix[i] + U1d.xix[i+1]);
			xiy     = 0.5*(U1d.xiy[i] + U1d.xiy[i+1]);
			pave    = 0.5*(U1d.p[i]   + U1d.p[i+1]);
			tave    = 0.5*(U1d.t[i]   + U1d.t[i+1]);
			gave    = 0.5*(U1d.gam[i] + U1d.gam[i+1]);
			qave[0] = 0.5*(U1d.rho[i] + U1d.rho[i+1]);
			qave[1] = 0.5*(U1d.u[i] + U1d.u[i+1]);
			qave[2] = 0.5*(U1d.v[i] + U1d.v[i+1]);
			qave[3] = 0.5*(U1d.e[i] + U1d.e[i+1]);

		/*---- 2. calculate eigenvector ----*/
			getEigenvector(qave,pave,tave,gave,xix,xiy,le,re);

		/*---- 3. calculate split Flux and central term from the six stencils cells ----*/
			maxLamda = 0.;
			for(k=0; k<=5; k++)
			{
				ik  = i-2 + k;
				te  = U1d.u[ik]*U1d.xix[ik] + U1d.v[ik]*U1d.xiy[ik];
				alf = sqrt(U1d.xix[ik]*U1d.xix[ik] + U1d.xiy[ik]*U1d.xiy[ik]);
				c   = sqrt(U1d.gam[ik]*U1d.p[ik]/U1d.rho[ik]*Upsilon);
				temp = fabs(te) + c*alf;
				if(temp > maxLamda)
				  	maxLamda = temp;
			}

			for(iv=0; iv<neqn; iv++)
				LF[iv] = 0.;
			for(k=0; k<=5; k++)
			{
				ik = i-2 + k;
				/*  i-2, i-1, i, i+1, i+2, i+3.
				 * The corresponding cells begin with 2-2 = 0,
				 * end with (NI+2)+3 = NI+5 */

				te     = U1d.u[ik]*U1d.xix[ik] + U1d.v[ik]*U1d.xiy[ik];
				f06[0] = U1d.yas[ik] * U1d.rho[ik]*te;
				f06[1] = U1d.yas[ik] * (U1d.rho[ik]*U1d.u[ik]*te + U1d.xix[ik]*U1d.p[ik]*Upsilon);
				f06[2] = U1d.yas[ik] * (U1d.rho[ik]*U1d.v[ik]*te + U1d.xiy[ik]*U1d.p[ik]*Upsilon);
				f06[3] = U1d.yas[ik] * (U1d.rho[ik]*U1d.e[ik] + U1d.p[ik]*Upsilon)*te;

				UU[0] = U1d.rho[ik];
				UU[1] = U1d.rho[ik]*U1d.u[ik];
				UU[2] = U1d.rho[ik]*U1d.v[ik];
				UU[3] = U1d.rho[ik]*U1d.e[ik];

				for(iv=0; iv<neqn; iv++)
				{
					Fplus[k][iv]  = 0.5*(f06[iv] + maxLamda*UU[iv]*U1d.yas[ik]);
					Fminus[k][iv] = 0.5*(f06[iv] - maxLamda*UU[iv]*U1d.yas[ik]);

					LF[iv] = LF[iv] + lf[k]*f06[iv];
				}
			}

			/*---- 4. calculate delta flux ----*/
			for(k=0; k<=4; k++)
				for(iv=0; iv<neqn; iv++)
				{
					dFplus[k][iv]  = Fplus[k+1][iv]  - Fplus[k][iv];
					dFminus[k][iv] = Fminus[k+1][iv] - Fminus[k][iv];
				}

			/*---- 5. Approximate the fluxes in local characteristic field ----*/
			for(s=0; s<neqn; s++)
			{
				for(k=0; k<=4; k++)
				{
					sum1 = 0.;
					sum2 = 0.;
					for (iv = 0; iv<neqn; iv++)
					{
						sum1 = sum1 + dFplus[k][iv]*le[s][iv];
						sum2 = sum2 + dFminus[k][iv]*le[s][iv];
					}
					dsp[k] = sum1;
					dsm[k] = sum2;
				}

				phip = phin(dsp[0], dsp[1], dsp[2], dsp[3]);
				phim = phin(dsm[4], dsm[3], dsm[2], dsm[1]);

				ph[s] = -phip + phim;
			}
			/*---- 6. Project back to the component space, get the upwind flux ----*/
			for(s=0; s<neqn; s++)
			{
				sum3 = 0.;
				for(iv=0; iv<neqn; iv++)
					sum3 = sum3 + ph[iv]*re[s][iv];

				phi_N[s] = sum3;
			}
			/*---- 7. get the final flux ----*/
			for(iv=0; iv<neqn; iv++)
				U1d.flux[i][iv] = LF[iv] + phi_N[iv];
		}

		jj = j-config1.Ng;
		for(i=config1.Ng; i<ir; i++)
		{
			ii = i - config1.Ng;
			ic = ii*config1.nj + jj;
			for(iv=0; iv<neqn; iv++)
				rhs[ic][iv] = - (U1d.flux[i][iv] - U1d.flux[i-1][iv])/dxc;
		}
	}
	freeFlux(I0, &U1d);
}


/*---------------------------------------------------
 * Calculate the inviscid flux in y direction
 * ------------------------------------------------*/
void fluxG(double **rhs)
{
	int    i, ir, ii, j, jr, jj, jl, ic, s, k, jk, iv;
	double ph[maxeqn], phi_N[maxeqn], dsm[5], dsp[5], Fplus[6][maxeqn], UU[maxeqn],
		   Fminus[6][maxeqn], dFplus[5][maxeqn] , dFminus[5][maxeqn], LF[maxeqn],
		   qave[maxeqn], f06[maxeqn], le[maxeqn][maxeqn], re[maxeqn][maxeqn],
		   phip, phim, sum1, sum2, sum3, c, alf, maxLamda, pave, tave, gave,
		   te, temp, etx, ety;

	double lf[6] = {0., -1./12., 7./12., 7./12., -1./12., 0.};

	double phin(double fa, double fb, double fc, double fd);
	void boundY();
	void allocateFlux(int nlen, struct strct_flux *f);
	void freeFlux(int nlen, struct strct_flux *f);
	void getEigenvector(double qave[], double p, double t, double ga, double kx, 
		                double ky, double (*le)[maxeqn], double (*re)[maxeqn]);


	ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;
	jl = config1.Ng - 1;

	boundY();

	allocateFlux(J0,&U1d);

	for(i=config1.Ng; i<ir; i++)
	{
		for(j=0; j<J0; j++)
		{
			/*---- convert to 1D-array ----*/
			ic = i*J0 + j;

			U1d.yas[j] =  mesh.yaks[ic];
			U1d.etx[j] = -mesh.y_xi[ic]/U1d.yas[j];
			U1d.ety[j] =  mesh.x_xi[ic]/U1d.yas[j];
			U1d.rho[j] =  Ug.q[ic][0];
			U1d.u[j]   =  Ug.q[ic][1];
			U1d.v[j]   =  Ug.q[ic][2];
			U1d.e[j]   =  Ug.q[ic][3];
			U1d.p[j]   =  Ug.pre[ic];
			U1d.t[j]   =  Ug.tem[ic];
			U1d.gam[j] =  Ug.gam[ic];
		}

		for(j=jl; j<jr; j++) // loop for j faces get the value at lines of fluid without ghost cell
		{
			etx     = 0.5*(U1d.etx[j]  + U1d.etx[j+1]);
			ety     = 0.5*(U1d.ety[j]  + U1d.ety[j+1]);
			pave    = 0.5*(U1d.p[j]    + U1d.p[j+1]);
			tave    = 0.5*(U1d.t[j]    + U1d.t[j+1]);
			gave    = 0.5*(U1d.gam[j]  + U1d.gam[j+1]);
			qave[0] = 0.5*(U1d.rho[j]  + U1d.rho[j+1]);
			qave[1] = 0.5*(U1d.u[j]    + U1d.u[j+1]);
			qave[2] = 0.5*(U1d.v[j]    + U1d.v[j+1]);
			qave[3] = 0.5*(U1d.e[j]    + U1d.e[j+1]);

			getEigenvector(qave,pave,tave,gave,etx,ety,le,re);

			maxLamda = 0.;
			for(k=0; k<=5; k++)
			{
				jk  = j-2 + k;//jk 0-5
				te  = U1d.u[jk]*U1d.etx[jk] + U1d.v[jk]*U1d.ety[jk];//v
				alf = sqrt(U1d.etx[jk]*U1d.etx[jk] + U1d.ety[jk]*U1d.ety[jk]);//1
				c   = sqrt(U1d.gam[jk]*U1d.p[jk]/(U1d.rho[jk])*Upsilon);//c
				temp = fabs(te) + c*alf;//v + c
				if(temp > maxLamda)
				  	maxLamda = temp;
			}

			for(iv=0; iv<neqn; iv++)//neqn=4 number of solution without chemical terms
				LF[iv] = 0.;
			for(k=0; k<=5; k++)//j = 0
			{
				jk = j-2 + k;//jk 0-5

				te     = U1d.u[jk]*U1d.etx[jk] + U1d.v[jk]*U1d.ety[jk];//v
				f06[0] = U1d.yas[jk] * U1d.rho[jk]*te;//rho*v
				f06[1] = U1d.yas[jk] * (U1d.rho[jk]*U1d.u[jk]*te + U1d.etx[jk]*U1d.p[jk]*Upsilon);//rho*uv
				f06[2] = U1d.yas[jk] * (U1d.rho[jk]*U1d.v[jk]*te + U1d.ety[jk]*U1d.p[jk]*Upsilon);//rho*vv
				f06[3] = U1d.yas[jk] * (U1d.rho[jk]*U1d.e[jk] + U1d.p[jk]*Upsilon)*te;//rho*e + p

				UU[0] = U1d.rho[jk];
				UU[1] = U1d.rho[jk]*U1d.u[jk];
				UU[2] = U1d.rho[jk]*U1d.v[jk];
				UU[3] = U1d.rho[jk]*U1d.e[jk];

				for(iv=0; iv<neqn; iv++)//0-3
				{
					Fplus[k][iv]  = 0.5*(f06[iv] + maxLamda*UU[iv]*U1d.yas[jk]);//rho*v+lamda*rho change
					Fminus[k][iv] = 0.5*(f06[iv] - maxLamda*UU[iv]*U1d.yas[jk]);//rho*v-lamda*rho change

					LF[iv] = LF[iv] + lf[k]*f06[iv];//1/12*
				}
			}

			for(k=0; k<=4; k++)//0-4
				for(iv = 0; iv<neqn; iv++)
				{
					dFplus[k][iv]  = Fplus[k+1][iv]  - Fplus[k][iv];//1-5  -  0-4  deltaflux
					dFminus[k][iv] = Fminus[k+1][iv] - Fminus[k][iv];
				}

			for(s=0; s<neqn; s++)//0-3
			{
				for(k=0; k<=4; k++)//0-4
				{
					sum1 = 0.;
					sum2 = 0.;
					for (iv=0; iv<neqn; iv++)
					{
						sum1 = sum1 + dFplus[k][iv]*le[s][iv];//delta*left engiv
						sum2 = sum2 + dFminus[k][iv]*le[s][iv];
					}
					dsp[k] = sum1;
					dsm[k] = sum2;
				}

				phip = phin(dsp[0], dsp[1], dsp[2], dsp[3]);
				phim = phin(dsm[4], dsm[3], dsm[2], dsm[1]);

				ph[s] = -phip + phim;//upwind part phi
			}

			for(s=0; s<neqn; s++)
			{
				sum3 = 0.;
				for(iv = 0; iv<neqn; iv++)
					sum3 = sum3 + ph[iv]*re[s][iv];//phi*right engiv
				phi_N[s] = sum3;
			}

				for(iv=0; iv<neqn; iv++)
					U1d.flux[j][iv] = LF[iv] + phi_N[iv];//final part fi+1/2

		}

		ii = i - config1.Ng;//i = 3 -- N+2    ii  0 - N-1
		for(j=config1.Ng; j<jr; j++)// j = 3 - N+2
		{
			jj = j - config1.Ng;//0 - N-1
			ic = ii*config1.nj + jj;//point in fulid part
			for(iv=0; iv<neqn; iv++)
				rhs[ic][iv] = rhs[ic][iv] -(U1d.flux[j][iv] - U1d.flux[j-1][iv])/dyc;

			// the rhs at j=config1.Ng is actually not used
		}

	}
	freeFlux(J0, &U1d);
}

/*---------------------------------------------------
 * Calculate inviscid flux of real-gas flow in x direction
 * ------------------------------------------------*/
void fluxchemF(double **rhs)
{
	int    i, ir, ii, il, j, jj, jr, ic, s, k, ik, ns, iv, m;
	double ph[maxeqn], phi_N[maxeqn], dsm[5], dsp[5], Fplus[6][maxeqn], Fminus[6][maxeqn],
	       UU[maxeqn], dFplus[5][maxeqn], dFminus[5][maxeqn], LF[maxeqn], f06[maxeqn],
	       le[maxeqn][maxeqn], re[maxeqn][maxeqn], qave[maxeqn], qsave[maxeqn],phip, phim,
	       sum1, sum2, c, alf, maxLamda, pave, tave, gave, te, temp, xix, xiy;

	double lf[6] = {0., -1./12., 7./12., 7./12., -1./12., 0.};

	void boundX();
	void allocateFlux(int nlen, struct strct_flux *f);
	void freeFlux(int nlen, struct strct_flux *f);
	double phin(double fa, double fb, double fc, double fd);
	void getEigenvectorChem(double qave[], double qsave[], double p, double t, double ga, 
		                    double kx, double ky, double (*le)[maxeqn], double (*re)[maxeqn]);


    ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;

    boundX();

    m = config1.nspec;
    allocateFlux(I0, &U1d);

	for(j=config1.Ng; j<jr; j++)
	{
		for(i=0; i<I0; i++)
		{
		/*---- convert to 1D-array ----*/
			ic = i*J0 + j;

			U1d.yas[i] =  mesh.yaks[ic];
			U1d.xix[i] =  mesh.y_et[ic]/mesh.yaks[ic];
			U1d.xiy[i] = -mesh.x_et[ic]/mesh.yaks[ic];
			U1d.rho[i] =  Ug.q[ic][0];
			U1d.u[i]   =  Ug.q[ic][1];
			U1d.v[i]   =  Ug.q[ic][2];
			U1d.e[i]   =  Ug.q[ic][3];
			U1d.p[i]   =  Ug.pre[ic];
			U1d.t[i]   =  Ug.tem[ic];
			U1d.gam[i] =  Ug.gam[ic];
			for(ns=0; ns<m; ns++)
				U1d.qs[i][ns] = Ug.qs[ic][ns];
		}

		il = config1.Ng - 1;
		for(i=il; i<ir; i++)
		{
		/*---- 1. obtain the averaged value ----*/
			xix     = 0.5*(U1d.xix[i] + U1d.xix[i+1]);
			xiy     = 0.5*(U1d.xiy[i] + U1d.xiy[i+1]);
			pave    = 0.5*(U1d.p[i]   + U1d.p[i+1]);
			tave    = 0.5*(U1d.t[i]   + U1d.t[i+1]);
			gave    = 0.5*(U1d.gam[i] + U1d.gam[i+1]);
			qave[0] = 0.5*(U1d.rho[i] + U1d.rho[i+1]);
			qave[1] = 0.5*(U1d.u[i] + U1d.u[i+1]);
			qave[2] = 0.5*(U1d.v[i] + U1d.v[i+1]);
			qave[3] = 0.5*(U1d.e[i] + U1d.e[i+1]);
			for(ns=0; ns<m; ns++)
				qsave[ns] = 0.5*(U1d.qs[i][ns] + U1d.qs[i+1][ns]);

		/*---- 2. calculate eigenvector ----*/
			getEigenvectorChem(qave, qsave, pave, tave, gave, xix, xiy, le, re);

		/*---- 3. calculate split Flux and central term from the six stencils cells ----*/
			maxLamda = 0.;
			for(k=0; k<=5; k++)
			{
				ik  = i-2 + k;
				alf = sqrt(U1d.xix[ik]*U1d.xix[ik] + U1d.xiy[ik]*U1d.xiy[ik]);
				te  = U1d.u[ik]*xix + U1d.v[ik]*xiy;
				c   = sqrt(U1d.gam[ik]*U1d.p[ik]/U1d.rho[ik]*Upsilon);
				temp = fabs(te) + c*alf;
				if(temp > maxLamda)
				  	maxLamda = temp;
			}

			for(iv=0; iv<neqn; iv++)
				LF[iv] = 0.;
			for(k=0; k<=5; k++)
			{
				ik = i-2 + k;

				te  = U1d.u[ik]*U1d.xix[ik] + U1d.v[ik]*U1d.xiy[ik];
				for(ns=0; ns<m; ns++)
				f06[ns] = U1d.yas[ik]*(U1d.rho[ik]*U1d.qs[ik][ns]*te);
				f06[m]   = U1d.yas[ik]*(U1d.rho[ik]*U1d.u[ik]*te + U1d.xix[ik]*U1d.p[ik]*Upsilon);
				f06[m+1] = U1d.yas[ik]*(U1d.rho[ik]*U1d.v[ik]*te + U1d.xiy[ik]*U1d.p[ik]*Upsilon);
				f06[m+2] = U1d.yas[ik]*(U1d.rho[ik]*U1d.e[ik] + U1d.p[ik]*Upsilon)*te;

				for(ns=0; ns<m; ns++)
				UU[ns] = U1d.rho[ik]*U1d.qs[ik][ns];
				UU[m]   = U1d.rho[ik]*U1d.u[ik];
				UU[m+1] = U1d.rho[ik]*U1d.v[ik];
				UU[m+2] = U1d.rho[ik]*U1d.e[ik];

				for(iv=0; iv<neqn; iv++)
				{
					Fplus[k][iv]  = 0.5*(f06[iv] + maxLamda*UU[iv]*U1d.yas[ik]);
					Fminus[k][iv] = 0.5*(f06[iv] - maxLamda*UU[iv]*U1d.yas[ik]);

					LF[iv] = LF[iv] + lf[k]*f06[iv];
				}
			}

			/*---- 4. calculate delta flux ----*/
			for(k=0; k<=4; k++)
				for(iv=0; iv<neqn; iv++)
				{
					dFplus[k][iv]  = Fplus[k+1][iv]  - Fplus[k][iv];
					dFminus[k][iv] = Fminus[k+1][iv] - Fminus[k][iv];
				}

			/*---- 5. Approximate the fluxes in local characteristic field ----*/
			for(s=0; s<neqn; s++)
			{
				for(k=0; k<=4; k++)
				{
					sum1 = 0.;
					sum2 = 0.;
					for (iv = 0; iv<neqn; iv++)
					{
						sum1 = sum1 + dFplus[k][iv]*le[s][iv];
						sum2 = sum2 + dFminus[k][iv]*le[s][iv];
					}
					dsp[k] = sum1;
					dsm[k] = sum2;
				}

				phip = phin(dsp[0], dsp[1], dsp[2], dsp[3]);
				phim = phin(dsm[4], dsm[3], dsm[2], dsm[1]);

				ph[s] = -phip + phim;
			}
			/*---- 6. Project back to the component space, get the upwind fluxes ----*/
			for(s=0; s<neqn; s++)
			{
				sum1 = 0.;
				for(iv = 0; iv<neqn; iv++)
					sum1 = sum1 + ph[iv]*re[s][iv];

				phi_N[s] = sum1;
			}

			/*---- 7. get the final flux ----*/

			for(iv=0; iv<neqn; iv++)
				U1d.flux[i][iv] = LF[iv] + phi_N[iv];
		}

		/*------ Near wall boundary, reduced to up-wind scheme ------*/
		if(MyID == NMAXproc)
		{
			for(k=0; k<config1.Ng; k++)
			{
				ik = config1.ni + k;

				te  = U1d.u[ik]*U1d.xix[ik] + U1d.v[ik]*U1d.xiy[ik];
				alf = sqrt(U1d.xix[ik]*U1d.xix[ik] + U1d.xiy[ik]*U1d.xiy[ik]);
				c   = sqrt(U1d.gam[ik]*U1d.p[ik]/U1d.rho[ik]*Upsilon);
				temp = fabs(te) + c*alf;
				if(temp > maxLamda)
				  	maxLamda = temp;
			}
			for(k=0; k<config1.Ng; k++)
			{
				ik = config1.ni + k;

				te  = U1d.u[ik]*U1d.xix[ik] + U1d.v[ik]*U1d.xiy[ik];
				for(ns=0; ns<m; ns++)
					f06[ns] = U1d.yas[ik]*(U1d.rho[ik]*U1d.qs[ik][ns]*te);
				f06[m]   = U1d.yas[ik]*(U1d.rho[ik]*U1d.u[ik]*te + U1d.xix[ik]*U1d.p[ik]*Upsilon);
				f06[m+1] = U1d.yas[ik]*(U1d.rho[ik]*U1d.v[ik]*te + U1d.xiy[ik]*U1d.p[ik]*Upsilon);
				f06[m+2] = U1d.yas[ik]*(U1d.rho[ik]*U1d.e[ik] + U1d.p[ik]*Upsilon)*te;

				for(ns=0; ns<m; ns++)
					UU[ns] = U1d.rho[ik]*U1d.qs[ik][ns];
				UU[m]   = U1d.rho[ik]*U1d.u[ik];
				UU[m+1] = U1d.rho[ik]*U1d.v[ik];
				UU[m+2] = U1d.rho[ik]*U1d.e[ik];

				for(iv=0; iv<neqn; iv++)
				{
					Fplus[k][iv]  = 0.5*(f06[iv] + maxLamda*UU[iv]*U1d.yas[ik]);
					Fminus[k][iv] = 0.5*(f06[iv] - maxLamda*UU[iv]*U1d.yas[ik]);
				}
			}
			for(k=0; k<(config1.Ng-1); k++)
			{
				ii = config1.ni + k;
				for(iv=0; iv<neqn; iv++)
					U1d.flux[ii][iv] = Fplus[k][iv] + Fminus[k+1][iv];
			}
		}
		/*------ End wall boundary correction ------*/

		jj = j-config1.Ng;
		for(i=config1.Ng; i<ir; i++)
		{
			ii = i - config1.Ng;
			ic = ii*config1.nj + jj;
			for(iv=0; iv<neqn; iv++)
				rhs[ic][iv] = - (U1d.flux[i][iv] - U1d.flux[i-1][iv])/dxc;

		}
	}
	freeFlux(I0, &U1d);
}


/*---------------------------------------------------
 * Calculate inviscid flux of real-gas flow in x direction
 * ------------------------------------------------*/
void fluxchemG(double **rhs)
{
	int    i, ir, ii, j, jr, jj, jl, ic, s, k, jk, ns, iv, m;
	double ph[maxeqn], phi_N[maxeqn], dsm[5], dsp[5], Fplus[6][maxeqn], Fminus[6][maxeqn],
	       UU[maxeqn], dFplus[5][maxeqn], dFminus[5][maxeqn], LF[maxeqn], f06[maxeqn],
	       le[maxeqn][maxeqn], re[maxeqn][maxeqn], qave[maxeqn], qsave[maxeqn], phip, phim,
	       sum1, sum2, sum3, c, alf, maxLamda, pave, tave, gave, te, temp, etx, ety;

	double lf[6] = {0., -1./12., 7./12., 7./12., -1./12., 0.};

	void boundY();
	double phin(double fa, double fb, double fc, double fd);
	void allocateFlux(int nlen, struct strct_flux *f);
	void freeFlux(int nlen, struct strct_flux *f);
	void getEigenvectorChem(double qave[], double qsave[], double p, double t, double ga, 
		                    double kx, double ky, double (*le)[maxeqn], double (*re)[maxeqn]);

	ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;

	boundY();

	m = config1.nspec;
	allocateFlux(J0,&U1d);

	for(i=config1.Ng; i<ir; i++)
	{
		for(j=0; j<J0; j++)
		{
			/*---- convert to 1D-array ----*/
			ic = i*J0 + j;

			U1d.yas[j] =  mesh.yaks[ic];
			U1d.etx[j] = -mesh.y_xi[ic]/U1d.yas[j];
			U1d.ety[j] =  mesh.x_xi[ic]/U1d.yas[j];
			U1d.rho[j] =  Ug.q[ic][0];
			U1d.u[j]   =  Ug.q[ic][1];
			U1d.v[j]   =  Ug.q[ic][2];
			U1d.e[j]   =  Ug.q[ic][3];
			U1d.p[j]   =  Ug.pre[ic];
			U1d.t[j]   =  Ug.tem[ic];
			U1d.gam[j] =  Ug.gam[ic];
			for(ns=0; ns<m; ns++)
				U1d.qs[j][ns] = Ug.qs[ic][ns];
		}

		jl = config1.Ng - 1;
		for(j=jl; j<jr; j++) // loop for j faces
		{
			//change the sequence of mats
			etx  = 0.5*(U1d.etx[j] + U1d.etx[j+1]);
			ety  = 0.5*(U1d.ety[j] + U1d.ety[j+1]);
			pave = 0.5*(U1d.p[j]   + U1d.p[j+1]);
			tave = 0.5*(U1d.t[j]   + U1d.t[j+1]);
			gave = 0.5*(U1d.gam[j] + U1d.gam[j+1]);
			qave[0] = 0.5*(U1d.rho[j]  + U1d.rho[j+1]);
			qave[1] = 0.5*(U1d.u[j]    + U1d.u[j+1]);
			qave[2] = 0.5*(U1d.v[j]    + U1d.v[j+1]);
			qave[3] = 0.5*(U1d.e[j]    + U1d.e[j+1]);
			for(ns=0; ns<m; ns++)
				qsave[ns] = 0.5*(U1d.qs[j][ns] + U1d.qs[j+1][ns]);

			getEigenvectorChem(qave, qsave, pave, tave, gave, etx, ety, le, re);

			maxLamda = 0.;
			for(k=0; k<=5; k++)
			{
				jk  = j-2 + k;
				te  = U1d.u[jk]*U1d.etx[jk] + U1d.v[jk]*U1d.ety[jk];
				alf = sqrt(U1d.etx[jk]*U1d.etx[jk] + U1d.ety[jk]*U1d.ety[jk]);
				c   = sqrt(U1d.gam[jk]*U1d.p[jk]/(U1d.rho[jk])*Upsilon);
				temp = fabs(te) + c*alf;
				if(temp > maxLamda)
				  	maxLamda = temp;
			}
			for(iv=0; iv<neqn; iv++)
				LF[iv] = 0.;

			for(k=0; k<=5; k++)
			{
				jk = j-2 + k;

				te  = U1d.u[jk]*U1d.etx[jk] + U1d.v[jk]*U1d.ety[jk];
				for(ns=0; ns<m; ns++)
					f06[ns] = U1d.yas[jk] * (U1d.rho[jk]*U1d.qs[jk][ns]*te);
				f06[m]   = U1d.yas[jk] * (U1d.rho[jk]*U1d.u[jk]*te + U1d.etx[jk]*U1d.p[jk]*Upsilon);
				f06[m+1] = U1d.yas[jk] * (U1d.rho[jk]*U1d.v[jk]*te + U1d.ety[jk]*U1d.p[jk]*Upsilon);
				f06[m+2] = U1d.yas[jk] * (U1d.rho[jk]*U1d.e[jk] + U1d.p[jk]*Upsilon)*te;

				for(ns=0; ns<m; ns++)
					UU[ns] = U1d.rho[jk]*U1d.qs[jk][ns];
				UU[m]   = U1d.rho[jk]*U1d.u[jk];
				UU[m+1] = U1d.rho[jk]*U1d.v[jk];
				UU[m+2] = U1d.rho[jk]*U1d.e[jk];

				for(iv=0; iv<neqn; iv++)
				{
					Fplus[k][iv]  = 0.5*(f06[iv] + maxLamda*UU[iv]*U1d.yas[jk]);
					Fminus[k][iv] = 0.5*(f06[iv] - maxLamda*UU[iv]*U1d.yas[jk]);

					LF[iv] = LF[iv] + lf[k]*f06[iv];

				}
			}

			for(k=0; k<=4; k++)
				for(iv = 0; iv<neqn; iv++)
				{
					dFplus[k][iv]  = Fplus[k+1][iv]  - Fplus[k][iv];
					dFminus[k][iv] = Fminus[k+1][iv] - Fminus[k][iv];
				}

			for(s=0; s<neqn; s++)
			{
				for(k=0; k<=4; k++)
				{
					sum1 = 0.;
					sum2 = 0.;
					for(iv=0; iv<neqn; iv++)
					{
						sum1 = sum1 + dFplus[k][iv]*le[s][iv];
						sum2 = sum2 + dFminus[k][iv]*le[s][iv];
					}
					dsp[k] = sum1;
					dsm[k] = sum2;
				}

				phip = phin(dsp[0], dsp[1], dsp[2], dsp[3]);
				phim = phin(dsm[4], dsm[3], dsm[2], dsm[1]);

				ph[s] = -phip + phim;
			}

			for(s=0; s<neqn; s++)
			{
				sum3 = 0.;
				for(iv = 0; iv<neqn; iv++)
					sum3 = sum3 + ph[iv]*re[s][iv];

				phi_N[s] = sum3;
			}

			for(iv = 0; iv<neqn; iv++)
				U1d.flux[j][iv] = LF[iv] + phi_N[iv];

		}

		ii = i - config1.Ng;
		for(j=config1.Ng; j<jr; j++)
		{
			jj = j - config1.Ng;
			ic = ii*config1.nj + jj;
			for(iv=0; iv<neqn; iv++)
				rhs[ic][iv] = rhs[ic][iv] -(U1d.flux[j][iv] - U1d.flux[j-1][iv])/dyc;

		}
	}
	freeFlux(J0, &U1d);
}

/*---------------------------------------------------
 * Copy solution vector U to Ug with ghost cells
 * ------------------------------------------------*/
void adGhost()
{
	int i, j, ii, jj, ic, ic1, ns;

	for(i = 0; i<config1.ni; i++)
	{
		ii = i + config1.Ng;
		for(j = 0; j<config1.nj; j++)
		{
			jj = j + config1.Ng;
			ic = i*config1.nj + j;
			ic1 = ii*J0 + jj;

    		Ug.q[ic1][0] =  U.q[ic][0];
    		Ug.q[ic1][1] =  U.q[ic][1];
    		Ug.q[ic1][2] =  U.q[ic][2];
    		Ug.q[ic1][3] =  U.q[ic][3];
			if(config1.gasModel != 0)
				for(ns=0; ns<config1.nspec; ns++)
				{
					Ug.qs[ic1][ns] = U.qs[ic][ns];
					Ug.di[ic1][ns] = U.di[ic][ns];
				}
    		Ug.tem[ic1] =  U.tem[ic];
    		Ug.pre[ic1] =  U.pre[ic];
        	Ug.gam[ic1] =  U.gam[ic];
        	Ug.rgas[ic1]=  U.rgas[ic];
    		Ug.mu[ic1]  =  U.mu[ic];
    		Ug.kt[ic1]  =  U.kt[ic];
		}
	}
}
/*---------------------------------------------------
 * Calculate the Phi function of WENO scheme
 * ------------------------------------------------*/
double phin(double fa, double fb, double fc, double fd)
{
	double is0, is1, is2, epss, al0, al1, al2, om0, om2, fz;

	epss = 1.e-13;

	is0 = 13.*(fa - fb)*(fa - fb) + 3.*(fa - 3.*fb)*(fa - 3.*fb);
	is1 = 13.*(fb - fc)*(fb - fc) + 3.*(fb + fc)*(fb + fc);
	is2 = 13.*(fc - fd)*(fc - fd) + 3.*(3.*fc - fd)*(3.*fc - fd);

	al0 = 1./pow((is0 + epss),2.0);
	al1 = 6./pow((is1 + epss),2.0);
	al2 = 3./pow((is2 + epss),2.0);

	om0 = al0/(al0 + al1 + al2);  om2 = al2/(al0 + al1 + al2);

	fz = 1./3.*om0*(fa - 2.*fb + fc) + 1./6.*(om2 - 0.5)*(fb - 2.*fc + fd);
	return (fz);
}

/*---------------------------------------------------
 * Calculate the Eigen-vectors of perfect gas
 * ------------------------------------------------*/
void getEigenvector(double q[], double p, double t, double ga, double kx, double ky,
					 double (*leftEigenvector)[maxeqn], double (*rightEigenvector)[maxeqn])
{
	double rho, u, v, kx1, ky1, c, alf, beta, term1, r2c,
		   te, H, u2, tiny;

	tiny = 1.e-12;

	rho = q[0];
	u   = q[1];
	v   = q[2];
    beta = ga - 1.;
	u2   = 0.5*(u*u + v*v);
    c    = sqrt(ga*p/(rho + tiny)*Upsilon);
    H    = (q[3] + p/(rho + tiny)*Upsilon);

	term1 = beta/(c*c + tiny);
	r2c = 1./(c*1.414214 + tiny);

	alf = sqrt(kx*kx + ky*ky);
	kx1 = kx/alf;
	ky1 = ky/alf;
	te  = u*kx1 + v*ky1;

	/*Modified from  Volume 2, C.Hirsch, Page 183 */

	leftEigenvector[0][0] =  1. - term1*u2;
	leftEigenvector[0][1] =  u*term1;
    leftEigenvector[0][2] =  v*term1;
    leftEigenvector[0][3] = -term1;
	leftEigenvector[1][0] =  v*kx1 - u*ky1;
	leftEigenvector[1][1] =  ky1;
	leftEigenvector[1][2] = -kx1;
	leftEigenvector[1][3] =  0.;
	leftEigenvector[2][0] =  r2c*(beta*u2 - c*te);
	leftEigenvector[2][1] =  r2c*(kx1*c - beta*u);
	leftEigenvector[2][2] =  r2c*(ky1*c - beta*v);
	leftEigenvector[2][3] =  r2c*beta;
	leftEigenvector[3][0] =  r2c*(beta*u2 + c*te);
	leftEigenvector[3][1] = -r2c*(kx1*c + beta*u);
	leftEigenvector[3][2] = -r2c*(ky1*c + beta*v);
	leftEigenvector[3][3] =  r2c*beta;

	rightEigenvector[0][0] =  1.;
	rightEigenvector[0][1] =  0.;
	rightEigenvector[0][2] =  r2c;
	rightEigenvector[0][3] =  r2c;
	rightEigenvector[1][0] =  u;
	rightEigenvector[1][1] =  ky1;
	rightEigenvector[1][2] =  r2c*(u + c*kx1);
	rightEigenvector[1][3] =  r2c*(u - c*kx1);
	rightEigenvector[2][0] =  v;
	rightEigenvector[2][1] = -kx1;
	rightEigenvector[2][2] =  r2c*(v + c*ky1);
	rightEigenvector[2][3] =  r2c*(v - c*ky1);
	rightEigenvector[3][0] =  u2;
	rightEigenvector[3][1] =  u*ky1 - v*kx1;
	rightEigenvector[3][2] =  r2c*(H + c*te);
	rightEigenvector[3][3] =  r2c*(H - c*te);
}

/*---------------------------------------------------
 * Calculate the Eigen-vectors of real-gas
 * ------------------------------------------------*/
void getEigenvectorChem(double q[], double qsin[], double p, double t, double ga, double kx, double ky,
					  double (*leftEigenvector)[maxeqn], double (*rightEigenvector)[maxeqn])
{
	int i, j, ns, m;
	double rho, u, v, alf, k1x, k1y, c, rc2, beta, te1, H, u2, tiny;
	double Rs[maxspec], es[maxspec], qs[maxspec], pqs[maxspec];
	double getes(int ns, double t);

	tiny = 1.e-20;
	m = config1.nspec;

	rho = q[0];
	u   = q[1];
	v   = q[2];
    beta = ga - 1.;
    u2   = u*u + v*v;
    c    = sqrt(ga*p/(rho + tiny)*Upsilon);
    H    = (q[3] + p/(rho + tiny)*Upsilon);

	alf  = sqrt(kx*kx + ky*ky);
	k1x = kx/alf;
	k1y = ky/alf;
	te1 = u*k1x + v*k1y;

	rc2 = 1./(c*c + tiny);

	/* 1. Calculate pressure derivatives */
    for(ns=0; ns<m; ns++)
    {
		qs[ns] = qsin[ns];
        es[ns] = getes(ns,t);
        Rs[ns] = (ru/specData[ns].wm)/rgasRef;
    }
	for(ns = 0; ns<m; ns++)
		pqs[ns] = 0.5*beta*u2 - beta*es[ns] + Rs[ns]*t*Upsilon;

    /* 2. Calculate left and right Eigenvectors */

	//2.1.a. Left-eigenvectors - upper part
	for(i=0; i<m; i++)
	{
	    for(j=0; j<m; j++)
	    {
	    	if(j == i)
	    		leftEigenvector[i][j] = c*c - qs[i]*pqs[j];
	    	else
	    		leftEigenvector[i][j] =     - qs[i]*pqs[j];
	    }
	    leftEigenvector[i][m]   =  qs[i]*beta*u;
	    leftEigenvector[i][m+1] =  qs[i]*beta*v;
	    leftEigenvector[i][m+2] = -qs[i]*beta;
	}

	//2.1.b. Left-eigenvectors - bottom part
	for(j=0; j<m; j++)
	{
		leftEigenvector[m  ][j] =  v*k1x - u*k1y;
		leftEigenvector[m+1][j] =  pqs[j] - te1*c;
		leftEigenvector[m+2][j] =  pqs[j] + te1*c;
	}

    leftEigenvector[m][m]     =  k1y;
    leftEigenvector[m][m+1]   = -k1x;
    leftEigenvector[m][m+2]   =  0.;
    leftEigenvector[m+1][m]   =  c*k1x - beta*u;
    leftEigenvector[m+1][m+1] =  c*k1y - beta*v;
    leftEigenvector[m+1][m+2] =  beta;
    leftEigenvector[m+2][m]   = -c*k1x - beta*u;
    leftEigenvector[m+2][m+1] = -c*k1y - beta*v;
    leftEigenvector[m+2][m+2] =  beta;

    //2.2.a. Right-eigenvectors - upper part
    for(i=0; i<m; i++)
    {
	    for(j=0; j<m; j++)
	    {
	    	if(j == i)
	    		rightEigenvector[i][j] = rc2;
	    	else
	    		rightEigenvector[i][j] = 0.;
	    }
		rightEigenvector[i][m]   = 0.;
		rightEigenvector[i][m+1] = 0.5*qs[i]*rc2;
		rightEigenvector[i][m+2] = 0.5*qs[i]*rc2;
	}
    //2.2.b. Right-eigenvectors - bottom part
	for(j=0; j<m; j++)
	{
		rightEigenvector[m  ][j] = u*rc2;
		rightEigenvector[m+1][j] = v*rc2;
		rightEigenvector[m+2][j] = (u*u + v*v)*rc2 - pqs[j]/(beta*c*c + tiny);
	}
    rightEigenvector[m][m]     =  k1y;
    rightEigenvector[m][m+1]   =  0.5*(u + c*k1x)*rc2;
    rightEigenvector[m][m+2]   =  0.5*(u - c*k1x)*rc2;
    rightEigenvector[m+1][m]   = -k1x;
    rightEigenvector[m+1][m+1] =  0.5*(v + c*k1y)*rc2;
    rightEigenvector[m+1][m+2] =  0.5*(v - c*k1y)*rc2;
    rightEigenvector[m+2][m]   =  u*k1y - v*k1x;
    rightEigenvector[m+2][m+1] =  0.5*(H + c*te1)*rc2;
    rightEigenvector[m+2][m+2] =  0.5*(H - c*te1)*rc2;
}

/*---------------------------------------------------
 * allocate memory for calculation inviscid flux
 * ------------------------------------------------*/
void allocateFlux(int nlen, struct strct_flux *f)
{
	int i;

	f->xix  = (double*)malloc(sizeof(double)*nlen);
	f->xiy  = (double*)malloc(sizeof(double)*nlen);
	f->etx  = (double*)malloc(sizeof(double)*nlen);
	f->ety  = (double*)malloc(sizeof(double)*nlen);
	f->yas  = (double*)malloc(sizeof(double)*nlen);
	f->rho  = (double*)malloc(sizeof(double)*nlen);
	f->u    = (double*)malloc(sizeof(double)*nlen);
	f->v    = (double*)malloc(sizeof(double)*nlen);
	f->e    = (double*)malloc(sizeof(double)*nlen);
	f->p    = (double*)malloc(sizeof(double)*nlen);
	f->t    = (double*)malloc(sizeof(double)*nlen);
	f->gam  = (double*)malloc(sizeof(double)*nlen);
	f->qs   = (double**)malloc(sizeof(double*)*nlen);
	f->flux = (double**)malloc(sizeof(double*)*nlen);
	for(i=0; i<nlen; i++)
	{
		f->qs[i]     = (double*)malloc(sizeof(double)*config1.nspec);
		f->flux[i]   = (double*)malloc(sizeof(double)*neqn);
	}
}
/*---------------------------------------------------
 * free the memory
 * ------------------------------------------------*/
void freeFlux(int nlen, struct strct_flux *f)
{
	int i;

	free(f->xix);
	free(f->xiy);
	free(f->etx);
	free(f->ety);
	free(f->yas);
	free(f->rho);
	free(f->u);
	free(f->v);
	free(f->e);
	free(f->p);
	free(f->t);

	free(f->gam);

	for(i=0; i<nlen; i++)
	{
		free(f->qs[i]);
		free(f->flux[i]);
	}
	free(f->qs);
	free(f->flux);
}