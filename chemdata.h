/*
 * chemdata.h
 *
 *  Created on: Jan 21, 2014
 *
 */

#ifndef CHEMDATA_H_
#define CHEMDATA_H_

#define maxspec 6
#define maxreac 11
#define maxmix 2

#ifndef SCARF_DECLSPEC
#define SCARF_DECLSPEC extern
#endif

SCARF_DECLSPEC char spname[6][15];
SCARF_DECLSPEC struct strct_specData
{
	int  ntemrng, mdof;
	double wm, hof, tempmax, tempmin, temrng[3][2], acoef[3][10], qsin[2], vis[3];
} specData[6];

SCARF_DECLSPEC struct strct_reacData
{
	int thirdbody, backrate;
	double af, nf, thetaf, ab, nb, thetab, Br[maxspec],
			 thrdeff[maxspec], nrsr[maxspec], npsr[maxspec];
} reacData[11];

#endif /* CHEMDATA_H_ */