#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include <fftw3.h>
#include "random_number.h"

#define INICON  0.28       
#define LEN1 64   
#define LEN2 64
#define LEN3 64     
#define N 64
#define M 1

#define LOG1 6     
#define LOG2  6
#define STEP1 500        
#define STEP 30          
#define T    0.02      
#define G    0.1       
#define FRA  0.1           
#define FORBEGIN	for(i=0;i<LEN1;i++){  	for(j=0;j<LEN2;j++)		{ for(k=0;k<LEN2;k++)			{
#define FOREND }}}

#define 	XYZ  (i*LEN2*LEN3+j*LEN2+k)

//yxz
#define 	Za  (i*LEN2*LEN3+j*LEN2+(k+1)%LEN1)
#define 	Zd  (i*LEN2*LEN3+j*LEN2+(k-1+LEN1)%LEN1)

#define 	Xa  (i*LEN2*LEN3+((j+1)%LEN1)*LEN2+k )
#define 	Xd  (i*LEN2*LEN3+((j-1+LEN1)%LEN1)*LEN2+k)

#define 	Ya	(((i+1)%LEN1)*LEN2*LEN3+j*LEN2+k) 
#define 	Yd	(((i-1+LEN1)%LEN1)*LEN2*LEN3+j*LEN2+k) 


double	matrix[LEN1][LEN2][LEN3], matrix1[LEN1][LEN2][LEN3], hydro[LEN1][LEN2][LEN3],

dynamicX[LEN1][LEN2][LEN3 * 2], dynamicY[LEN1][LEN2][LEN3 * 2], dynamicZ[LEN1][LEN2][LEN3 * 2],
dynamic1X[LEN1][LEN2][LEN3], dynamic1Y[LEN1][LEN2][LEN3], dynamic1Z[LEN1][LEN2][LEN3],
dynamic2X[LEN1][LEN2][LEN3], dynamic2Y[LEN1][LEN2][LEN3], dynamic2Z[LEN1][LEN2][LEN3],

mindynamicX[LEN1][LEN2][LEN3], mindynamicY[LEN1][LEN2][LEN3], mindynamicZ[LEN1][LEN2][LEN3],

assidynamic1X[LEN1][LEN2][LEN3], assidynamic1Y[LEN1][LEN2][LEN3], assidynamic1Z[LEN1][LEN2][LEN3],
assidynamic2X[LEN1][LEN2][LEN3], assidynamic2Y[LEN1][LEN2][LEN3], assidynamic2Z[LEN1][LEN2][LEN3],

gradvelocity[LEN1][LEN2][LEN3][3][3], trangradvelocity[LEN1][LEN2][LEN3][3][3],

tensor[LEN1][LEN2][LEN3][3][3],
tensor1[LEN1][LEN2][LEN3][3][3], tensor1shear[LEN1][LEN2][LEN3][3][3], tensor1bulk[LEN1][LEN2][LEN3][3][3],
tensor2[LEN1][LEN2][LEN3][3][3], tensor2shear[LEN1][LEN2][LEN3][3][3], tensor2bulk[LEN1][LEN2][LEN3][3][3],

assistensor[LEN1][LEN2][LEN3][3][3],
gradtensor[LEN1][LEN2][LEN3][3][3],

forceos, forcesh, forcebu, forcea;

fftw_plan planX, planY, planZ, planRX, planRY, planRZ;


int InitialMtrix()
{
	int i, j, k;

	mtinit(8527);

	memset(tensor1shear, 0, sizeof(tensor1shear));
	memset(tensor1bulk, 0, sizeof(tensor1bulk));
	memset(tensor2shear, 0, sizeof(tensor2shear));
	memset(tensor2bulk, 0, sizeof(tensor2bulk));
	memset(dynamic1X, 0, sizeof(dynamic1X));
	memset(dynamic1Y, 0, sizeof(dynamic1Y));	//z 
	memset(dynamic2X, 0, sizeof(dynamic2X));
	memset(dynamic2Y, 0, sizeof(dynamic2Y));	//z

	FORBEGIN
		matrix[i][j][k] = INICON + 0.01 * nrand();
	FOREND

		return 0;
}

//K \tau_b  bulk time; Kgb bulk modulus  
void tensorbulknext(double* tens, double* dynaX, double* dynaY, double* dynaZ, double K, double Kgb, double tstep)
{
	int i, j, k;
	double relax, tb, Gb, tensorbulk[LEN1][LEN2][LEN3], meta[LEN1][LEN2][LEN3];
	//	memset(assistensor,0,sizeof(assistensor));	
	memset(meta, 0, sizeof(meta));


	FORBEGIN
		tensorbulk[i][j][k] = tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 0];
	FOREND

		FORBEGIN
		meta[i][j][k] = dynaX[i * LEN2 * LEN2 + j * LEN2 + k] * (tensorbulk[i][(j + 1) % LEN2][k] - tensorbulk[i][(LEN2 + j - 1) % LEN2][k]) / (-2);
	meta[i][j][k] += dynaY[i * LEN2 * LEN2 + j * LEN2 + k] * (tensorbulk[(i + 1) % LEN1][j][k] - tensorbulk[(LEN1 + i - 1) % LEN1][j][k]) / (-2);
	meta[i][j][k] += dynaZ[i * LEN2 * LEN2 + j * LEN2 + k] * (tensorbulk[i][j][(k + 1) % LEN1] - tensorbulk[i][j][(LEN1 + k - 1) % LEN1]) / (-2);
	FOREND


		FORBEGIN
		if (matrix[i][j][k] > INICON)
			Gb = Kgb;
		else
			Gb = 0;
	meta[i][j][k] += Gb * ((dynaX[(i * LEN2 + (j + 1) % LEN2) * LEN2 + k] - dynaX[(i * LEN2 + (LEN2 + j - 1) % LEN2) * LEN2 + k]) / 2);
	meta[i][j][k] += Gb * (dynaY[(((i + 1) % LEN1) * LEN2 + j) * LEN2 + k] - dynaY[(((LEN1 + i - 1) % LEN1) * LEN2 + j) * LEN2 + k]) / 2;
	meta[i][j][k] += Gb * (dynaZ[(i * LEN2 + j) * LEN2 + (k + 1) % LEN2] - dynaZ[(i * LEN2 + j) * LEN2 + (LEN2 + k - 1) % LEN2]) / 2;
	FOREND

		FORBEGIN
		tb = K * matrix[i][j][k] * matrix[i][j][k];
	relax = pow(2.71828, -tstep / tb);
	if (abs(K) < 0.01) relax = 0.;


	meta[i][j][k] *= tstep;
	tensorbulk[i][j][k] *= relax;
	tensorbulk[i][j][k] += meta[i][j][k];
	FOREND

		FORBEGIN
		tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 0] = tensorbulk[i][j][k];
	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 1] = 0;
	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 2] = 0;

	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 3] = 0;
	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 4] = tensorbulk[i][j][k];
	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 5] = 0;

	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 6] = 0;
	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 7] = 0;
	tens[(i * LEN2 * LEN2 + j * LEN2 + k) * 9 + 8] = tensorbulk[i][j][k];
	FOREND
}

void divtens(double* tens1, double* velX, double* velY, double* velZ)
{
	int i, j, k;
	FORBEGIN
		velX[XYZ] = (tens1[(Xa) * 9 + 0] - tens1[(Xd) * 9 + 0]) / 2 + (tens1[(Ya) * 9 + 3] - tens1[(Yd) * 9 + 3]) / 2 + (tens1[(Za) * 9 + 6] - tens1[(Zd) * 9 + 6]) / 2;
	velY[XYZ] = (tens1[(Xa) * 9 + 1] - tens1[(Xd) * 9 + 1]) / 2 + (tens1[(Ya) * 9 + 4] - tens1[(Yd) * 9 + 4]) / 2 + (tens1[(Za) * 9 + 7] - tens1[(Zd) * 9 + 7]) / 2;
	velZ[XYZ] = (tens1[(Xa) * 9 + 2] - tens1[(Xd) * 9 + 2]) / 2 + (tens1[(Ya) * 9 + 5] - tens1[(Yd) * 9 + 5]) / 2 + (tens1[(Za) * 9 + 8] - tens1[(Zd) * 9 + 8]) / 2;
	FOREND
}

void tensplus(double* tens1, double* tens2)
{
	int i, j, k;
	FORBEGIN
		tens1[(XYZ) * 9 + 0] = tens1[(XYZ) * 9 + 0] + tens2[(XYZ) * 9 + 0];
	tens1[(XYZ) * 9 + 1] = tens1[(XYZ) * 9 + 1] + tens2[(XYZ) * 9 + 1];
	tens1[(XYZ) * 9 + 2] = tens1[(XYZ) * 9 + 2] + tens2[(XYZ) * 9 + 2];
	tens1[(XYZ) * 9 + 3] = tens1[(XYZ) * 9 + 3] + tens2[(XYZ) * 9 + 3];
	tens1[(XYZ) * 9 + 4] = tens1[(XYZ) * 9 + 4] + tens2[(XYZ) * 9 + 4];
	tens1[(XYZ) * 9 + 5] = tens1[(XYZ) * 9 + 5] + tens2[(XYZ) * 9 + 5];
	tens1[(XYZ) * 9 + 6] = tens1[(XYZ) * 9 + 6] + tens2[(XYZ) * 9 + 6];
	tens1[(XYZ) * 9 + 7] = tens1[(XYZ) * 9 + 7] + tens2[(XYZ) * 9 + 7];
	tens1[(XYZ) * 9 + 8] = tens1[(XYZ) * 9 + 8] + tens2[(XYZ) * 9 + 8];
	FOREND
}

void gradtens(double* tens, double* dynaX, double* dynaY, double* dynaZ)
{
	int i, j, k;
	FORBEGIN
		gradtensor[i][j][k][0][0] = (tens[(Xa) * 9 + 0] - tens[(Xd) * 9 + 0]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 0] - tens[(Yd) * 9 + 0]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 0] - tens[(Zd) * 9 + 0]) / 2 * (dynaZ[XYZ]);

	gradtensor[i][j][k][0][1] = (tens[(Xa) * 9 + 1] - tens[(Xd) * 9 + 1]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 1] - tens[(Yd) * 9 + 1]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 1] - tens[(Zd) * 9 + 1]) / 2 * (dynaZ[XYZ]);

	gradtensor[i][j][k][0][2] = (tens[(Xa) * 9 + 2] - tens[(Xd) * 9 + 2]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 2] - tens[(Yd) * 9 + 2]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 2] - tens[(Zd) * 9 + 2]) / 2 * (dynaZ[XYZ]);



	gradtensor[i][j][k][1][0] = (tens[(Xa) * 9 + 3] - tens[(Xd) * 9 + 3]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 3] - tens[(Yd) * 9 + 3]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 3] - tens[(Zd) * 9 + 3]) / 2 * (dynaZ[XYZ]);

	gradtensor[i][j][k][1][1] = (tens[(Xa) * 9 + 4] - tens[(Xd) * 9 + 4]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 4] - tens[(Yd) * 9 + 4]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 4] - tens[(Zd) * 9 + 4]) / 2 * (dynaZ[XYZ]);

	gradtensor[i][j][k][1][2] = (tens[(Xa) * 9 + 5] - tens[(Xd) * 9 + 5]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 5] - tens[(Yd) * 9 + 5]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 5] - tens[(Zd) * 9 + 5]) / 2 * (dynaZ[XYZ]);



	gradtensor[i][j][k][2][0] = (tens[(Xa) * 9 + 6] - tens[(Xd) * 9 + 6]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 6] - tens[(Yd) * 9 + 6]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 6] - tens[(Zd) * 9 + 6]) / 2 * (dynaZ[XYZ]);

	gradtensor[i][j][k][2][1] = (tens[(Xa) * 9 + 7] - tens[(Xd) * 9 + 7]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 7] - tens[(Yd) * 9 + 7]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 7] - tens[(Zd) * 9 + 7]) / 2 * (dynaZ[XYZ]);

	gradtensor[i][j][k][2][2] = (tens[(Xa) * 9 + 8] - tens[(Xd) * 9 + 8]) / 2 * (dynaX[XYZ])
		+ (tens[(Ya) * 9 + 8] - tens[(Yd) * 9 + 8]) / 2 * (dynaY[XYZ])
		+ (tens[(Za) * 9 + 8] - tens[(Zd) * 9 + 8]) / 2 * (dynaZ[XYZ]);

	FOREND
}


void dot(double* tens1, double* tens2)
{
	int i, j, k;
	FORBEGIN
		assistensor[i][j][k][0][0] += tens1[(XYZ) * 9 + 0] * tens2[(XYZ) * 9 + 0] + tens1[(XYZ) * 9 + 1] * tens2[(XYZ) * 9 + 3] + tens1[(XYZ) * 9 + 2] * tens2[(XYZ) * 9 + 6];
	assistensor[i][j][k][0][1] += tens1[(XYZ) * 9 + 0] * tens2[(XYZ) * 9 + 1] + tens1[(XYZ) * 9 + 1] * tens2[(XYZ) * 9 + 4] + tens1[(XYZ) * 9 + 2] * tens2[(XYZ) * 9 + 7];
	assistensor[i][j][k][0][2] += tens1[(XYZ) * 9 + 0] * tens2[(XYZ) * 9 + 2] + tens1[(XYZ) * 9 + 1] * tens2[(XYZ) * 9 + 5] + tens1[(XYZ) * 9 + 2] * tens2[(XYZ) * 9 + 8];

	assistensor[i][j][k][1][0] += tens1[(XYZ) * 9 + 3] * tens2[(XYZ) * 9 + 0] + tens1[(XYZ) * 9 + 4] * tens2[(XYZ) * 9 + 3] + tens1[(XYZ) * 9 + 5] * tens2[(XYZ) * 9 + 6];
	assistensor[i][j][k][1][1] += tens1[(XYZ) * 9 + 3] * tens2[(XYZ) * 9 + 1] + tens1[(XYZ) * 9 + 4] * tens2[(XYZ) * 9 + 4] + tens1[(XYZ) * 9 + 5] * tens2[(XYZ) * 9 + 7];
	assistensor[i][j][k][1][2] += tens1[(XYZ) * 9 + 3] * tens2[(XYZ) * 9 + 2] + tens1[(XYZ) * 9 + 4] * tens2[(XYZ) * 9 + 5] + tens1[(XYZ) * 9 + 5] * tens2[(XYZ) * 9 + 8];

	assistensor[i][j][k][2][0] += tens1[(XYZ) * 9 + 6] * tens2[(XYZ) * 9 + 0] + tens1[(XYZ) * 9 + 7] * tens2[(XYZ) * 9 + 3] + tens1[(XYZ) * 9 + 8] * tens2[(XYZ) * 9 + 6];
	assistensor[i][j][k][2][1] += tens1[(XYZ) * 9 + 6] * tens2[(XYZ) * 9 + 1] + tens1[(XYZ) * 9 + 7] * tens2[(XYZ) * 9 + 4] + tens1[(XYZ) * 9 + 8] * tens2[(XYZ) * 9 + 7];
	assistensor[i][j][k][2][2] += tens1[(XYZ) * 9 + 6] * tens2[(XYZ) * 9 + 2] + tens1[(XYZ) * 9 + 7] * tens2[(XYZ) * 9 + 5] + tens1[(XYZ) * 9 + 8] * tens2[(XYZ) * 9 + 8];
	FOREND
}

void transpose(double* tensin, double* tensout)
{
	int i, j, k;
	FORBEGIN
		tensout[(XYZ) * 9 + 0] = tensin[(XYZ) * 9 + 0];
	tensout[(XYZ) * 9 + 1] = tensin[(XYZ) * 9 + 3];
	tensout[(XYZ) * 9 + 2] = tensin[(XYZ) * 9 + 6];

	tensout[(XYZ) * 9 + 3] = tensin[(XYZ) * 9 + 1];
	tensout[(XYZ) * 9 + 4] = tensin[(XYZ) * 9 + 4];
	tensout[(XYZ) * 9 + 5] = tensin[(XYZ) * 9 + 7];

	tensout[(XYZ) * 9 + 6] = tensin[(XYZ) * 9 + 2];
	tensout[(XYZ) * 9 + 7] = tensin[(XYZ) * 9 + 5];
	tensout[(XYZ) * 9 + 8] = tensin[(XYZ) * 9 + 8];

	FOREND
}

void gradvelo(double* dynamic1, double* dynamic2, double* dynamic3)
{
	int i, j, k;
	FORBEGIN
		gradvelocity[i][j][k][0][0] = (dynamic1[Xa] - dynamic1[Xd]) / 2;
	gradvelocity[i][j][k][0][1] = (dynamic2[Xa] - dynamic2[Xd]) / 2;
	gradvelocity[i][j][k][0][2] = (dynamic3[Xa] - dynamic3[Xd]) / 2;

	gradvelocity[i][j][k][1][0] = (dynamic1[Ya] - dynamic1[Yd]) / 2;
	gradvelocity[i][j][k][1][1] = (dynamic2[Ya] - dynamic2[Yd]) / 2;
	gradvelocity[i][j][k][1][2] = (dynamic3[Ya] - dynamic3[Yd]) / 2;

	gradvelocity[i][j][k][2][0] = (dynamic1[Za] - dynamic1[Zd]) / 2;
	gradvelocity[i][j][k][2][1] = (dynamic2[Za] - dynamic2[Zd]) / 2;
	gradvelocity[i][j][k][2][2] = (dynamic3[Za] - dynamic3[Zd]) / 2;

	FOREND
}

void tensorshearnext(double* tens, double* dynaX, double* dynaY, double* dynaZ, double K, double Kgs, double tstep)
{
	int i, j, k;
	double relax, Gs, ts, add;

	memset(assistensor, 0, sizeof(assistensor));


	gradtens(tens, dynaX, dynaY, dynaZ);

	FORBEGIN
		gradtensor[i][j][k][0][0] *= -1;
	gradtensor[i][j][k][0][1] *= -1;
	gradtensor[i][j][k][0][2] *= -1;

	gradtensor[i][j][k][1][0] *= -1;
	gradtensor[i][j][k][1][1] *= -1;
	gradtensor[i][j][k][1][2] *= -1;

	gradtensor[i][j][k][2][0] *= -1;
	gradtensor[i][j][k][2][1] *= -1;
	gradtensor[i][j][k][2][2] *= -1;
	FOREND


		gradvelo(dynaX, dynaY, dynaZ);
	transpose(&gradvelocity[0][0][0][0][0], &trangradvelocity[0][0][0][0][0]);
	dot(&gradvelocity[0][0][0][0][0], tens);
	dot(tens, &trangradvelocity[0][0][0][0][0]);

	tensplus(&gradvelocity[0][0][0][0][0], &trangradvelocity[0][0][0][0][0]);

	FORBEGIN
		Gs = Kgs * (matrix[i][j][k] * matrix[i][j][k]);

	gradvelocity[i][j][k][0][0] *= Gs;
	gradvelocity[i][j][k][0][1] *= Gs;
	gradvelocity[i][j][k][0][2] *= Gs;

	gradvelocity[i][j][k][1][0] *= Gs;
	gradvelocity[i][j][k][1][1] *= Gs;
	gradvelocity[i][j][k][1][2] *= Gs;

	gradvelocity[i][j][k][2][0] *= Gs;
	gradvelocity[i][j][k][2][1] *= Gs;
	gradvelocity[i][j][k][2][2] *= Gs;

	FOREND


		tensplus(&assistensor[0][0][0][0][0], &gradvelocity[0][0][0][0][0]);
	tensplus(&assistensor[0][0][0][0][0], &gradtensor[0][0][0][0][0]);

	FORBEGIN
		assistensor[i][j][k][0][0] *= tstep;
	assistensor[i][j][k][0][1] *= tstep;
	assistensor[i][j][k][0][2] *= tstep;

	assistensor[i][j][k][1][0] *= tstep;
	assistensor[i][j][k][1][1] *= tstep;
	assistensor[i][j][k][1][2] *= tstep;

	assistensor[i][j][k][2][0] *= tstep;
	assistensor[i][j][k][2][1] *= tstep;
	assistensor[i][j][k][2][2] *= tstep;

	FOREND

		FORBEGIN

		ts = K * matrix[i][j][k] * matrix[i][j][k];
	relax = pow(2.71828, -tstep / ts);
	if (abs(K) < 0.01) relax = 0.;

	tens[(XYZ) * 9 + 0] *= relax;
	tens[(XYZ) * 9 + 1] *= relax;
	tens[(XYZ) * 9 + 2] *= relax;

	tens[(XYZ) * 9 + 3] *= relax;
	tens[(XYZ) * 9 + 4] *= relax;
	tens[(XYZ) * 9 + 5] *= relax;

	tens[(XYZ) * 9 + 6] *= relax;
	tens[(XYZ) * 9 + 7] *= relax;
	tens[(XYZ) * 9 + 8] *= relax;
	FOREND

		tensplus(tens, &assistensor[0][0][0][0][0]);

	FORBEGIN
		add = 0;
	add = tens[(XYZ) * 9 + 0] + tens[(XYZ) * 9 + 4] + tens[(XYZ) * 9 + 8];
	add /= 3;
	tens[(XYZ) * 9 + 0] -= add;
	tens[(XYZ) * 9 + 4] -= add;
	tens[(XYZ) * 9 + 8] -= add;
	FOREND

}

void tensor1next()// p
{
	//no shear k = 0 kgs=0
	tensorshearnext(&tensor1shear[0][0][0][0][0], &dynamic1X[0][0][0], &dynamic1Y[0][0][0], &dynamic1Z[0][0][0], 0, 0.0, T);
	//tensorshearnext(&tensor1shear[0][0][0][0][0],&dynamic1X[0][0][0],&dynamic1Y[0][0][0],&dynamic1Z[0][0][0],10,0.2,T) ;

	//no bulk k=0 kgb =0
	tensorbulknext(&tensor1bulk[0][0][0][0][0], &dynamic1X[0][0][0], &dynamic1Y[0][0][0], &dynamic1Z[0][0][0], 0, 0, T);
	//tensorbulknext(&tensor1bulk[0][0][0][0][0],&dynamic1X[0][0][0],&dynamic1Y[0][0][0],&dynamic1Z[0][0][0],10,5,T);


	memset(tensor1, 0, sizeof(tensor1));
	tensplus(&tensor1[0][0][0][0][0], &tensor1shear[0][0][0][0][0]);
	tensplus(&tensor1[0][0][0][0][0], &tensor1bulk[0][0][0][0][0]);
}



int  HydroDynamic()
{
	double k2, complex[2], m, n, l;
	int i, j, k;
	memset(dynamicX, 0, sizeof(dynamicX));
	memset(dynamicY, 0, sizeof(dynamicY));
	memset(dynamicZ, 0, sizeof(dynamicZ));
	///
	divtens(&tensor1[0][0][0][0][0], &assidynamic1X[0][0][0], &assidynamic1Y[0][0][0], &assidynamic1Z[0][0][0]);

	forceos = 0;

	FORBEGIN
		dynamicX[i][j][2 * k] = matrix[i][j][k] * (hydro[i][(j + 1) % LEN2][k] - hydro[i][(LEN2 + j - 1) % LEN2][k]) / 2;
	dynamicY[i][j][2 * k] = matrix[i][j][k] * (hydro[(i + 1) % LEN1][j][k] - hydro[(LEN1 + i - 1) % LEN1][j][k]) / 2;
	dynamicZ[i][j][2 * k] = matrix[i][j][k] * (hydro[i][j][(k + 1) % LEN1] - hydro[i][j][(LEN1 + k - 1) % LEN1]) / 2;

	forceos += pow(dynamicX[i][j][2 * k] * dynamicX[i][j][2 * k] + dynamicY[i][j][2 * k] * dynamicY[i][j][2 * k] +
		dynamicZ[i][j][2 * k] * dynamicZ[i][j][2 * k], 0.5);



	mindynamicX[i][j][k] = -(1 - matrix[i][j][k]) / FRA * (dynamicX[i][j][2 * k] - assidynamic1X[i][j][k]);
	dynamicX[i][j][2 * k] = assidynamic1X[i][j][k] - dynamicX[i][j][2 * k];

	mindynamicY[i][j][k] = -(1 - matrix[i][j][k]) / FRA * (dynamicY[i][j][2 * k] - assidynamic1Y[i][j][k]);
	dynamicY[i][j][2 * k] = assidynamic1Y[i][j][k] - dynamicY[i][j][2 * k];

	mindynamicZ[i][j][k] = -(1 - matrix[i][j][k]) / FRA * (dynamicZ[i][j][2 * k] - assidynamic1Z[i][j][k]);
	dynamicZ[i][j][2 * k] = assidynamic1Z[i][j][k] - dynamicZ[i][j][2 * k];

	dynamicX[i][j][2 * k] *= pow(-1.0, i + j + k);
	dynamicY[i][j][2 * k] *= pow(-1.0, i + j + k);
	dynamicZ[i][j][2 * k] *= pow(-1.0, i + j + k);
	FOREND

		fftw_execute(planX);
	fftw_execute(planY);
	fftw_execute(planZ);

	forceos /= LEN1 * LEN2 * LEN3;

	FORBEGIN
		m = (i - LEN1 / 2.0) * 6.28 / (LEN1) * 2;
	n = (j - LEN2 / 2.0) * 6.28 / (LEN1) * 2;
	l = (k - LEN2 / 2.0) * 6.28 / (LEN1) * 2;

	k2 = (m * m + n * n + l * l);


	if (k2 == 0) {
		complex[0] = 0;
		complex[1] = 0; k2 = 1;
	}
	else {
		complex[0] = (n * dynamicX[i][j][2 * k] + m * dynamicY[i][j][2 * k] + l * dynamicZ[i][j][2 * k]) / k2;
		complex[1] = (n * dynamicX[i][j][2 * k + 1] + m * dynamicY[i][j][2 * k + 1] + l * dynamicZ[i][j][2 * k + 1]) / k2;
	}

	dynamicX[i][j][2 * k] -= n * complex[0];
	dynamicX[i][j][2 * k + 1] -= n * complex[1];
	dynamicX[i][j][2 * k] /= G * k2;
	dynamicX[i][j][2 * k + 1] /= G * k2;

	dynamicY[i][j][2 * k] -= m * complex[0];
	dynamicY[i][j][2 * k + 1] -= m * complex[1];
	dynamicY[i][j][2 * k] /= G * k2;
	dynamicY[i][j][2 * k + 1] /= G * k2;

	dynamicZ[i][j][2 * k] -= l * complex[0];
	dynamicZ[i][j][2 * k + 1] -= l * complex[1];
	dynamicZ[i][j][2 * k] /= G * k2;
	dynamicZ[i][j][2 * k + 1] /= G * k2;
	FOREND
		//  printf("all %f \n",dynamicX[0][0][1]*n+dynamicY[0][0][1]*m+dynamicZ[0][0][1]*l);
		fftw_execute(planRX);
	fftw_execute(planRY);
	fftw_execute(planRZ);

	FORBEGIN
		dynamicX[i][j][2 * k] *= pow(-1.0, i + j + k) / pow(N, 3.0);
	dynamicY[i][j][2 * k] *= pow(-1.0, i + j + k) / pow(N, 3.0);
	dynamicZ[i][j][2 * k] *= pow(-1.0, i + j + k) / pow(N, 3.0);

	dynamic1X[i][j][k] = dynamicX[i][j][2 * k] + (1 - matrix[i][j][k]) * mindynamicX[i][j][k];
	dynamic2X[i][j][k] = dynamicX[i][j][2 * k] - matrix[i][j][k] * mindynamicX[i][j][k];


	dynamic1Y[i][j][k] = dynamicY[i][j][2 * k] + (1 - matrix[i][j][k]) * mindynamicY[i][j][k];
	dynamic2Y[i][j][k] = dynamicY[i][j][2 * k] - matrix[i][j][k] * mindynamicY[i][j][k];

	dynamic1Z[i][j][k] = dynamicZ[i][j][2 * k] + (1 - matrix[i][j][k]) * mindynamicZ[i][j][k];
	dynamic2Z[i][j][k] = dynamicZ[i][j][2 * k] - matrix[i][j][k] * mindynamicZ[i][j][k];
	FOREND
		//	   printf("%f %f\n",dynamicX[0][0][0],dynamicX[0][0][1]/pow(N,3.0));
		//	 	printf("x%f ",-0.5*(dynamicX[1][0][2]-dynamicX[1][2][2]+dynamicY[0][1][2]-dynamicY[2][1][2]+dynamicZ[1][1][0]-dynamicZ[1][1][4]));
			  //	printf("y%f ",G*(2*dynamicX[1][1][2]-dynamicX[0][1][2]-dynamicX[2][1][2]+2*dynamicY[1][1][2]-dynamicY[0][1][2]-dynamicY[2][1][2]+2*dynamicZ[1][1][2]-dynamicZ[0][1][2]-dynamicZ[2][1][2]));
				//printf("z%f \n",G*(2*dynamicX[1][1][2]-dynamicX[1][1][4]-dynamicX[1][1][0]+2*dynamicY[1][1][2]-dynamicY[1][1][4]-dynamicY[1][1][0]+2*dynamicZ[1][1][2]-dynamicZ[1][1][4]-dynamicZ[1][1][0]));
		return 0;
}




int StepRun()
{
	int i, j, k;
	double add;
	memset(matrix1, 0, sizeof(matrix1));
	memset(hydro, 0, sizeof(hydro));

	FORBEGIN
		add = matrix[(i + 1) % LEN1][j][k] + matrix[(LEN1 + i - 1) % LEN1][j][k] + matrix[i][(j + 1) % LEN2][k]
		+ matrix[i][(LEN2 + j - 1) % LEN2][k] + matrix[i][j][(k + 1) % LEN2] + matrix[i][j][(LEN2 + k - 1) % LEN2];
	add = add - 6 * matrix[i][j][k];
	hydro[i][j][k] = 1.3 * (log(matrix[i][j][k] / (1 - matrix[i][j][k])) + 2.7 * (1 - 2 * matrix[i][j][k])) - add;
	FOREND

		HydroDynamic();

	FORBEGIN
		add = matrix[i][((j + 1) % LEN2)][k] * dynamic1X[i][((j + 1) % LEN2)][k] - matrix[i][((LEN2 + j - 1) % LEN2)][k] * dynamic1X[i][((LEN2 + j - 1) % LEN2)][k];
	add += matrix[(i + 1) % LEN1][j][k] * dynamic1Y[(i + 1) % LEN1][j][k] - matrix[(LEN1 + i - 1) % LEN1][j][k] * dynamic1Y[(LEN1 + i - 1) % LEN1][j][k];
	add += matrix[i][j][(k + 1) % LEN1] * dynamic1Z[i][j][((k + 1) % LEN1)] - matrix[i][j][(LEN1 + k - 1) % LEN1] * dynamic1Z[i][j][((LEN1 + k - 1) % LEN1)];
	add *= T / -2.0;
	matrix1[i][j][k] = add;
	FOREND

		FORBEGIN
		matrix[i][j][k] += matrix1[i][j][k];
	FOREND

		return 0;
}

double size(double* condense)
{
	int i, j, k;
	double q1, q2, q3, q4, sta[2 * N * M][2];
	fftw_complex data[N * M][N * M][N * M];
	fftw_plan plan;

	memset(data, 0, sizeof(data));
	memset(sta, 0, sizeof(sta));
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			for (k = 0; k < N; k++)
			{
				data[i + (M * N - N) / 2][j + (M * N - N) / 2][k + (M * N - N) / 2][0] = condense[i * N * N + j * N + k];
			}

	for (i = 0; i < N * M; i++)
		for (j = 0; j < N * M; j++)
			for (k = 0; k < N * M; k++)
			{
				data[i][j][k][0] *= pow(-1.0, i + j + k);
				data[i][j][k][1] = 0.0;
			}

	plan = fftw_plan_dft_3d(N * M, N * M, N * M, &data[0][0][0], &data[0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);


	for (i = 0; i < N * M; i++)
		for (j = 0; j < N * M; j++)
			for (k = 0; k < N * M; k++)
			{
				data[i][j][k][0] = pow(data[i][j][k][0], 2.0) + pow(data[i][j][k][1], 2.0);
				q1 = i - N * M / 2;
				q2 = j - N * M / 2;
				q3 = k - N * M / 2;
				q4 = q1 * q1 + q2 * q2 + q3 * q3;
				q4 = pow(q4, 0.5);

				sta[(int)(q4)][0] += data[i][j][k][0];
				sta[(int)(q4)][1] += 1;
			}
	q1 = 0; q2 = 0; q3 = 0; q4 = 0;

	for (i = M; i < 2 * M * N; i++)
	{
		if (sta[i][1] != 0)
		{
			sta[i][0] = sta[i][0] / sta[i][1];
			q3 += i / M * 1.0 * sta[i][0];
			q4 += sta[i][0];
			//	q1+=i/M*1.0*data[N*M/2][N*M/2][i+M*N/2][0];
				//q2+=data[N*M/2][N*M/2][i+M*N/2][0];   
		}

	}


	fftw_destroy_plan(plan);
	return q3 / q4;
}

void output(char* filename, double* out, double* v1x, double* v1y, double* v1z)
{
	FILE* file;
	int i, j, k;

	file = fopen(filename, "w");
	fprintf(file, "Title=VE\n Variables=x y z dens v1x v1y vz \n");
	fprintf(file, "Zone I=%d J=%d  K=%d f=point\n", LEN1, LEN2, LEN3);
	FORBEGIN
		fprintf(file, "%6.2f %6.2f %6.2f  %6.4f %6.4f %6.4f %6.4f\n",
			1.0 * i, 1.0 * j, 1.0 * k,
			out[XYZ], v1x[XYZ], v1y[XYZ], v1z[XYZ]);
	FOREND
		fflush(file);
	fclose(file);
}

int main(int argc, char* argv[])
{
	clock_t   Time1, Time2;
	int i, j, k, z;
	FILE* fp2;
	double op;

	srand((int)time(0));
	InitialMtrix();


	planX = fftw_plan_dft_3d(N * M, N * M, N * M, (fftw_complex*)&dynamicX[0][0][0], (fftw_complex*)&dynamicX[0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	planY = fftw_plan_dft_3d(N * M, N * M, N * M, (fftw_complex*)&dynamicY[0][0][0], (fftw_complex*)&dynamicY[0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	planZ = fftw_plan_dft_3d(N * M, N * M, N * M, (fftw_complex*)&dynamicZ[0][0][0], (fftw_complex*)&dynamicZ[0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	planRX = fftw_plan_dft_3d(N * M, N * M, N * M, (fftw_complex*)&dynamicX[0][0][0], (fftw_complex*)&dynamicX[0][0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	planRY = fftw_plan_dft_3d(N * M, N * M, N * M, (fftw_complex*)&dynamicY[0][0][0], (fftw_complex*)&dynamicY[0][0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	planRZ = fftw_plan_dft_3d(N * M, N * M, N * M, (fftw_complex*)&dynamicZ[0][0][0], (fftw_complex*)&dynamicZ[0][0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	fp2 = fopen("experiment.dat", "wt+");
	Time1 = clock();


	for (z = 1; z < 100001; z++)
	{

		printf("%d\t%f\n", z, matrix[0][0][0]);

		if (z % 50 == 0) {

			forcesh = 0;
			forcebu = 0;
			forcea = 0;

			divtens(&tensor1[0][0][0][0][0], &assidynamic1X[0][0][0], &assidynamic1Y[0][0][0], &assidynamic1Z[0][0][0]);
			FORBEGIN
				forcea += pow(assidynamic1X[i][j][k] * assidynamic1X[i][j][k] + assidynamic1Y[i][j][k] * assidynamic1Y[i][j][k] + assidynamic1Z[i][j][k] * assidynamic1Z[i][j][k], 0.5);
			FOREND
				forcea /= LEN1 * LEN2 * LEN3;

			divtens(&tensor1shear[0][0][0][0][0], &assidynamic1X[0][0][0], &assidynamic1Y[0][0][0], &assidynamic1Z[0][0][0]);
			FORBEGIN
				forcesh += pow(assidynamic1X[i][j][k] * assidynamic1X[i][j][k] + assidynamic1Y[i][j][k] * assidynamic1Y[i][j][k] + assidynamic1Z[i][j][k] * assidynamic1Z[i][j][k], 0.5);
			FOREND
				forcesh /= LEN1 * LEN2 * LEN3;

			divtens(&tensor1bulk[0][0][0][0][0], &assidynamic1X[0][0][0], &assidynamic1Y[0][0][0], &assidynamic1Z[0][0][0]);
			FORBEGIN
				forcebu += pow(assidynamic1X[i][j][k] * assidynamic1X[i][j][k] + assidynamic1Y[i][j][k] * assidynamic1Y[i][j][k] + assidynamic1Z[i][j][k] * assidynamic1Z[i][j][k], 0.5);
			FOREND
				forcebu /= LEN1 * LEN2 * LEN3;

			op = 0.;
			FORBEGIN
				op += fabs(matrix[i][j][k] - INICON);
			FOREND
				op /= LEN1 * LEN2 * LEN3;

			fprintf(fp2, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", z, size(&matrix[0][0][0]), forcea, forcebu, forcesh, forceos, op); fflush(fp2);

		}
		if (z % 1000 == 0)
		{
			char filename[50]; sprintf(filename, "dens_%d.dat", z);
			output(filename, &matrix[0][0][0], &dynamic1X[0][0][0], &dynamic1Y[0][0][0], &dynamic1Z[0][0][0]);
		}

		tensor1next();
		StepRun();

	}

	Time2 = clock(); printf("Run   code   time   is:     %lf h\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC / 3600.);
	fclose(fp2);

	fftw_destroy_plan(planX);
	fftw_destroy_plan(planY);
	fftw_destroy_plan(planZ);
	fftw_destroy_plan(planRX);
	fftw_destroy_plan(planRY);
	fftw_destroy_plan(planRZ);


	return 0;
}


