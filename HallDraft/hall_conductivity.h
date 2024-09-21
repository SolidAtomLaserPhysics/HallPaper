#ifndef HALL
#define HALL
#include "omp.h"
#include <pthread.h>
#include <complex> 
#include <vector>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstring>
using namespace std;
using Complex= complex<double>;
const int L =3;
//structure for the program parameters

struct Parameters{
  Parameters();
  //U - interaction strenght, mu - chemical potential, beta - inverse temperature, nAverage - electron density
  double U, mu, beta, nAverage;
  //Denominator and numerator of magnetic flux through the cell 
  int fluxDen, fluxNum;
  //Minimal and maximal value of omega
  double wmin,wmax,deltaW;
  //number of grid points
  int kxpoints, kypoints, wpoints;
  //treshholds value for imaginary part of sigma and denominator of Fermi prefactor of the integral.
  double treshSigma, treshFermi;
  //flags for calculation of density of states and sigma(0)
  bool denFlag, chernFlag;
  //order of pade fitting
  int pade_order;
  //array of pade coefficients
  Complex* pade_num,*pade_den;
};
//structure for matrices vx, vy (current operators) and e (energy in orbital space)

struct Matrices{
  Matrices();
   
  Complex vx[L*L], vy[L*L];
  double energy[L];
};

/* struct MatricesE{ */
/*   MatricesE(); */
/*   static int L; */
/*   MatrixXcd vx,vy; */
/*   VectorXd energy; */
/* }; */


void genEnMat(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky); 
void genVyMat(Complex* matrix, int fluxNum, int fluxDen, double ky);
void genVxMat(Complex* matrix, int fluxNum, int fluxDen, double kx);

/* void genEnMat(MatrixXcd& matrix, int fluxNum, int fluxDen, double kx, double ky); */
/* void genVyMat(MatrixXcd& matrix, int fluxNum, int fluxDen, double ky); */
/* void genVxMat(MatrixXcd& matrix, int fluxNum, int fluxDen, double kx); */

void genSigma(Complex* sigma);
void genDOS(double* density, Complex* sigma, Matrices* allMat);
//void genDOSE(double* density, Complex* sigma, MatricesE* allMat);
//double calcHallConductivityE(MatricesE* allMat, Complex* sigma);
double calcHallConductivity(Matrices* allMat, Complex* sigma,double* fermi);

void genEnMatT(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky); 
void genVyMatT(Complex* matrix, int fluxNum, int fluxDen, double ky);
void genVxMatT(Complex* matrix, int fluxNum, int fluxDen, double kx);

double fermi(double beta, double omega);
double fermi_pref(double beta, double w1, double w2);

//void fillMatricesEigen(MatricesE& matrices, double kx, double ky);
//void fillAllOMPEigen(MatricesE* allMat);
void fillAllOMPLapack(Matrices* allMat);
void precalculate_fermi_pref(double* fermi);
void printMatrix(Complex* matrix, int N, int M);
void printArray(Complex* array,int N);
void outDOS(double* dos);

extern "C" int zheev_(char *, char *, const  int *, Complex *,const int *, double *, Complex *,
		      int *, double *, int *);
extern "C" void zgemm_(char* transA,char* transB,const int* M, const int* N,const int* K, Complex* alpha,Complex* A,
		       const int* lda,Complex* B,const int* ldb, Complex* beta, Complex* C,const int* ldc);

#endif
