#include "hall_conductivity.h"


using namespace std;
using Clock=std::chrono::high_resolution_clock;


//Some notes after talking to Anton on 20.9:
  //TODO: want to read out Sigma in genSigma instead of reading in pade coefficients

  //genDOS is just postprocessing and do not need it/have not changed it so probably will not work
  //run by:
  //  make
  //  ./hall_conductivity
  //better run on node with large memory
  //only care about mu, not n



//define necessary constants;

const double t = 0.25;
const double tPrime = -0.05;
const double tPrimePrime = 0.025;
const Complex I = Complex(0.0, 1.0);
const double pi = acos(-1.0);


//read parameters and initialize the size of the orbital space 
Parameters parameters;



int main(int argv, char* argc[]){
  cout << "wmin:" << parameters.wmin << " \n"; 
  cout << "wmax:" << parameters.wmax << " \n"; 
  cout << "wpoints:" << parameters.wpoints << " \n"; 
  cout << "kxPoints:" << parameters.kxpoints << " \n"; 
  cout << "density:" << parameters.nAverage << " \n"; 
  cout << "beta:" << parameters.beta << " \n"; 
  cout << "pi:" << pi << " \n"; 
  Matrices* allMat = new Matrices[parameters.kxpoints*parameters.kypoints];
  Complex* sigma = new Complex[parameters.wpoints];
  double* fermi_pref = new double[parameters.wpoints*parameters.wpoints]; 
  genSigma(sigma);
  fillAllOMPLapack(allMat);
  precalculate_fermi_pref(fermi_pref);
  auto  start = Clock::now(); 
  double sigma2 = calcHallConductivity(allMat, sigma,fermi_pref);
  auto finish = Clock::now();
  auto time = chrono::duration_cast<std::chrono::seconds>(finish - start).count();  
  cout <<"Calculation took:" << time << " s \n"; 
  cout <<"\n" << sigma2 <<"\n";
  delete[] allMat;
  delete[] sigma;
  delete[] fermi_pref;
}

//============Generate matrices in orbital space. This is necessary to send a matrix inicialized with zero values!!!==========

//generate dispersion matrix
//multiplied with - since that is my convention
void genEnMat(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky){ 
  //int L = fluxDen;
  //for (int i = 1; i<L; i++){
  //  matrix[(i-1)*L + i] = -t;
  //  matrix[i*L + i - 1] = -t;
  //  matrix[i*L + i] = -2.0*t*cos(kx + 2 * pi*i*fluxNum / fluxDen);
  //}
  //matrix[0] = -2.0*t*cos(kx);
  //matrix[fluxDen - 1] = -t*exp(I*double(fluxDen)*ky);
  //matrix[(fluxDen - 1)*L] = -t*exp(-I*double(fluxDen)*ky);



  //write epsilon for L=3 to make it easier
  if (L == 3){
    matrix[0] = 2.0*t*cos(kx) + 2.0*tPrimePrime*cos(2.0*kx);
    matrix[1] = t + 2.0*tPrime*cos(kx + pi*fluxNum/ fluxDen) + tPrimePrime*exp(I*double(fluxDen)*ky);
    matrix[2] = t*exp(I*double(fluxDen)*ky) + 2.0*tPrime*cos(kx - pi*fluxNum/ fluxDen)*exp(I*double(fluxDen)*ky) + tPrimePrime;

    matrix[3] = t + 2.0*tPrime*cos(kx + pi*fluxNum/ fluxDen) + tPrimePrime*exp(-I*double(fluxDen)*ky);
    matrix[4] = 2.0*t*cos(kx + 2.0*pi*fluxNum / fluxDen) + 2.0*tPrimePrime*cos(2.0*(kx + 2.0*pi*fluxNum / fluxDen));
    matrix[5] = t + 2.0*tPrime*cos(kx + 3.0*pi*fluxNum/ fluxDen) + tPrimePrime*exp(I*double(fluxDen)*ky);

    matrix[6] = t*exp(-I*double(fluxDen)*ky) + 2.0*tPrime*cos(kx - pi*fluxNum/ fluxDen)*exp(-I*double(fluxDen)*ky) + tPrimePrime;
    matrix[7] = t + 2.0*tPrime*cos(kx + 3.0*pi*fluxNum/ fluxDen) + tPrimePrime*exp(-I*double(fluxDen)*ky);
    matrix[8] = 2.0*t*cos(kx + 2.0*pi*2.0*fluxNum / fluxDen) + 2.0*tPrimePrime*cos(2.0*(kx + 4.0*pi*fluxNum / fluxDen));
  }
}

//Generate velocity matrices 

//In y direction

//v = del epsilon/ del ky - i(position of l in y - position of lPrime in y)epsilon l lPrime
//only e to power ky survive del and 
//            unit cell high so that second term survives depending on l
//            diagonale is l=lPrime and term falls away
//            nextdiagonale survives with + or - dependent on where in matrix
//            edge should also have q as factor
//            others fall awy since eps = 0

void genVyMat(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky){
  ////int L = fluxDen;
  //for (int i = 1; i < L ; i++){
  //  matrix[(i - 1)*L + i] = I*t;
  //  matrix[i*L + i - 1] = -I*t;
  //}

  //matrix[fluxDen - 1] = -I*t*exp(I*double(fluxDen)*ky);
  //matrix[(fluxDen - 1)*L] = I*t*exp(-I*double(fluxDen)*ky);
  //write epsilon for L=3 to make it easier
  if (L == 3){
    matrix[0] = 0.0;
    matrix[1] = -I*t - 2.0*I*tPrime*cos(kx + pi*fluxNum/double(fluxDen)) + 2.0*I*tPrimePrime*exp(I*ky*double(fluxDen));
    matrix[2] = I*t*exp(I*ky*double(fluxDen)) + 2.0*I*tPrime*cos(kx - pi*fluxNum/double(fluxDen))*exp(I*ky*double(fluxDen)) - 2.0*I*tPrimePrime;

    matrix[3] = I*t + 2.0*I*tPrime*cos(kx + pi*fluxNum/double(fluxDen)) - 2.0*I*tPrimePrime*exp(-I*ky*double(fluxDen));
    matrix[4] = 0.0;
    matrix[5] = -I*t - 2.0*I*tPrime*cos(kx + 3*pi*fluxNum/double(fluxDen)) + 2.0*I*tPrimePrime*exp(I*ky*double(fluxDen));

    matrix[6] = -I*t*exp(-I*ky*double(fluxDen)) - 2.0*I*tPrime*cos(kx - pi*fluxNum/double(fluxDen))*exp(-I*ky*double(fluxDen)) + 2.0*I*tPrimePrime;
    matrix[7] = I*t + 2.0*I*tPrime*cos(kx + 3*pi*fluxNum/double(fluxDen)) - 2.0*I*tPrimePrime*exp(-I*ky*double(fluxDen));
    matrix[8] = 0.0;
  }
}

//And in x - direction

//v = del epsilon/ del kx - i(position of l in x - position of lPrime in x)epsilon l lPrime
//only sin(kx) survive del and 
//            since only one thick in x direction, second term always falls away

void genVxMat(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky){
  //int L = fluxDen;
  //for (int i = 0; i < fluxDen; i++){
  //  matrix[i*L + i] = 2.0 * t * sin(kx + 2 * pi*i*fluxNum / fluxDen);
  //}

  if (L == 3){
    matrix[0] = - 2.0*t*sin(kx) - 4.0*tPrimePrime*sin(2.0*kx);
    matrix[1] = - 2.0*tPrime*sin(kx + pi*fluxNum/double(fluxDen));
    matrix[2] = - 2.0*tPrime*sin(kx - pi*fluxNum/double(fluxDen)) * exp(I*ky*double(fluxDen));

    matrix[3] = - 2.0*tPrime*sin(kx + pi*fluxNum/double(fluxDen));
    matrix[4] = -2.0*t*sin(kx + 2.0*pi*fluxNum/double(fluxDen)) - 4.0*tPrimePrime*sin(2.0*kx + 4*pi*fluxNum/double(fluxDen));
    matrix[5] = -2.0*tPrime*sin(kx + 3*pi*fluxNum/double(fluxDen));

    matrix[6] = - 2.0*tPrime *sin(kx - pi*fluxNum/double(fluxDen))*exp(-I*ky*double(fluxDen));
    matrix[7] = - 2.0*tPrime*sin(kx + 3*pi*fluxNum/double(fluxDen));
    matrix[8] = - 2.0*t*sin(kx + 4*pi*fluxNum/double(fluxDen)) - 4.0*tPrimePrime*cos(2.0*kx + 8*pi*fluxNum/double(fluxDen));
  }
}
//====================And... Again gengerate matrices. This time transpose for lapack. For test.=====

//generate dispersion matrix
//only change exponents of exp(I*ky*q)

void genEnMatT(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky){ 
  //int L = fluxDen;
  //for (int i = 1; i<L; i++){
  //  matrix[(i-1)*L + i] = -t;
  //  matrix[i*L + i - 1] = -t;
  //  matrix[i*L + i] = -2.0*t*cos(kx + 2 * pi*i*fluxNum / fluxDen);
  //}
  //matrix[0] = -2.0*t*cos(kx);
  //matrix[fluxDen - 1] = -t*exp(-I*double(fluxDen)*ky);
  //matrix[(fluxDen - 1)*L] = -t*exp(I*double(fluxDen)*ky);

  //write epsilon for L=3 to make it easier
  if (L == 3){
    matrix[0] = 2.0*t*cos(kx) + 2.0*tPrimePrime*cos(2.0*kx);
    matrix[1] = t + 2.0*tPrime*cos(kx + pi*fluxNum/ fluxDen) + tPrimePrime*exp(-I*double(fluxDen)*ky);
    matrix[2] = t*exp(-I*double(fluxDen)*ky) + 2.0*tPrime*cos(kx - pi*fluxNum/ fluxDen)*exp(-I*double(fluxDen)*ky) + tPrimePrime;

    matrix[3] = t + 2.0*tPrime*cos(kx + pi*fluxNum/ fluxDen) + tPrimePrime*exp(I*double(fluxDen)*ky);
    matrix[4] = 2.0*t*cos(kx + 2.0*pi*fluxNum / fluxDen) + 2.0*tPrimePrime*cos(2.0*(kx + 2.0*pi*fluxNum / fluxDen));
    matrix[5] = t + 2.0*tPrime*cos(kx + 3.0*pi*fluxNum/ fluxDen) + tPrimePrime*exp(-I*double(fluxDen)*ky);

    matrix[6] = t*exp(I*double(fluxDen)*ky) + 2.0*tPrime*cos(kx - pi*fluxNum/ fluxDen)*exp(I*double(fluxDen)*ky) + tPrimePrime;
    matrix[7] = t + 2.0*tPrime*cos(kx + 3.0*pi*fluxNum/ fluxDen) + tPrimePrime*exp(I*double(fluxDen)*ky);
    matrix[8] = 2.0*t*cos(kx + 2.0*pi*2.0*fluxNum / fluxDen) + 2.0*tPrimePrime*cos(2.0*(kx + 4.0*pi*fluxNum / fluxDen));
  }
}

//Generate velocity matrices 

//In y direction
void genVyMatT(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky){
  //int L = fluxDen;
  //for (int i = 1; i < L ; i++){
  //  matrix[(i - 1)*L + i] = -I*t;
  //  matrix[i*L + i - 1] = I*t;
  //}

  if (L == 3){
    matrix[0] = 0.0;
    matrix[1] = I*t + 2.0*I*tPrime*cos(kx + pi*fluxNum/double(fluxDen)) - 2.0*I*tPrimePrime*exp(-I*ky*double(fluxDen));
    matrix[2] = -I*t*exp(-I*ky*double(fluxDen)) - 2.0*I*tPrime*cos(kx - pi*fluxNum/double(fluxDen))*exp(-I*ky*double(fluxDen)) + 2.0*I*tPrimePrime;

    matrix[3] = -I*t - 2.0*I*tPrime*cos(kx + pi*fluxNum/double(fluxDen)) + 2.0*I*tPrimePrime*exp(I*ky*double(fluxDen));
    matrix[4] = 0.0;
    matrix[5] = I*t + 2.0*I*tPrime*cos(kx + 3*pi*fluxNum/double(fluxDen)) - 2.0*I*tPrimePrime*exp(-I*ky*double(fluxDen));

    matrix[6] = I*t*exp(I*ky*double(fluxDen)) + 2.0*I*tPrime*cos(kx - pi*fluxNum/double(fluxDen))*exp(I*ky*double(fluxDen)) - 2.0*I*tPrimePrime;
    matrix[7] = -I*t - 2.0*I*tPrime*cos(kx + 3*pi*fluxNum/double(fluxDen)) + 2.0*I*tPrimePrime*exp(I*ky*double(fluxDen));
    matrix[8] = 0.0;
  }
}

//And in x - direction
void genVxMatT(Complex* matrix, int fluxNum, int fluxDen, double kx, double ky){
  //int L = fluxDen;
  //for (int i = 0; i < fluxDen; i++){
  //  matrix[i*L + i] = 2 * t * sin(kx + 2 * pi*i*fluxNum / fluxDen);
  //}

  if (L == 3){
    cout << "I am inside my vx:" << L << " \n"; 
    matrix[0] = - 2.0*t*sin(kx) - 4.0*tPrimePrime*sin(2.0*kx);
    matrix[1] = - 2.0*tPrime*sin(kx + pi*fluxNum/double(fluxDen));
    matrix[2] = - 2.0*tPrime*sin(kx - pi*fluxNum/double(fluxDen)) * exp(-I*ky*double(fluxDen));

    matrix[3] = - 2.0*tPrime*sin(kx + pi*fluxNum/double(fluxDen));
    matrix[4] = -2.0*t*sin(kx + 2.0*pi*fluxNum/double(fluxDen)) - 4.0*tPrimePrime*sin(2.0*kx + 4*pi*fluxNum/double(fluxDen));
    matrix[5] = -2.0*tPrime*sin(kx + 3*pi*fluxNum/double(fluxDen));

    matrix[6] = - 2.0*tPrime *sin(kx - pi*fluxNum/double(fluxDen))*exp(I*ky*double(fluxDen));
    matrix[7] = - 2.0*tPrime*sin(kx + 3*pi*fluxNum/double(fluxDen));
    matrix[8] = - 2.0*t*sin(kx + 4*pi*fluxNum/double(fluxDen)) - 4.0*tPrimePrime*cos(2.0*kx + 8*pi*fluxNum/double(fluxDen));
  }
}

//====================Calculation of matrices in energy-diagonal basis=========================
void fillAllOMPLapack(Matrices* allMat){
#pragma omp parallel 
{
  //define variables for lapack routines. For details go to LAPACK User's Guide
  char jobs = 'V';
  char uplo = 'U';
  char notrans = 'N'; 
  char conjTr = 'C';
 
  int info;
  //int L = parameters.fluxDen;
  double* w = new double[L]();
  int lwork = 2*L-1;  
  Complex* work = new Complex[lwork]();
  double* rwork = new double[3*L-2]();
 
  Complex* auxEnMat = new Complex[L*L]();
  Complex alpha = Complex(1,0);
  Complex beta = Complex(0,0);
  
  Complex* auxVmat = new Complex[L*L]();
  Complex* multMat = new Complex[L*L]();
 
  double deltaKx,kx,deltaKy,ky;
 
  deltaKx = (2*pi)/(parameters.kxpoints);
  deltaKy = (2*pi)/(parameters.fluxDen*parameters.kypoints);        //  2PI/qky
#pragma omp for 
  for(int i=0; i<parameters.kxpoints; i++){
    kx = -pi  + deltaKx*i;
    for(int j=0;j<parameters.kypoints; j++){
      ky =-pi/parameters.fluxDen + deltaKy*j;     
      //gengerate despersion matrix, find its eigenvalues and eigenvectors
      fill(auxEnMat,auxEnMat+L*L,0);
      genEnMatT(auxEnMat, parameters.fluxNum, parameters.fluxDen,kx, ky);   
      zheev_(&jobs, &uplo,&L ,auxEnMat, &L,  allMat[i*parameters.kypoints + j].energy, work, &lwork, rwork, &info);
      //transfrom vx and vy
     
      //gengerate vx matrix and transform it
      fill(auxVmat,auxVmat+L*L,0);
      genVxMatT(auxVmat,parameters.fluxNum,parameters.fluxDen,kx, ky);
      zgemm_(&conjTr,&notrans,&L, &L, &L, &alpha, auxEnMat,&L,auxVmat, &L, &beta, multMat,&L);
      zgemm_(&notrans,&notrans,&L, &L, &L, &alpha, multMat,&L,auxEnMat, &L, &beta, allMat[i*parameters.kypoints + j].vx, &L);  
  
      //gengerate vy matrix and transform it
      fill(auxVmat,auxVmat+L*L,0);
      genVyMatT(auxVmat,parameters.fluxNum,parameters.fluxDen,kx, ky);
      zgemm_(&conjTr,&notrans,&L, &L, &L, &alpha, auxEnMat,&L,auxVmat, &L, &beta, multMat,&L);
      zgemm_(&notrans,&notrans,&L, &L, &L, &alpha, multMat,&L,auxEnMat, &L, &beta, allMat[i*parameters.kypoints + j].vy, &L);
    }
  }

  //freem allocated memory
  delete[] w;
  delete[] work;
  delete[] rwork;
  delete[] auxEnMat;
  delete[] auxVmat;    
  delete[] multMat;
}
}
//Fermi distribution
double fermi(double beta, double omega){
  double result;
  if (beta == 0){
    result = (omega < 0) ? 1.0 : 0.0;
  }
  else{
    result = 1 / (exp(beta*omega) + 1);
  }
  return result;
}
//calculate temperature dependent factor in the integral.  
double fermi_pref(double beta, double w1, double w2){
  double fermi_pref;
  fermi_pref = (fermi(beta, w1) - fermi(beta, w2)) / (pow((w1 - w2), 2) + pow(parameters.treshFermi,2));
  return fermi_pref;
}


//read in self-energy
//for U=0 no changes, can therefore not be used
void genSigma(Complex* sigma){
  if(parameters.U==0){
    //in noninteracting case, fill array with treshSigma value
    fill(sigma,sigma + parameters.wpoints,parameters.treshSigma*I);
    cout << "sigmaU0:" << sigma << " \n"; 
  }
  else{
    //read file
    fstream Sigma_file("selfEnergyHall.dat",fstream::in);
      double temp;
      double trash;
      int wPoint = 0;

    std::string line;
    while (std::getline(Sigma_file, line)) {
        // Use stringstream to parse the line
        std::istringstream iss(line);
        // Variables to store the values from the second and third columns, which are real and imag part
        double real, imag;
        // Read values from the stringstream
        if (iss >> trash >> real >> imag) {
            // Add the values to the vector
          sigma[wPoint]= real + I*imag;
          wPoint++;
        }
      //while(!Sigma_file.eof()){
        //first line is frequency and do not want it
      //  if(Sigma_file){
          //then read in Im(Sigma)
      //    Sigma_file>>temp;
      //    sigma[wPoint]=temp;
          //last line is Re(Sigma) = const and do not need it
          

    }
    Sigma_file.close(); // Close the file when done
    cout << "sigma:" << sigma[3] << " \n"; 
  }
}


//Calcualte density of states 
//only postprocessing, do not need it
void genDOS(double* density, Complex* sigma, Matrices* allMat){
#pragma omp parallel
  {
  double densTotal;
  Complex dens;
  double omega, deltaW;
  double norm;
  deltaW = double(parameters.wmax - parameters.wmin)/double(parameters.wpoints - 1);
  #pragma omp for
  for(int w_count=0;w_count < parameters.wpoints;w_count++){
    omega = parameters.wmin + w_count*deltaW;
    densTotal = 0;
    for(int kx_count=0; kx_count < parameters.kxpoints;kx_count++){
      for(int ky_count = 0; ky_count <parameters.kypoints;ky_count++){
	dens = 0;
	for( int i_count=0; i_count < parameters.fluxDen;i_count++){
	  dens += -1.0/(omega - allMat[parameters.kypoints*kx_count+ky_count].energy[i_count]-sigma[w_count]
		     +parameters.mu-parameters.U*parameters.nAverage/2.0);
	}
	densTotal +=dens.imag(); 
      }
    }
    norm =double(4*pi)/double(parameters.kxpoints*parameters.kypoints*parameters.fluxDen*parameters.fluxDen);
    density[w_count] = densTotal*norm;
  }
}
}

//====================================Calculate hall conductivity==============================

double calcHallConductivity(Matrices* allMat, Complex* sigma, double* fermi_pref){
  //fstream kGrid;
  //kGrid.open("kGrid.dat", fstream::out);
  //kGrid <<  ;
  //kGrid.close();

  double norm;
  double result = 0;
   #pragma omp parallel
    {
    double deltaW, omega1, omega2;
    double aMat1[L], aMat2[L];
    double temp_result = 0;
    Matrices* mat;
    deltaW = double(parameters.wmax - parameters.wmin)/double(parameters.wpoints - 1);
    #pragma omp for reduction(+:result) 
    for(int kx_count=0;kx_count<parameters.kxpoints;kx_count++){     
      for(int ky_count=0;ky_count<parameters.kypoints;ky_count++){    
       
	mat = &allMat[parameters.kypoints*kx_count+ky_count];
     	for(int w1_count=0;w1_count<parameters.wpoints;w1_count++){
	  omega1 = parameters.wmin + w1_count*deltaW;

	  #pragma unroll
	  for(int k_count=0; k_count < L;k_count++){
	    aMat1[k_count] = imag( -1.0/(omega1 - mat->energy[k_count]-sigma[w1_count]+parameters.mu
					-parameters.U*parameters.nAverage/2.0));
	  }
	  for(int w2_count=0;w2_count<parameters.wpoints;w2_count++){
	    temp_result = 0;
	    omega2 = parameters.wmin + w2_count*deltaW;
	    #pragma unroll
	    for(int i_count=0; i_count < L;i_count++){
	      aMat2[i_count] = imag( -1.0/(omega2 - mat->energy[i_count]-sigma[w2_count] + parameters.mu
					   -parameters.U*parameters.nAverage/2.0));
	      #pragma unroll
	      for(int j_count=0;j_count < L;j_count++){ 
		temp_result += imag(mat->vx[L*j_count+i_count]*aMat1[j_count]*mat->vy[L*i_count+j_count]*aMat2[i_count]);
	      }
	    }
	    
	    temp_result *= fermi_pref[w1_count*parameters.wpoints + w2_count];
	    result +=temp_result/parameters.wpoints;  
	  }
	}
      }
    }
    norm = 2 * 2 * pow(parameters.wmax - parameters.wmin,2)/(pi*double(L*parameters.wpoints*
								       parameters.kxpoints*parameters.kypoints));
    }
  result*=norm;
  return result;  
}



//======================auxilary functions for classes, input and output =========================


//here we load parameters.dat but also coefficients.dat
Parameters::Parameters(){
  fstream param;
  //read from "parameters.dat". The character ">" should stay in line before parameters 
  param.open("parameters.dat", fstream::in);
  param.ignore(256, '>');
  param >> U >> mu >> beta >> nAverage;
  param.ignore(256, '>');
  param >> fluxNum >> fluxDen;
  param.ignore(256, '>');
  param >> wpoints >> wmin >> wmax;
  param.ignore(256, '>');
  param >> kxpoints >> kypoints;
  param.ignore(256, '>');
  param >> treshFermi >> treshSigma >> pade_order;
  param.close();
  deltaW = double(wmax - wmin)/double(wpoints - 1);
  //these two flags are to be read from terminal
  denFlag = false; chernFlag = false;

  


  //read pade coefficients from the file pade_coefficients.dat
  //can leave it in and then it will not do anything, but here I commented it out
  //if(U!=0){
  //  int order = 0;
  //  pade_num = new Complex[parameters.pade_order+1];
  //  pade_den = new Complex[parameters.pade_order+1];
  //  fstream pade_file("pade_coefficients.dat",fstream::in);
  //  double temp;
    //I think, reads in 
  //  while(!pade_file.eof()){
  //    pade_file>>temp;
  //    if(pade_file){
	//pade_den[order]=temp;
	//pade_file>>temp;
	//pade_den[order]+=I*temp;
	//pade_file>>temp; 
	//pade_num[order]=temp;
	//pade_file>>temp;
	//pade_num[order]+=I*temp;
	//order++;
  //    }
  //  }
  //}
}

// initialize all matrices with zero LxL  matrices
Matrices::Matrices(){  
   fill(energy,energy+L,0);
   fill(vx,vx+L*L,0);
   fill(vy,vy+L*L,0);    
}

// Matrices::~Matrices(){
//   delete[] vx;
//   delete[] vy;
//   delete[] energy;
// }
// MatricesE::MatricesE(){
//   vx = MatrixXcd::Zero(L,L);
//   vy = MatrixXcd::Zero(L,L);
//   energy = VectorXd::Zero(L);
// }

//parse command line arguments.
void getflags(int argc, char* argv[]){  
  for (int i = 0; i < argc; i++){
    string c = argv[i];
    if (c == "-den") parameters.denFlag = true;
    else if (c == "-chern") parameters.chernFlag = true;
    else{
      parameters.denFlag = false;
      parameters.chernFlag = false;
    }
  }
}

void printMatrix(Complex* matrix, int N, int M){
  for(int i=0;i < N; i++){
    for(int j=0;j < M;j++){
      cout << matrix[i*M+j]<<" ";
    }
    cout<< "\n";
  }
}


void printArray(Complex* array, int N){
  for(int i=0;i < N; i++){
      cout << array[i]<<"\n";
    }
}

void outDOS(double* dos){
  fstream dosFile("dos.dat",fstream::out);        //results/dos.dat
  double deltaW, omega;
  deltaW = double(parameters.wmax - parameters.wmin)/double(parameters.wpoints - 1);
  for(int w_count=0;w_count<parameters.wpoints;w_count++){
    omega =  parameters.wmin + w_count*deltaW;

    dosFile << omega <<"   "<<dos[w_count]<<"\n";
  }
  dosFile.close();
}
void precalculate_fermi_pref(double* fermi){ 
#pragma omp parallel
{ 
  double omega1,omega2,deltaW;
  #pragma omp for
  for( int w1_count=0; w1_count<parameters.wpoints;w1_count++){
    deltaW = double(parameters.wmax - parameters.wmin)/double(parameters.wpoints - 1);
    omega1 = parameters.wmin+deltaW*w1_count;
    //cout << "omega1:" << omega1 << " \n"; 
    //cout << "deltaW:" << deltaW << " \n"; 
    //cout << "parameters.deltaW:" << parameters.deltaW << " \n"; 
    for( int w2_count=0; w2_count<parameters.wpoints;w2_count++){
      omega2 = parameters.wmin+parameters.deltaW*w2_count;
      fermi[w1_count*parameters.wpoints+w2_count]=fermi_pref(parameters.beta,omega1,omega2);
    }
  }
}
}
