#define _USE_MATH_DEFINES
#include<fftw3.h>
#include<iostream>
#include<cmath>
#include <fstream>
#include <random>
//Macros for Real and imaginary parts
#define real 0
#define imag 1

double fun_1(double t)
{
  return pow(t,2)*cos(t);
}

double fun_2(double t)
{
  return sin(25*(2*M_PI)*t) + sin(80*(2*M_PI)*t) + sin(125*(2*M_PI)*t) + sin(240*(2*M_PI)*t) + sin(315*(2*M_PI)*t);  
}

void get_data(fftw_complex * x, int n, double * val, int problem)
{
 for(int i = 0;i<n;i++)
    {
      double temp = val[i];
      if(problem == 1)
      {
        temp = fun_1(temp);
      }
      if(problem == 2)
      {
        temp = fun_2(temp);
      }
      x[i][real] = temp;
      x[i][imag] = 0;
    }
}

void write_file(fftw_complex * y,fftw_complex * x,fftw_complex * xi,int n,const char S[],double * val)
{
    std::ofstream file;
    file.open(S,std::ios::out|std::ios::trunc);
      for(int i=n/2;i<n;i++)
      {   
          /*Normalize the result of IIFT*/
          xi[i][real] = xi[i][real]/n;
          xi[i][imag] = xi[i][imag]/n; 
          file<<val[i]<<"   "<<x[i][real]<<"   "<<y[i][real]<<"   "<<y[i][imag]<<"   "<<xi[i][real]<<"   "<<xi[i][imag]<<"\n";
      }
      file.close();
}


void write_file2(fftw_complex * y,fftw_complex * x,int n,const char S[],double * val)
{
    std::ofstream file;
    file.open(S,std::ios::out|std::ios::trunc);
   
      for(int i=0;i<n;i++)
      {
        file<<val[i]<<"   "<<x[i][real]<<"   "<<y[i][real]<<"   "<<y[i][imag]<<"\n";
      }
    file.close();
}

void get_grid(double * val, int n, double Period, int problem)
{
   if(problem == 1)
   {  
    for(int i = 0;i<n;i++)
      {
        double temp = -1*Period + (i*2*Period)/n;
        val[i] =  temp;
      }
   }
   if(problem == 2)
   {
    double temp = 1/(double)n; 
    for(int i = 0;i<n;i++)
      {
        val[i] =  i*temp;
      }
   } 
}


void add_normal_noise(fftw_complex * x,int n, double mean,double stddev)
{
    double noise; 
    std::default_random_engine generator;
    std::normal_distribution<double> normal(mean,stddev);
    for( int i = 0; i<n; i++){
        noise = normal(generator);
        x[i][real] = x[i][real] + noise;
    }
}


int main()
{
  //define list of different n
  int list_n[3]={16,32,64};
  int n  = 0,problem = 0;
  double Period =  0;


  for(int j=0;j<3;j++)
  {
      //Length of Complex Array 
      n = list_n[j];
      //Input array
      fftw_complex x[n]; // This is equivalent to: double x[n][2];
      //output array
      fftw_complex y[n];
      // To save the Inverse result
      fftw_complex xi[n];

      Period = M_PI;
      problem = 1;
      // fILL array with some data
      double val[n] = {0};
      //Linspacing between Period
      get_grid(val,n,Period,problem);
      //Value of Function at above linspacing
      get_data(x,n,val,problem);
      //Plant the FFT and execute it
      fftw_plan plan = fftw_plan_dft_1d(n,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(plan);
      //Do some cleaning
      fftw_destroy_plan(plan);
      
      // Inverse FFT
      fftw_plan plan_inverse = fftw_plan_dft_1d(n, y, xi, FFTW_BACKWARD, FFTW_ESTIMATE); 
      fftw_execute(plan_inverse);
      //Display the results and write to file
      
      if(j==0)
      {
        char const * name = "16.txt";
        write_file(y,x,xi,n,name,val);
      }
      if(j==1)
      {
        char const * name1 = "32.txt";
        write_file(y,x,xi,n,name1,val);
      }
      if(j==2)
      {
        char const * name2 = "64.txt";
        write_file(y,x,xi,n,name2,val);
      }

  }

//***********************Problem2,Part 2**************************************************
  n = pow(2,10);
  problem = 2;
  Period = 1;
  double val_2[n] = {0};
  fftw_complex x2[n];
  fftw_complex y2[n];
  get_grid(val_2,n,Period,problem);
  get_data(x2,n,val_2,problem);
  fftw_plan plan_2 = fftw_plan_dft_1d(n,x2,y2,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan_2);
  fftw_destroy_plan(plan_2);
  /*add result to file */
  char const * name4 = "problem_2_a.txt";
  write_file2(y2,x2,n,name4,val_2);


  //Let's add Noise
  double mean = 0, stddev = 2;
  add_normal_noise(x2,n,mean,stddev);
  //DFT
  fftw_plan plan_3 = fftw_plan_dft_1d(n,x2,y2,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan_3);
  fftw_destroy_plan(plan_3);

  /*add result to file */
  char const * name5 = "problem_2_b.txt";
  write_file2(y2,x2,n,name5,val_2);

  fftw_cleanup();
  return 0;
}
