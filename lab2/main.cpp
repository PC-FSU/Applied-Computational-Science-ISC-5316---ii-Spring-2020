#include<iostream>
#include<fstream>
#include<math.h>
#include<cmath>
#include <random>
#define epsilon  0.0001
#define lambda  1
#define g_y 0.1
// theoretical result
#define distribution_mean 1.253
#define  distribution_variance 0.42920367320510344



// To read the file
void read_file(double * list,int n, const char S[])
{
 std::ifstream file;
 file.open(S,std::ios::in);
    if(file.is_open())
    {
        double temp; 
        int i = 0;
        while(file>>temp)
        {
            list[i] = temp;
            i++;
        }
        file.close();
    }
    else
    {
        std::cout<<"unable to open the file"<<"\n";
    }

}

// To write to the file
void write_file(double * list,int n, const char  S[])
{
 std::ofstream file;
 file.open(S,std::ios::out);
    if(file.is_open())
    {
        double temp; 
        int i = 0;
        for(int i=0;i<n;i++)
        {
            file<<list[i]<<" ";
        }
        file.close();
    }
    else
    {
        std::cout<<"unable to open the file"<<"\n";
    }

}


// To print a list of Number
void print_num(double *list, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << list[i] << " ";
    }
    std::cout << "\n";
}


double mean(double * list, int n)
{ 
    double  sum = 0;
    for(int i =0; i<n;i++)
    {
        sum += list[i];
    }
 return sum/n;
}

double variance(double * list, int n)
{
    double MEAN = mean(list,n);
    double var = 0;
    for(int i =0; i<n;i++)
    {
        var += pow(list[i] - MEAN,2);
    }
   return var/(n-1);
}


//**************************************************Supportive function for problem 1********************************************

// CDF for problem 1.
double CDF(double X)
{
    return 1 - exp(-(lambda*X));
}

// CDF derivative or PDF for problem 1.
double PDF(double X)
{
    
    return lambda*exp(-lambda*X);
}


//*P_CDF is pointer to function CDF.
//*P_PDF is pointer to function PDF.
// Inverse transform sampling routine for problem 1.
void inverse_transform_sample(double (*P_CDF)(double X),double (*P_PDF)(double U), double * number_list, double * inverse_random_list,int n)
{
 
    for(int i = 0;i<n;i++)
    {
        // Newton Raphson Routine
        double X = 0;
        double U = number_list[i];
        double h = ((*P_CDF)(X) - U) / (*P_PDF)(X);
        double prev = -3;
        double curr = X;
        for(;fabs(curr - prev) > epsilon;)
            {
                prev  = X;
                h = ((*P_CDF)(X) - U) / (*P_PDF)(X);
                X = X-h; 
                curr = X;
            }
        inverse_random_list[i] = X;
    }

}


//**************************************************Supportive function for problem 2********************************************


double P( double X)
{
    return X*exp(-((std::pow(X,2))/2));
}


double G()
{
  static thread_local std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  
  double number = distribution(generator);
  return number;
}


double PDF_Piecwise_G( double X)
{ 
    if(X <= 5 && X > 0)
    {
        return 0.19;
    }
    else if( X > 5 && X <= 10 )
    { 
        return 0.01;
    }
}

double CDF_Piecwise_G(double X)
{
    if(X <= 5 && X > 0)
    {
        return 0.19*X;
    }
    else if( X > 5 && X <= 10 )
    {
        return 0.9 + 0.01*X;
    }
}

double inverse_transform_sample_2(double (*CDF)(double X),double (*PDF)(double U), double U )
{
        // Newton Raphson Routine
        double X = 1;
        double h = ((*CDF)(X) - U)/ (*PDF)(X);
        double prev = -3;
        double curr = X;
        for(;fabs(curr - prev) > epsilon;)
            {
                prev  = X;
                h = ((*CDF)(X) -U) / (*PDF)(X);
                X = X-h; 
                curr = X;
            }
        return X;
}


void uniform_sampling(double M,double * number_list, int number_of_sample,bool disp)
{
    int total_count = 0;
    int accept_count = 0;
    double u = 0;
    double y = 0; 
    // For uniform Distribution
    while(accept_count < number_of_sample )
    {
      y = G()*10;
      u = G();
      if( P(y)/(M*g_y) >= u)
      {
         number_list[accept_count] = y;
         accept_count++;
      }
      total_count++;
    }
    if(disp)
    {
    std::cout<<"probability of acceptence,Uniform (Problem 2,part 1): " << double(accept_count)/double(total_count)<<"\n"<<"\n";
    }
}


void Piecewise_sampling(double M,double * number_list, int number_of_sample)
{
    int total_count = 0;
    int accept_count = 0;
    double u = 0;
    double y = 0; 
    while(accept_count < number_of_sample )
    {
      y = G();
      y = inverse_transform_sample_2(CDF_Piecwise_G,PDF_Piecwise_G,y);
      u =  G();
      if( P(y)/(M*PDF_Piecwise_G(y)) >= u)
      {
        number_list[accept_count] = y;
        accept_count++;
      }
       total_count++;
    } 
    std::cout<<"probability of acceptence,Piecewise(Problem 2,part 2) : " << double(accept_count)/double(total_count)<<"\n"<<"\n";
}



//**************************************************Supportive function for problem 3********************************************


void CMT(double * Mean, double * Stats, int number_of_sample, int number_of_iteration)
{
    double M = 6.5;
    double temp = 0;
    double temp_list[number_of_sample] = {};
    bool disp = false;
    for(int i = 0;i<number_of_iteration;i++)
    {
         uniform_sampling(M,temp_list,number_of_sample,disp);
         Mean[i] = mean(temp_list,number_of_sample);
         Stats[i] = sqrt(number_of_sample)*(distribution_mean - Mean[i]);
    }
}


// Main Function
int main()
{
    // **************************************Problem_1***************************************************************
    double number_list[1000] = {};
    double inverse_random_list[1000] = {};
    // Read 1000 number Generated using Python
    char const * s1 = "Number_list_1.txt";
    read_file(number_list,1000,s1);
    // Inverse Sample
    inverse_transform_sample(CDF,PDF,number_list,inverse_random_list,1000); 
    // Write the transformed  Number to file
    char const * s2 = "Inverse_random_sample.txt";
    write_file(inverse_random_list,1000,s2);
    
    // *************************************Problem 2****************************************************************
    // **********************First Part*************************
  
    double number_list_2[10000] = {};
    int number_of_sample = 10000;
    double M = 6;
    bool disp = true;
    // Generate from uniform distribution
    uniform_sampling(M,number_list_2,number_of_sample,true);
    // Write to file
    char const * s3 = "Number_list_2.txt";
    write_file(number_list_2,10000,s3);
    double mean2_1 = mean(number_list_2,10000);
    double var2_1  = variance(number_list_2,10000);
    std::cout<<"Theoretical mean and variance : "<<distribution_mean<<", "<<distribution_variance<<"\n"<<"\n";
    std::cout<<"Calculated mean and variance  : "<<mean2_1<<", "<<var2_1<<"\n"<<"\n";
    
    // **************Second Part**************************
    // For Piecewise
    double number_list_3[10000] = {};
    M = 3.2;
    // Generate from piecewise distribution
    Piecewise_sampling(M,number_list_3,number_of_sample);
    char const * s4 = "Number_list_3.txt";
    write_file(number_list_3,10000,s4);
    double mean2_2 = mean(number_list_3,10000);
    double var2_2  = variance(number_list_3,10000);
    std::cout<<"Theoretical mean and variance : "<<distribution_mean<<", "<<distribution_variance<<"\n"<<"\n";
    std::cout<<"Calculated mean and variance  : "<<mean2_2<<", "<<var2_2<<"\n"<<"\n";

    //*************************************Problem 3****************************************************************
    double mean_1[1000] = {};
    double mean_2[1000] = {};
    double stat_1[1000] = {};
    double stat_2[1000] = {};
    CMT(mean_1,stat_1,100,1000);
    CMT(mean_2,stat_2,1000,1000);
    char const * s5 = "stats_sample_100.txt";
    write_file(stat_1,1000,s5);
    char const * s6 = "stats_sample_1000.txt";
    write_file(stat_2,1000,s6);
}






