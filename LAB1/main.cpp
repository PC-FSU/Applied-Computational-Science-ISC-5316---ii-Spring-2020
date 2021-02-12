#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>
#include <omp.h>
#include <random>

// To generate n Random number and then map them to [0,1]
void get_number(double *res, int n)
{
 
 #pragma omp parallel 
 {
  static thread_local std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  #pragma omp for 
  for (int i=0; i<n; ++i) 
  {
    double number = distribution(generator);
    res[i] = number;
  }
 }
}

//********************************************************//
// To print the number
void print_num(double *res, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << res[i] << " ";
    }
    std::cout << "\n";
}

//********************************************************//
// To calculate chi_square
void chi_square(double *res, double *O, int n, int m)
{

    int e = 0;
    e = n / m;
    int i = 0;
    {
#pragma omp parallel for default(none) shared(n, O, res, m) private(i)
        // To calculate O
        for (i = 0; i < n; i++)
        {
            O[int(res[i] * m)] += 1;
        }
    }
    // To calculate chi_square
    double chi_square = 0;
    for (int i = 0; i < m; i++)
    {
        chi_square += pow((O[i] - e), 2);
    }
    chi_square = chi_square / e;
    std::cout << "The value of chi_square : " << chi_square << "\n";
}

//********************************************************//
// Function to write in a file.

void write_file(double *res, double *O, int n, int m)
{
    // The address of working directory of python code developed in first part.

    // To write the numbers.
    //const char *path = "C://Users//18503//Documents//ACS2//LAB//Number.txt";
    std::ofstream myfile("Number_100000.txt");
    {
        if (myfile.is_open())
        {
            for (int i = 0; i < n; i++)
            {
                if (i % 10 == 0)
                {
                    myfile << "\n";
                }
                myfile << res[i] << " ";
            }
            myfile.close();
        }
        else
        {
            std::cout << "Unable to open file."
                      << "\n";
        }
    }

    // To write the observed frequency in each subinterval
    //const char *path1 = "C://Users//18503//Documents//ACS2//LAB//Observe_frequency.txt";
    std::ofstream myfile1("Observe_frequency_100000.txt");
    {
        if (myfile1.is_open())
        {
            for (int i = 0; i < m; i++)
            {
                if (i % 10 == 0)
                {
                    myfile1 << "\n";
                }
                myfile1 << O[i] << " ";
            }
            myfile1 << n / m << " ";
            myfile1.close();
        }
        else
        {
            std::cout << "Unable to open file."
                      << "\n";
        }
    }
}

//********************************************************//
// The  main function
int main()
{
    int n = 0;
    std::cout << "Total Numbers you wish to Generate,N : "
              << "\n";
    std::cin >> n;
    double res[n] = {};

    double wtime = 0;
    wtime = omp_get_wtime();
    get_number(res, n);
    wtime = omp_get_wtime()-wtime;

    //print_num(res,n);
    int m = 0;
    for (;;)
    {
        std::cout << "Enter the number of subinterval,M (Note: M should divide N) : "
                  << "\n";
        std::cin >> m;
        if (n % m == 0)
        {
            break;
        }
    }
    double O[m] = {};
    double wtime1 = 0;
    wtime1 = omp_get_wtime();
    chi_square(res, O, n, m);
    wtime1 = omp_get_wtime()-wtime1;
    //write_file(res, O, n, m);
    std::cout<<"Time for Random Number Generation : "<<wtime<<" and time for counting the freuency : "<<wtime1<<", Total Time = "<<wtime+wtime1<<"\n";
}
