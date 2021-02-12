#define _USE_MATH_DEFINES
#include<iostream>
#include<armadillo>
#include<cmath>
#include<math.h>

double f_prime(double x)
{
    return -10*x;
}

double Newton_Solver(double y_n,double h)
{
   double eps = 0.00001;
   double x_1 = 0, x_0 = y_n;
   while(1)
   {
       x_1 = x_0 - (x_0 + h*10*x_0 - y_n)/(1+10*h);
       if(x_1-x_0<eps)
       {
           break;
       }
       x_1 =  x_0;
   }
return x_1;
}


void exact(arma::mat & result,double h)
{
    int n = result.n_cols;
    for(int i=0;i<n-1;i++)
    {
       result(0,i) = exp(-10*(i*h));
    }
    result(0,n-1) = exp(-10);
}

void Forward_Euler(arma::mat & result,double h)
{
    double n = result.n_cols;
    result(1,0) = 1;
    double f_p = 0;
    for(int i=0;i<n-1;i++)
    {
     f_p = f_prime(result(1,i));   
     result(1,i+1) = result(1,i) + h*f_p;
    }
}


void Backward_Euler(arma::mat & result,double h)
{
    double n = result.n_cols;
    result(2,0) = 1;
    double f_p = 0;
    for(int i=0;i<n-1;i++)
    {   
     result(2,i+1) = Newton_Solver(result(2,i),h);
    }
}

void Forward_Euler_2(int Nt,int Ns,double dt,arma::mat & Solution)
{
    double dx  = 1/((double)Ns),Nu = 0.002,ul=0,ur=0;
    arma::mat    u0(1,Ns,arma::fill::zeros);
    arma::mat  unew(1,Nt,arma::fill::zeros);
    arma::mat     x(1,Ns,arma::fill::zeros);
    //Calculate x
    for(int i=0;i<Ns;i++)
    {
        x(0,i) = i*dx;
    }
    //x.print();
    //Updata Value of X with SIN(2*PI*X)
    u0 = arma::sin(2*M_PI*x);
    
    unew = u0;
    
    for(int i=0;i<Nt;i++)
    {
        for(int j=0;j<Ns;j++)
        {
            // Find left and right pair
            if(j==0){
                ul = u0(0,Ns-1);   // For First Element
                ur = u0(0,j+1);
            }
            else if(j == Ns-1){    // For Last Element
                ul = u0(0,j-1);
                ur = u0(0,0);
            }
            else{
                ul = u0(0,j-1);   
                ur = u0(0,j+1);
            }
            // Updatae Unew
            unew(0,j) = u0(0,j) + dt*(-sin(2*M_PI*x(0,j))*(ur-ul)/(2*dx) + Nu*(ur -2*u0(0,j) + ul)/(dx*dx));
        }
        u0 = unew;
    }
    Solution = u0;
}


void Backward_Euler_2(int Nt,int Ns,double dt,arma::mat & Solution)
{
    double dx  = 1/((double)Ns),Nu = 0.002;

    arma::mat    u0(Ns,1,arma::fill::zeros);
    arma::mat     x(Ns,1,arma::fill::zeros);
    arma::mat     A(Ns,Ns,arma::fill::zeros);
    arma::mat     b(Ns,1,arma::fill::zeros);

    //Calculate x
    for(int i=0;i<Ns;i++)
    {
        x(i,0) = i*dx;
    }
    //Updata Value of X with SIN(2*PI*X)
    u0 = arma::sin(2*M_PI*x);


    for(int i=0;i<Nt;i++)
    {
        for(int j=0;j<Ns;j++)
        {
            
            if(j==0){              // For First Element
                A(j,j+1)   = -dt*(-0.5*sin(2*M_PI*x(j,0))/dx + Nu/pow(dx,2))/(1+2*Nu*dt/pow(dx,2));
                A(j,j)     =  1;
                A(j,Ns-1)  = -dt*(0.5*sin(2*M_PI*x(j,0))/dx + Nu/pow(dx,2))/(1+2*Nu*dt/pow(dx,2));
            }
            else if(j == Ns-1){    // For Last Element
                A(j,0)     = -dt*( -0.5*sin(2*M_PI*x(j,0))/dx + Nu/pow(dx,2) )/ (1+2*Nu*dt/pow(dx,2)) ;
                A(j,j)     =  1;
                A(j,j-1)  = -dt*( 0.5*sin(2*M_PI*x(j,0))/dx + Nu/pow(dx,2)  )/ (1+2*Nu*dt/pow(dx,2)) ;       
            }
            else{
                A(j,j+1)   = -dt*(-0.5*sin(2*M_PI*x(j,0))/dx + Nu/pow(dx,2))/(1+2*Nu*dt/pow(dx,2));
                A(j,j)     =  1;
                A(j,j-1)   = -dt*(0.5*sin(2*M_PI*x(j,0))/dx + Nu/pow(dx,2))/(1+2*Nu*dt/pow(dx,2));   
            }
            
            b(j,0) = u0(j,0)/(1+2*Nu*dt/pow(dx,2));
        }
        //Linear Solver
        u0 = arma::solve(A,b);
    }
    Solution = u0;
}


void CNLF(int Nt,int Ns,double dt,arma::mat & Solution)
{
    double dx  = 1/((double)Ns),Nu = 0.002;

    arma::mat    uNminusOne(Ns,1,arma::fill::zeros);
    arma::mat    uN(Ns,1,arma::fill::zeros);
    arma::mat    uNplusOne(Ns,1,arma::fill::zeros);
    arma::mat    x(Ns,1,arma::fill::zeros);
    arma::mat    A(Ns,Ns,arma::fill::zeros);
    arma::mat    b(Ns,1,arma::fill::zeros);

    //Calculate x
    for(int i=0;i<Ns;i++)
    {
        x(i,0) = i*dx;
    }
    //Initialize
    uNminusOne = arma::sin(2*M_PI*x);
    //Get First Solution from BE
    Backward_Euler_2(1,Ns,dt,uN);
    //Define Fraction term for matrix A and B
    double FractionTerm = 1 + Nu*dt/pow(dx,2);

    for(int i=0;i<Nt;i++)
    {
        for(int j=0;j<Ns;j++)
        {
            if(j==0){  // For First Element

                A(j,j+1)  = -dt*0.5*Nu/pow(dx,2)/FractionTerm;
                A(j,j)    = 1;
                A(j,Ns-1) = -dt*0.5*Nu/pow(dx,2)/FractionTerm;
                b(j,0)    = ( uNminusOne(j,0) - dt*sin(2*M_PI*x(j,0))*(uN(j+1,0) - uN(Ns-1,0))/dx + 0.5*Nu*dt*(uNminusOne(j+1,0) - 2*uNminusOne(j,0) + uNminusOne(Ns-1,0))/(dx*dx) ) / FractionTerm;
            }
            else if(j==Ns-1){  // For Last Element
  
                A(j,0)    = -dt*0.5*Nu/pow(dx,2)/FractionTerm;
                A(j,j)    = 1;
                A(j,j-1)  = -dt*0.5*Nu/pow(dx,2)/FractionTerm;
                b(j,0)    = ( uNminusOne(j,0) - dt*sin(2*M_PI*x(j,0))*(uN(0,0) - uN(j-1,0))/dx + 0.5*Nu*dt*(uNminusOne(0,0) - 2*uNminusOne(j,0) + uNminusOne(j-1,0))/(dx*dx) ) / FractionTerm;
            }
            else{
                A(j,j+1)  = -dt*0.5*Nu/pow(dx,2)/FractionTerm;
                A(j,j)    = 1;
                A(j,j-1)  = -dt*0.5*Nu/pow(dx,2)/FractionTerm;
                b(j,0)    = ( uNminusOne(j,0) - dt*sin(2*M_PI*x(j,0))*(uN(j+1,0) - uN(j-1,0))/dx + 0.5*Nu*dt*(uNminusOne(j+1,0) - 2*uNminusOne(j,0) + uNminusOne(j-1,0))/(dx*dx) ) / FractionTerm;
            }
        }
        //Linear Solver
        uNplusOne = arma::solve(A,b);
        // Reassign for next iterations
        uNminusOne = uN;
        uN = uNplusOne;
    }
    //Return the solution
    Solution = uNplusOne;
}



int main()
{
    double h[4] = {0.21,0.2,0.05,0.001};
    int dimension[4] = {0};
    
    for(int i = 0;i<4;i++)
    {
    dimension[i] = std::round(1/h[i]) + 1;
    }
    
    /*Each Matrix Contain 3 rows, 1st row for exact solution, 2nd for Forward_Euler solution and 3rd for Backward Euler */
    //           Defien Matrix                                 get Exact Result            Get FE result                       Get BE result                          Save Results
    arma::mat Result_h_1(3,dimension[0],arma::fill::zeros);  exact(Result_h_1,h[0]);  Forward_Euler(Result_h_1,h[0]); Backward_Euler(Result_h_1,h[0]); Result_h_1.save("FE_h1.txt",arma::raw_ascii);
    arma::mat Result_h_2(3,dimension[1],arma::fill::zeros);  exact(Result_h_2,h[1]);  Forward_Euler(Result_h_2,h[1]); Backward_Euler(Result_h_2,h[1]); Result_h_2.save("FE_h2.txt",arma::raw_ascii);
    arma::mat Result_h_3(3,dimension[2],arma::fill::zeros);  exact(Result_h_3,h[2]);  Forward_Euler(Result_h_3,h[2]); Backward_Euler(Result_h_3,h[2]); Result_h_3.save("FE_h3.txt",arma::raw_ascii);
    arma::mat Result_h_4(3,dimension[3],arma::fill::zeros);  exact(Result_h_4,h[3]);  Forward_Euler(Result_h_4,h[3]); Backward_Euler(Result_h_4,h[3]); Result_h_4.save("FE_h4.txt",arma::raw_ascii);

    //********************************   Part 2 ************************
    //******* Forward Euler ************
    //Nt = Number of discrete points in time domain
    //Ns = Number of discrete points in space domain
    int Nt = 100,Ns = 64;
    double dt = 0;
    dt = 1/(double)Nt;

    //Forward_Eular_2
    arma::mat CalculatedSolutionForwardEular(Ns,1,arma::fill::zeros);
    Forward_Euler_2(Nt,Ns,dt,CalculatedSolutionForwardEular);
    CalculatedSolutionForwardEular.save("FE_2.txt",arma::raw_ascii);

    //Backward_Eular_2
    arma::mat CalculatedSolutionBackwardEular(Ns,1,arma::fill::zeros);
    Backward_Euler_2(Nt,Ns,dt,CalculatedSolutionBackwardEular);
    CalculatedSolutionBackwardEular.save("BE_2.txt",arma::raw_ascii);
    
    //CNLF METHOD
    arma::mat CalculatedSolutionCNLP(Ns,1,arma::fill::zeros);
    CNLF(Nt,Ns,dt,CalculatedSolutionCNLP);
    CalculatedSolutionCNLP.save("CNLF.txt",arma::raw_ascii);

}