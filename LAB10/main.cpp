#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<armadillo>


double eps = 0.00001;

void write_file(std::vector<double> &x,const char S[])
{
    std::ofstream file;
    file.open(S,std::ios::out|std::ios::trunc);
    int n=x.size();
    for(int i=0;i<n;i++)
    {
        file<<x[i]<<"\n";
    }
}

double f(double x, double y)
{
    return -(20*pow(x,3)*pow(y,3) + 6*pow(x,5)*y);
}

double g(double x , double y)
{
    if(x == 0 || y == 0)
    {
        return 0;
    }
    else if(y==1)
    {
        return pow(x,5);
    }
    else if(x == 1)
    {
        return pow(y,3);
    }
}

double u(double x, double y)
{
    return pow(x,5)*pow(y,3);
}


void true_result(arma::mat & TrUe,int n,double dx, double dy)
{
    for(int i = 0;i<=n;i++)
    {
        for(int j = 0;j<=n;j++) 
        {
             TrUe(i,j) = pow(i*dx,5)*pow(j*dy,3);
        }  
    }
}


double solve(std::vector<double> &residual,arma::mat & vec,arma::mat & TrUe, int n, double dx, double dy)
{

    double tol = 1;
    int size = n+1;
    double h = 1/(double)n;
    double error = 0;
    arma:: mat U_0(size,size,arma::fill::zeros);
    arma:: mat U_1(size,size,arma::fill::zeros); 

    while(tol > eps)
    {
        for(int i = 0;i<=n;i++)
        {
            for(int j = 0;j<=n;j++) 
            {
                if(i == 0 || i == n || j == 0 || j == n)
                {
                    U_1(i,j) = vec(i,j);
                    U_0(i,j) = vec(i,j);
                    //std::cout<<vec(i,j)<<"\n";
                }
                else
                {
                    U_1(i,j) = (U_0(i-1,j) + U_0(i+1,j) + U_0(i,j-1) + U_0(i,j+1) + pow(h,2)*vec(i,j))/4;
                }
            }  
        }
        tol = arma::norm(U_1-U_0,"fro");
        U_0 = U_1;
        residual.push_back(tol);
       // std::cout<<tol<<"\n";
    }
   // U_1.print("u_1:");
    error = arma::norm(U_1-TrUe,"fro");
    return error*h;

}


int main()
{
    int N[3]= {4,8,16};
    int n = 0;

    for(int i = 0;i<3;i++)
    { 
        n = N[i];
        std::vector<double> residual;
        double error = 0;
        double dx = 0,dy=0;
        dx = 1/(double)n; dy = 1/(double)n;
        //set_boundary((double *)vec,n[0],dx,dy);
        arma::mat vec(n+1,n+1,arma::fill::zeros); 

        for(int i = 0;i<=n;i++)
        {
            for(int j=0;j<=n;j++)
            {
                if(i == n && j>0)
                {   
                    vec(i,j) = pow(j*dy,3);
                }
                else if(j == n && i>0)
                {
                    vec(i,j) = pow(i*dx,5);
                }
                else if(i == 0 || j == 0)
                {
                    vec(i,j) = 0;
                }
                else
                {
                    double Fun_val  = 0;
                    Fun_val = f(i*dx,j*dy);
                    vec(i,j) = Fun_val;
                }
            }
        }
        //vec.print("vec:");
        arma::mat TrUe(n+1,n+1,arma::fill::zeros);
        true_result(TrUe,n,dx,dy);
        //TrUe.print("true :");
        error = solve(residual,vec,TrUe,n,dx,dy);
        std::cout<<error<<"\n";
        
        if(i==0)
        {
            char const * name = "4.txt";
            write_file(residual,name);
        }
        else if (i == 1)
        {
            char const * name1 = "8.txt";
            write_file(residual,name1);
        }
        else
        {
            char const * name2 = "16.txt";
            write_file(residual,name2);
        }
    }
}











































/*


    for ( const std::vector<double> &v : vec )
    {
        for ( double x : v ) std::cout << x <<' ';
        std::cout << std::endl;
    }



void set_boundary(double * vec, int n, double dx, double dy)
{
    int index = 0;
    for(int i = 0;i<=n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            
            if(i == n && j>0)
            {   
                //std::cout<<i*n + index<<"\n";
                //std::cout<<" AT I : "<<i<<" and j : "<<j<<" "<<pow(j*dy,3)<<"  loop : "<<1<<"\n";
                vec[i*n + j+i] = pow(j*dy,3);
            }
            else if(j == n && i>0)
            {
                //std::cout<<i*n + index<<"\n";
                //std::cout<<" AT I : "<<i<<" and j : "<<j<<" "<<pow(i*dx,5)<<"  loop : "<<2<<"\n";
                vec[i*n + i+j] = pow(i*dx,5);
            }
            else if(i == 0 || j == 0)
            {
               // std::cout<<i*n + index<<"\n";
                //std::cout<<" AT I : "<<i<<" and j : "<<j<<" "<<0<<"  loop : "<<3<<"\n";
                vec[i*n + i+j] = 0;
            }
            else
            {
                //std::cout<<i*n + index<<"\n";
                double Fun_val  = 0;
                Fun_val = f(i*dx,j*dy);
                //std::cout<<" AT I : "<<i<<" and j : "<<j<<" "<<Fun_val<<"  loop : "<<4<<"\n";
                vec[i*n + i+j] = Fun_val;
            }
        }   
    }

}

*/