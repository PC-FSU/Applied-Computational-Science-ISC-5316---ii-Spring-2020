#define _USE_MATH_DEFINES
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



void solve(std::vector<double> &residual,arma::mat & vec, int n)
{

    double tol = 1;
    int size = n+1;
    arma:: mat U_0(size,size,arma::fill::zeros);
    arma:: mat U_1(size,size,arma::fill::zeros); 
    
    
    U_0 = vec;
    while(tol > eps)
    {
        for(int i = 0;i<n;i++)
        {
            for(int j = 0;j<=n;j++) 
            {
                if(j == n || j == 0)
                {
                    U_1(i,j) = vec(i,j);
                    U_0(i,j) = vec(i,j);
                    //std::cout<<vec(i,j)<<"\n";
                }
                else
                {
                    if(i==0)
                    {   
                            //std::cout<<"1"<<"\n";
                            U_1(i,j) = (U_0(n-1,j) + U_0(i+1,j) + U_0(i,j-1) + U_0(i,j+1))/4;
                    }
                    else if(i == n-1)
                    {
                            //std::cout<<"2"<<"\n";
                            U_1(i,j) = (U_0(i-1,j) + U_0(0,j) + U_0(i,j-1) + U_0(i,j+1))/4;
                    }
                    else
                    {
                            //std::cout<<"3"<<"\n";
                            U_1(i,j) = (U_0(i-1,j) + U_0(i+1,j) + U_0(i,j-1) + U_0(i,j+1))/4;
                    }
                }
            }  
        }
        
       // U_0.print("u_0:");
       // U_1.print("u_1:");

        tol = arma::norm(U_1-U_0,"fro");
        U_0 = U_1;
       // std::cout<<tol<<"\n";
        residual.push_back(tol);
    }
    
    for(int i = 0;i<n;i++)
    {
        U_1(n,i) = U_1(0,i);
    }
    vec = U_1;
}


int main()
{
    int N[3]= {4,8,16};
    int n = 0;

    for(int i = 0;i<3;i++)
    {
        n = N[i];
        std::vector<double> residual;
        double dx = 0,dy=0;
        dx = 1/(double)n; dy = 1/(double)n;
        //set_boundary((double *)vec,n[0],dx,dy);
        arma::mat vec(n+1,n+1,arma::fill::zeros); 

        for(int i = 0;i<=n;i++)
        {
            vec(i,n) = pow(sin(M_PI*i*dx),2);
        }
        
        solve(residual,vec,n);
        vec.print("vec:");

        if(i==0)
        {
            char const * name = "4b.txt";
            write_file(residual,name);
            vec.save("4bres.txt", arma::raw_ascii);
        }
        else if (i == 1)
        {
            char const * name1 = "8b.txt";
            write_file(residual,name1);
            vec.save("8bres.txt", arma::raw_ascii);
        }
        else
        {
            char const * name2 = "16b.txt";
            write_file(residual,name2);
            vec.save("16bres.txt", arma::raw_ascii);
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