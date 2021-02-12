#include<armadillo>
#include<iostream>
using namespace arma;

void hs(vec &  x)
{
    std::cout<<x;
    vec x1 = x;
    x1 = x1+x1;
    std::cout<<x1;
    x = x1;
}

int main()
{
    //mat x_0(2,1);
   // x_0(1,1) =1;
   // x_0(2,1) = 10;
   // vec x_1 = {2,3};
   // x_1 =  x_1.t()*x_0;
    //
   // hs(x_1);
    mat x_0=mat({{3.5,4.0},{1,21}}).t();
    x_0(0,0) = 2;
    std::cout<<x_0;
}