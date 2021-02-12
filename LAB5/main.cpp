#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include <iterator>


double gamma = 1.4;
double eps = 0.0000001;

void write_file(std::vector< std::vector<double> > &x,const char S[])
{
   std::ofstream file;
   file.open(S,std::ios::out|std::ios::trunc);
   for(const auto& vt : x) 
   {
     std::copy(vt.cbegin(), vt.cend(),std::ostream_iterator<double>(file, " "));
     file<<"\n";
   }
}


double C_s(double P,double rho)
{ 
    double C = 0;
    C =  sqrt(gamma*P*rho);
    std::cout<<"C : "<<C<<"\n";
    return C;
    
}

double W_s(double P_star,double P,double rho)
{
    double C = 0;
    C = C_s(P,rho);
    double W =0;
    //std::cout<<"c : "<<C<<"\n";
    W = C*sqrt(1 + ((gamma+1)/(2*gamma))*((P_star - P)/P));
    //std::cout<<"Square term : "<<1 + ((gamma+1)/(2*gamma))*((P_star - P)/P)<<"\n";
    return W;
}

double U_l(double U,double P_star,double P,double rho)
{
    double temp = 0;
    double W = 0;
    W = W_s(P_star,P,rho);
    temp = U - (P_star - P)/W;
    return temp;
}

double U_r(double U,double P_star,double P,double rho)
{
    double temp = 0;
    double W = 0;
    W = W_s(P_star,P,rho);
    temp = U + (P_star - P)/W;
    return temp;
}

double Q(double P_star,double P,double rho)
{
    double q = 0;
    double W = 0;
    double C = 0;
    W = W_s(P_star,P,rho);
    C = C_s(P,rho);
    q = 2*pow(W,3)/(pow(W,2) + pow(C,2));
    return q;
}

double Fun(double P,double rho_l,double u_l,double P_l,double rho_r,double u_r,double P_r)
{
    double temp = 0 ,W_l = 0,W_r = 0,res=0;
    W_l = W_s(P,P_l,rho_l);
    W_r = W_s(P,P_r,rho_r);
    temp = (W_l / (W_l + W_r)) * (P_r - P_l - W_r*(u_r - u_l));
    res = P-P_l-temp;
    return res;
}


double press(double P,double rho_l,double u_l,double P_l,double rho_r,double u_r,double P_r)
{
    double temp = 0 ,W_l = 0,W_r = 0,res=0;
    W_l = W_s(P,P_l,rho_l);
    W_r = W_s(P,P_r,rho_r);
    temp = (W_l / (W_l + W_r)) * (P_r - P_l - W_r*(u_r - u_l));
    res = P + temp;
    return res;
}



void Bisection(double rho_l,double u_l,double P_l,double rho_r,double u_r,double P_r,std::vector<double> &vec,std::vector<double> &vec2)
{
    double res =0,r_k = 1,a=0,b=0,c=0,F_a = 0,F_b = 0,F_c =0,temp=0 ;
    a = 0;
    //b = 2000;
    if(u_l == 0 || u_r == 0)
    {
        b =( P_l + P_r )/ 2;
    }
    else{
        b = ( u_l*P_l + u_r*P_r )/ 2;
    }
    
    //To check Initial Guess
    F_a = Fun(a,rho_l,u_l,P_l,rho_r,u_r,P_r);
    F_b = Fun(b,rho_l,u_l,P_l,rho_r,u_r,P_r);
    if(F_a > 0 && F_b < 0 || F_a < 0 && F_b > 0)
       {
            // MAIN loop
            while(r_k > eps)
            {
                F_a = Fun(a,rho_l,u_l,P_l,rho_r,u_r,P_r);
                F_b = Fun(b,rho_l,u_l,P_l,rho_r,u_r,P_r);
                c = (a + b)/2;
                F_c = Fun(c,rho_l,u_l,P_l,rho_r,u_r,P_r);
                
                if(F_a <0 && F_c > 0)
                {
                b = c;
                }
                else if(F_c < 0 && F_b>0)
                {
                    a = c;
                }
                r_k  = fabs(c - temp)/fabs(c);
                temp = c;
                vec.push_back(r_k);
                vec2.push_back(c-F_c);
            }
        std::cout<<"Bisection Pressure :"<<(c-F_c)<<"\n";   
        }
    else{
        std::cout<<" Wrong Initial Guess "<<"\n";
    }
}


void newton(double rho_l,double u_l,double P_l,double rho_r,double u_r,double P_r,std::vector<double> &vec, std::vector<double> &vec2)
{
    double res = 0;
    double P_k = 0,P_k_1 = 0,r_k=1,Q_r=0,Q_l=0,U_star_l=0,U_star_r=0;
    while(r_k > eps)
    {
        Q_r =  Q(P_k,P_r,rho_r);
        Q_l =  Q(P_k,P_l,rho_l);
        std::cout<<"Q_r : "<<Q_r<<"\n";
        std::cout<<"Q_l : "<<Q_l<<"\n";
        U_star_l = U_l(u_l,P_k,P_l,rho_l);
        U_star_r = U_r(u_r,P_k,P_r,rho_r);
        P_k_1 = P_k - ((Q_r*Q_l)/(Q_r + Q_l))*(U_star_r - U_star_l); 
        r_k = fabs(P_k_1 - P_k)/fabs(P_k_1);
        vec.push_back(r_k);
        P_k = P_k_1;
        vec2.push_back(P_k_1);
        
        //std::cout<<r_k;
    }
    res = P_k_1;
    std::cout<<"Newton Pressure :"<<P_k_1<<"\n";
   
}


void secent(double rho_l,double u_l,double P_l,double rho_r,double u_r,double P_r,std::vector<double> &vec,std::vector<double> &vec2 )
{
double P_k = 0,P_k_1 = 0,P_k_2 = 0,r_k = 1,F_k=0,F_k_1=0;
P_k = 20;
P_k_1 = 2000;

    while(r_k > eps)
        {
            F_k = Fun(P_k,rho_l,u_l,P_l,rho_r,u_r,P_r);
            F_k_1 = Fun(P_k_1,rho_l,u_l,P_l,rho_r,u_r,P_r);
            P_k_2 = P_k_1 - F_k_1*((P_k_1 - P_k)/(F_k_1 - F_k));
            r_k = fabs(P_k_2 - P_k_1)/fabs(P_k_2);
            P_k = P_k_1;
            P_k_1 = P_k_2;
            vec.push_back(r_k);
            vec2.push_back(P_k_2);
        }
    std::cout<<"Secant Pressure :"<<P_k_2<<"\n";
}



int main()
{
    double rho_l[4] = {1,1,1,5.99924};
    double u_l[4] = {0,0,0,19.5975};
    double P_l[4] = {1,1000,0.01,460.894};
    double rho_r[4] = {0.125,1.0,1.0,5.99242};
    double u_r[4] = {0,0,0,-6.19633};
    double P_r[4] = {0.1,0.01,100,46.0950};
/*
//************************************Bisection***************************************************
    std::vector< std::vector<double> > residual_bisection;
    std::vector< std::vector<double> > pressure_bisection;
    
//Calculate the residual and Pressure
    for(int i = 0 ;i<4;i++)
    {
        std::vector<double> vec;
        std::vector<double> vec2;
        Bisection(rho_l[i],u_l[i],P_l[i],rho_r[i],u_r[i],P_r[i],vec,vec2);
        residual_bisection.push_back(vec);
        pressure_bisection.push_back(vec2);
        //std::cout<<presure_bisection[i]<<"\n";
    }
//Write to a file
    char const * name = "residual_bisection.txt";
    char const * name_1 = "pressure_bisection.txt";
    write_file(residual_bisection,name);
    write_file(pressure_bisection,name_1);
*/

//*****************************************Newton***********************************************
    std::vector< std::vector<double> > residual_newton;
    std::vector< std::vector<double> > pressure_newton;
    

//Calculate the residual and Pressure    
    for(int i = 0 ;i<1;i++)
    {
        std::vector<double> vec;
        std::vector<double> vec2;
        newton(rho_l[i],u_l[i],P_l[i],rho_r[i],u_r[i],P_r[i],vec,vec2);
        residual_newton.push_back(vec);
        pressure_newton.push_back(vec2);
        //std::cout<<presure_newton[i]<<"\n";
    }
   
//Write to a file
    char const * name1 = "residual_newton.txt";
    char const * name1_1 = "pressure_newton.txt";

    write_file(residual_newton,name1);
    write_file(pressure_newton,name1_1);

/*
//*****************************************Secent*********************************************
    std::vector< std::vector<double> > residual_secent;
    std::vector< std::vector<double> > pressure_secent;
//Calculate the residual and Pressure       
    for(int i = 0 ;i<4;i++)
    {
        std::vector<double> vec;
        std::vector<double> vec2;
        secent(rho_l[i],u_l[i],P_l[i],rho_r[i],u_r[i],P_r[i],vec,vec2);
        residual_secent.push_back(vec);
        pressure_secent.push_back(vec2);
    }
// Write to a file

    char const * name2 = "residual_secent.txt";
    char const * name2_1 = "pressure_secent.txt";

    write_file(residual_secent,name2);
    write_file(pressure_secent,name2_1);
   */ 
}

































/*

//Write to a file
   std::ofstream output_file("secent_residual.txt");
   for(const auto& vt : residual_secent) 
   {
     std::copy(vt.cbegin(), vt.cend(),std::ostream_iterator<double>(output_file, " "));
     output_file<<"\n";
   }
   */

  // To print out the residual
    /*
    for ( const std::vector<double> &v : residual_bisection )
    {
        for ( double x : v ) std::cout << x << ' ';
        std::cout << std::endl;
    }
    
*/
