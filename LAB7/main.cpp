#include<iostream>
#include<armadillo>
#include<vector>
#include<math.h>
#include<omp.h>



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


double get_fun(arma::mat & x)
{
    double temp = 0,x_1 = 0,x_2 =0;
    x_1 = x(0,0), x_2  = x(1,0);
    temp =  100 * pow(x_2 - pow(x_1,2),2) + pow(1-x_1,2);
    return temp;
}

void get_grad(arma::mat & res, arma::mat & x)
{
    double x_1 = 0,x_2 =0;
    x_1 = x(0,0), x_2  = x(1,0);
    res(0,0) = 400*(pow(x_1,2) - x_2)*x_1 + 2*x_1 - 2;
    res(1,0) = 200*(x_2 - pow(x_1,2));
}



void get_hessian(arma::mat & H, arma::mat & x)
{
    double factor = 0,temp = 0; 
    //Hessian  
    H(0,0) = 400*(3*pow(x(0,0),2) - x(1,0)) + 2;
    H(0,1) = -1*400*x(0,0);
    H(1,0) = H(0,1);
    H(1,1) = 200;

    //Hessian Inv
    factor = H(0,0)*H(1,1) - H(1,0)*H(0,1);
    H(0,1) = -H(0,1);
    H(1,0) = -H(1,0);
    
    temp = H(0,0);
    H(0,0) = H(1,1);
    H(1,1) =  temp;
    H = H/factor;

}


void get_p(arma::mat & res,arma::mat & x,arma::mat & grad,double method)
{    
    arma::mat temp = arma::mat({0,0}).t();  
    if (method == 0)
    {   
        temp = -1*grad;
    }
    
    else 
    { 
        arma::mat B =  arma::mat({{0,0},{0,0}});
        get_hessian(B,x); 
        temp = B*grad;
        temp = -temp;
    }
    
res = temp;
}


void get_p_BFGS(arma::mat & p,arma::mat & B_k_1,arma::mat & B_k,arma::mat & x_k_1,arma::mat & x_k,arma::mat & grad_k_1,arma::mat & grad_k, double count)
{
    double alpha=0,beta=0,factor=0;
    arma::mat S = arma::mat({0,0}).t();
    arma::mat U = arma::mat({0,0}).t();
    arma::mat V = arma::mat({0,0}).t();
    arma::mat Prod;

      if(count  == 0)
      {

          get_hessian(B_k_1,x_k_1);
          B_k_1 = B_k;  
                
      }
      else
      {
          if(count < 3)
            {
                get_hessian(B_k_1,x_k_1);
                //std::cout<<B_k_1<<"\n";
                
            }
          S =  x_k_1 - x_k;
          U = grad_k_1 - grad_k;
          V = B_k*S;
          
          // To calculate alpha
          Prod = U.t()*S;
          alpha = 1/Prod(0,0);

          // To calculate beta 
          Prod = S.t()*V;
          beta = -1/Prod(0,0);
          
          // To calculate the denominator value in Morrison's formula
          Prod = S.t()*U;
          factor = Prod(0,0);
         
          // To calculate the first paranthesis of second term in Morrison's formula
          Prod =  S.t()*U + U.t()*(B_k*U);

          B_k_1 = B_k +  Prod(0,0)*(S*S.t())/ pow(factor,2) -  ( B_k*(U*S.t()) + (S*U.t())*B_k ) / factor;

      }   

      if(count < 3)
            {
                std::cout<<B_k_1<<"\n";   
            }
  
      p = -1*B_k_1*grad_k_1;
}


double get_alpha(arma::mat & x,arma::mat & grad,arma::mat & p)
{
    double alpha = 1, rho = 0.1 , gamma_1 = 0.01,r_val = 0, l_val = 0,temp2=0;

  
   // L.H.S of Equation
    arma::mat temp = arma::mat({0,0}).t();
    temp = x + alpha*p;    
    l_val = get_fun(temp);

   // R.H.S of Equation
    arma::mat prod;
    prod  = grad.t()*p; 
    temp2 = get_fun(x);
    r_val  = temp2 + alpha*gamma_1*prod(0,0); 

    while( l_val > r_val )
    {
      alpha = alpha * rho;
      temp = x + alpha*p;
      l_val = get_fun(temp);
      r_val =  temp2 + gamma_1*alpha*prod(0,0);
    } 
    return alpha;
}

void solver(arma::mat & x_0, std::vector<double> &x,double method)
{
    
    // count:count the number of iterations 
    double ratio  = 1,alpha = 0,count = 0,fun_val=0;
    arma::mat x_k_1 = arma::mat({0,0}).t();
    arma::mat p = arma::mat({0,0}).t();
    arma::mat grad = arma::mat({0,0}).t();
    //Only for BFGS
    arma::mat x_k = arma::mat({0,0}).t();
    arma::mat grad_k = arma::mat({0,0}).t();
    arma::mat B_k =  arma::mat({{1,0},{0,1}});
    arma::mat B_k_1 =  arma::mat({{0,0},{0,0}});
 
    x_k_1 = x_0;
    get_grad(grad,x_k_1);
    
    //Only for BFGS 
    x_k =  x_0;
    grad_k =  grad;
    

    while(ratio > 0.000001)
    {
        x_k_1 =  x_k + alpha*p; 
        //For gradient
        get_grad(grad,x_k_1);

        // To get P
        if(method ==  2)
            {
                get_p_BFGS(p,B_k_1,B_k,x_k_1,x_k,grad,grad_k,count);
            }
        else
            {
                get_p(p,x_k_1,grad,method);
            }
        //To get alpha
        alpha = get_alpha(x_k_1,grad,p);

        //To check that convergence criteria
        fun_val = get_fun(x_k_1);
        ratio  = sqrt(pow(grad(0,0),2) + pow(grad(1,0),2)) / (1 + fabs(fun_val)); 
        //Save Calculated value of X1,X2
        x.push_back(x_k_1(0,0));
        x.push_back(x_k_1(1,0));

        x_k = x_k_1;
            
        if(count > 2000)
        {
            break;
        }
        //For BFGS
        if(method == 2)
        {
            grad_k = grad;
            B_k = B_k_1;
        }
        count++;
    }
    std::cout<<"Number of Iteration : "<<count<<"\n";
    
}


int main()
{
    arma::mat x_0 = arma::mat({-3,-4}).t();
   //method represent the index for solver, 0:Stepest Descent, 1: Newton's , 2: BFGS
    double method = 0;
    double wtime=0,wtime1=0,wtime2=0;
    int temp;

// *********************************************Stepest Descent*********************************************
    std::cout<<"\nFor stepest Descent  "<<"\n\n";
    std::vector<double> x_stepest;
    wtime = omp_get_wtime();
    solver(x_0,x_stepest,method);
    wtime = omp_get_wtime()-wtime;
    std::cout<<"Time (in ms):"<<wtime*1000<<"\n";
    temp = x_stepest.size();
    std::cout<<"Value of X(x1,x2) : "<<x_stepest[temp-1]<< ","<<x_stepest[temp-2]<<"\n"<<"\n";

    char const * name = "Stepest.txt";
    write_file(x_stepest,name);


// *********************************************Newton's method*********************************************
    std::cout<<"For Newton's Method  "<<"\n\n";
    std::vector<double> x_newton;
    method = 1;
    wtime1 = omp_get_wtime();
    solver(x_0,x_newton,1);
    wtime1 = omp_get_wtime()-wtime1;
    std::cout<<"Time (in ms):"<<wtime1*1000<<"\n";
    temp = x_newton.size();
    std::cout<<"Value of X(x1,x2) : "<<x_newton[temp-1]<< ","<<x_newton[temp-2]<<"\n"<<"\n";

    char const * name1 = "newton.txt";
    write_file(x_newton,name1);

 // *********************************************BFGS method*********************************************
    std::cout<<"For BFGS Method  "<<"\n\n";
    std::vector<double> x_BFGS;
    method = 2;
    wtime2 = omp_get_wtime();
    solver(x_0,x_BFGS,2);   
    wtime2 = omp_get_wtime()-wtime2;
    std::cout<<"Time (in ms):"<<wtime2*1000<<"\n";
    temp = x_BFGS.size();
    std::cout<<"Value of X(x1,x2) : "<<x_BFGS[temp-1]<< ","<<x_BFGS[temp-2]<<"\n"<<"\n";
    
    char const * name2 = "BFGS.txt";
    write_file(x_BFGS,name2);

/*
   for(auto t=x_BFGS.begin(); t!=x_BFGS.end(); ++t)
   {
	 std::cout << *t <<"\n";
   }
*/

}


