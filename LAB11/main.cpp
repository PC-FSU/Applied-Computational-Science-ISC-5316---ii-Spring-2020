#include<iostream>
#include<armadillo>
#include<math.h>


double w = 2/3;

double f(double x)
{
    double temp=0;
    temp = 4*exp(2*x);
    return temp;
}


double f1(double x)
{
    double temp=0;
    temp = exp(2*x);
    return temp;
}

void get_exact(arma::mat & exact, double (*f)(double x))
{
   int n = 0;
   n = exact.n_cols;
   double h=0;
   h = 1/(double) n;

   for(int i=0;i<n;i++)
   {
       exact(0,i) = (*f)(i*h);
   } 
}


void add_noise(arma::mat & u_0)
{
    int n = 0;
    n = u_0.n_cols;
    double noise; 
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-10,10);
    for( int i = 1; i<n; i++)
    {    
        u_0(0,i) = distribution(generator);
    }
}


void weighted_jacobi_step(arma::mat & u_0, arma::mat & u_1, double (*f)(double x))
{
    int n = 0;
    double h=0;
    n = u_0.n_cols;
    h = 1/(double) n;
    double temp = 0;


    for(int i=1;i<n-1;i++)
    {
            temp = 0.66666*(( u_0(0,i-1) + u_0(0,i+1)   -  (*f)(h*i)*pow(h,2) )/2  )  + (0.33333)*u_0(0,i);        
            u_1(0,i) = temp;
    }

}

void function_restriction(arma::mat & u_h,arma::mat & u_2h )
{
    int n = 0;
    n = u_h.n_cols/2;
   
    for(int i = 0;i<=n;i++)
    {
        if(i == 0 || i == n)
        {
            u_2h(0,i)  = u_h(0,2*i);
        }
        else
        {
            u_2h(0,i)  = 0.25*(u_h(0,2*i-1) + 2*u_h(0,2*i) + u_h(0,2*i+1));
        }
    }
}

void function_prolongation(arma:: mat & a_2h,arma::mat & a_h)
{
    int n = 0;
    n = a_2h.n_cols - 1;

    for(int i = 0;i<2*n+1;i++)
    {
        if(i%2 == 0 || i == 2*n)
        {
            a_h(0,i)  = a_2h(0,i/2);
        }
        else
        {
            a_h(0,i)  = 0.5*(a_2h(0,(i-1)/2) + a_2h(0,(i+1)/2));
        }
    }
}

void thomas_solver( arma:: mat & x,arma:: mat & a,arma:: mat & b,arma:: mat & c,arma:: mat & d )
{
    int n=0;
    n = x.n_cols;
    arma::mat d_prime(1,n,arma::fill::zeros);
    arma::mat c_prime(1,n-1,arma::fill::zeros);

    for(int i=0;i<n;i++)
    {
        if(i==0)
        {
            d_prime(0,i) = d(0,i)/b(0,i);
            c_prime(0,i) = c(0,i)/b(0,i);
        }
        else
        {
            if(i < n-1)
            {
              c_prime(0,i) = c(0,i)/( b(0,i) - a(0,i)*c_prime(0,i-1) );
            }
            d_prime(0,i) = (d(0,i) - a(0,i)*d_prime(0,i-1))/(b(0,1) - a(0,1)*c_prime(0,i-1));
        }

    }

    for(int i = n-1;i>=0;i--)
    {
        
        if(i == n-1)
        {
            x(0,i) = d_prime(0,i);
        }
        else
        {
            x(0,i) = d_prime(0,i) - c_prime(0,i)*x(0,i+1);
        }
    }

}


void get_matrix(arma::mat & A)
{
    int n=0;
    n = A.n_cols-1;
    double h=0;
    h = 1/(double) n;

    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            if(i==j)
            {
                A(i,j) = -2;
            }
            if( i-j == 1 || j-i == 1)
            {
                A(i,j) = 1;
            }
        }
    }
    A = A/pow(h,2);
}



int jacobi_solver(arma::mat & u,int n, int method, double tol)
{
  
  double h = 1/(double) n;
  arma::mat u_0(1,n+1,arma::fill::zeros); 
  arma::mat u_1(1,n+1,arma::fill::zeros);
  arma::mat exact(1,n+1,arma::fill::zeros);
  get_exact(exact,f1);
   
  // Add Noise
  add_noise(u_0);
  u_0(0,0) = 1; u_1(0,0) = 1;
  u_0(0,n) = exp(2); u_1(0,n) = exp(2);
  u = u_0;
  
  if (method == 1)
  { 
    int jacobi_iteration =  (int) tol;
    arma::mat Norm(1,jacobi_iteration+1,arma::fill::zeros);
    Norm(0,0) = arma::norm(u_0-exact,"fro");

    for(int i = 0;i<jacobi_iteration;i++)
    {
        weighted_jacobi_step(u_0,u_1,f);
        Norm(0,i+1) = arma::norm(u_1-exact,"fro");
        u_0 = u_1;
    }
    Norm.save("Norm.txt", arma::raw_ascii );
    return 0;
  }
  

  else if(method == 2)
  {
      double stop = 1;
      int counter = 0;
      std::cout<<" Hey in "<<method<<" with tol : "<<tol<<"\n";
      while(stop > tol)
      {
       //u_0.print("U_0");   
       weighted_jacobi_step(u_0,u_1,f);
       u_0 = u_1;
       //u_0.print("U_0");  
       stop =  arma::norm(u_1-exact,"fro");
       //stop =  stop/1.1;
       counter++;
       std::cout<<stop<<"\n";
      }
    return counter;
  }
  
}



int multigrid_solver( arma::mat & u,int n, int method, double tol)
{


  arma::mat u_0(1,n+1,arma::fill::zeros);    u_0 = u;
  arma::mat u_1(1,n+1,arma::fill::zeros);    u_1(0,0) = u_0(0,0); u_1(0,n) = u_0(0,n);
  arma::mat r_h(1,n+1,arma::fill::zeros);
  arma::mat r_2h(1,n/2+1);
  arma::mat eh(1,n+1,arma::fill::zeros);  
  arma::mat exact(1,n+1,arma::fill::zeros);
  get_exact(exact,f1);


  //for thomas solver 
  double h2 = 0;
  h2 = 1/(double) (0.5*n);
  arma::mat a(1,n/2+1,arma::fill::zeros); a.fill( 1);  a(0,0) = 0;          a = a/pow(h2,2);
  arma::mat b(1,n/2+1,arma::fill::zeros); b.fill(-2);                         b = b/pow(h2,2);
  arma::mat c(1,n/2+1,arma::fill::zeros); c.fill( 1);  c(0,n/2) = 0;        c = c/pow(h2,2);
  arma::mat e2h(1,n/2+1,arma::fill::zeros);
  
  arma::mat A(n+1,n+1,arma::fill::zeros);
  get_matrix(A);
 


    if(method == 1)
    {
    int multigrid_cycles = (int) tol;   
    arma::mat Norm(1,multigrid_cycles+1,arma::fill::zeros);    
    Norm(0,0) = arma::norm(u_0-exact,"fro"); 

        for(int i=0;i<multigrid_cycles;i++)
        {
            
            // 1 step
            for(int j =0;j<2;j++)
            {
                weighted_jacobi_step(u_0,u_1,f);      
                u_0 = u_1;
            }
            
            //2nd step
            r_h = exact -  (A*u_1.t()).t(); 
            r_h(0,0) = 0; r_h(0,n) = 0;
            //r_h.print("r_h");

            //3rd step
            function_restriction(r_h,r_2h);
            //r_2h.print("r_2h");
            

        
            //4th step
            thomas_solver(e2h,a,b,c,r_2h);
            e2h(0,0) = 0; e2h(0,n/2)=0;
            // e2h.print("e_2h");
            
            
            //5th step
            function_prolongation(e2h,eh);
            eh(0,0) = 0; eh(0,n)=0;
            //eh.print("e_h");  
            
            //6th step
            u_1 = u_1 + eh;

            //7th step
            u_0 = u_1;

            Norm(0,i+1) = arma::norm(u_1-exact,"fro");
        }
        //Norm_2.print("Norm:");
        Norm.save("Norm_2.txt", arma::raw_ascii );
        return 0;
    }

    else if(method == 2)
    {
        double stop = 1;
        int counter = 0;
        while(stop > tol)
        {
            // 1 step
            for(int j =0;j<2;j++)
            {
                weighted_jacobi_step(u_0,u_1,f);      
                u_0 = u_1;
            }
            
            //2nd step
            r_h = exact -  (A*u_1.t()).t(); 
            r_h(0,0) = 0; r_h(0,n) = 0;
            //r_h.print("r_h");

            //3rd step
            function_restriction(r_h,r_2h);
            //r_2h.print("r_2h");
            

        
            //4th step
            thomas_solver(e2h,a,b,c,r_2h);
            e2h(0,0) = 0; e2h(0,n/2)=0;
            // e2h.print("e_2h");
            
            
            //5th step
            function_prolongation(e2h,eh);
            eh(0,0) = 0; eh(0,n)=0;
            //eh.print("e_h");  
            
            //6th step
            u_1 = u_1 + eh;

            //7th step
            u_0 = u_1;

            stop = arma::norm(u_1-exact,"fro");
            counter++;
        }
        return counter;
    }
}




int main()
{
  int n = 0;
  n = 128;
  arma::mat u_0(1,n+1,arma::fill::zeros); 
  arma::mat u_1(1,n+1,arma::fill::zeros);
  arma::mat residual(1,n+1,arma::fill::zeros);
  arma::mat exact(1,n+1,arma::fill::zeros);

 
  //Get exact number
  get_exact(exact,f);
 
  // Add Noise
  add_noise(u_0);
  u_0(0,0) = 1; u_1(0,0) = 1;
  u_0(0,n) = exp(2); u_1(0,n) = exp(2);

  
  residual = arma::abs(exact - u_0);
  residual.save("0.txt",arma::raw_ascii);

  for(int i = 0;i<2;i++)
  {
      weighted_jacobi_step(u_0,u_1,f);
      u_0 = u_1;
  }
  

  residual = arma::abs(exact - u_1);
  residual.save("1.txt",arma::raw_ascii);
  
  //exact.print("exact:");

  //************** 2nd part ***************************
  arma::mat u_h(1,n+1,arma::fill::zeros);
  arma::mat u_2h(1,n/2+1,arma::fill::zeros);
  add_noise(u_h);
  function_restriction(u_h,u_2h);
  u_h.save("h.txt",  arma::raw_ascii);
  u_2h.save("2h.txt",arma::raw_ascii);
  
  
  //************** 3rd part ***************************
  n =  32;
  arma::mat a_2h(1,n+1, arma::fill::zeros);
  arma::mat a_h(1,2*n+1,arma::fill::zeros);
  add_noise(a_2h);
  function_prolongation(a_2h,a_h);
  a_2h.save("a_2h.txt",arma::raw_ascii);
  a_h.save("a_h.txt",  arma::raw_ascii);


  //************* 4TH PART ****************************
  n = 10; 
  arma::mat a(1,n,arma::fill::zeros); a.fill(-1); a(0,0) = 0;
  arma::mat b(1,n,arma::fill::zeros); b.fill(3);
  arma::mat c(1,n,arma::fill::zeros); c.fill(-1); c(0,n-1) = 0;
  arma::mat d(1,n,arma::fill::zeros); d.fill(1);  d(0,0) =  2; d(0,n-1) = 2;
  arma::mat x(1,n,arma::fill::zeros);  

  thomas_solver(x,a,b,c,d);
  
  x.print(" Value of X by Thomas Solver, x:");
  

  //************** 5th Part *****************
  int method = 0,     res = 0;    // if method == 1, solve by number of iteration, if method == 2, solve by error tolerance 
  double tol = 0;

  //*************** Part a ******************

  n = 256 ; method = 1 ;tol = 2000; arma::mat u_0_1;   // To get the intial vector to use in part B("Use the exact same initial vector as above")
  res = jacobi_solver(u_0_1,n,method,tol);
  

  //*************** Part B ******************

  method = 1;tol = 1000;
  res = multigrid_solver(u_0_1,n,method,tol);
  
   //************** Part C ******************
   /*
  n = 64; method = 2;     
  double tolerance_list[5] = {0.1,0.01,0.001,0.0001,0.00001};
  arma::mat    jacobi_iterate(0,5,arma::fill::zeros);
  arma::mat multigrid_iterate(0,5,arma::fill::zeros);
  
  arma::mat ref_mat(n+1,n+1,arma::fill::zeros); 
  add_noise(ref_mat);
  //ref_mat(0,0) = 1; ref_mat(n,n) = exp(2);
  
  for(int i=0;i<1;i++)
  {
      arma::mat ref_mat; 
      tol = tolerance_list[i];
      //res   =    jacobi_solver(ref_mat,n,method,tol);   
      //ref_mat.print("ref_1");

      res = multigrid_solver(ref_mat,n,method,tol);
      //ref_mat.print("ref_2");
  }    
*/
}

