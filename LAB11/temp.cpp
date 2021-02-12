












  /*
        double temp = 0;
        for(int i = 1;i<n;i++)
        {  
                 temp = f1(i*h);
                 r_h(0,i) = temp - ( b_1(0,i)*u_1_2(0,i) + a_1(0,i)*u_1_2(0,i-1) + c_1(0,i)*u_1_2(0,i+1) );    
        }
        r_h.print("r_h"); 
        
    







   
  
  n = 256;
  double h = 1/(double) n;
  int jacobi_iteration = 2000;
  arma::mat u_0_1(1,n+1,arma::fill::zeros); 
  arma::mat u_1_1(1,n+1,arma::fill::zeros);
  arma::mat Norm(1,jacobi_iteration+1,arma::fill::zeros);
  arma::mat exact_1(1,n+1,arma::fill::zeros);
  arma::mat temp_vec;
  get_exact(exact_1,f1);
  
  // Add Noise
  add_noise(u_0_1);
  u_0_1(0,0) = 1; u_1_1(0,0) = 1;
  u_0_1(0,n) = exp(2); u_1_1(0,n) = exp(2);
  temp_vec = u_0_1; 
  Norm(0,0) = arma::norm(u_0_1-exact_1,"fro");

  for(int i = 0;i<jacobi_iteration;i++)
  {
      weighted_jacobi_step(u_0_1,u_1_1,f);
      Norm(0,i+1) = arma::norm(u_1_1-exact_1,"fro");
      u_0_1 = u_1_1;
  }
  
  Norm.save("Norm.txt", arma::raw_ascii );



  int multigrid_cycles = 1000;

  arma::mat u_0_2(1,n+1,arma::fill::zeros);                               u_0_2 = temp_vec;
  arma::mat u_1_2(1,n+1,arma::fill::zeros);                               u_1_2(0,0) = u_0_2(0,0); u_1_2(0,n) = u_0_2(0,n);
  arma::mat Norm_2(1,multigrid_cycles+1,arma::fill::zeros);               Norm_2(0,0) = arma::norm(u_0_2-exact_1,"fro");
  arma::mat r_h(1,n+1,arma::fill::zeros);
  arma::mat r_2h(1,n/2+1);
  arma::mat eh(1,n+1,arma::fill::zeros);  

  //for thomas solver 
  double h2 = 0;
  h2 = 1/(double) (0.5*n);
  arma::mat a_1(1,n/2+1,arma::fill::zeros); a_1.fill( 1);  a_1(0,0) = 0;          a_1 = a_1/pow(h2,2);
  arma::mat b_1(1,n/2+1,arma::fill::zeros); b_1.fill(-2);                         b_1 = b_1/pow(h2,2);
  arma::mat c_1(1,n/2+1,arma::fill::zeros); c_1.fill( 1);  c_1(0,n/2) = 0;        c_1 = c_1/pow(h2,2);
  arma::mat e2h(1,n/2+1,arma::fill::zeros);
  
  arma::mat A(n+1,n+1,arma::fill::zeros);
  arma::mat A2h(n/2+1,n/2+1,arma::fill::zeros);
  get_matrix(A2h);  
  get_matrix(A);
 
  
    for(int i=0;i<multigrid_cycles;i++)
    {
        
        // 1 step
        for(int j =0;j<2;j++)
        {
            weighted_jacobi_step(u_0_2,u_1_2,f);      
            u_0_2 = u_1_2;
        }
        
        //2nd step
      
        r_h = exact_1 -  (A*u_1_2.t()).t(); 
        r_h(0,0) = 0; r_h(0,n) = 0;
        //r_h.print("r_h");

        //3rd step
        function_restriction(r_h,r_2h);
        //r_2h.print("r_2h");
        

       
        //4th step
        thomas_solver(e2h,a_1,b_1,c_1,r_2h);
        e2h(0,0) = 0; e2h(0,n/2)=0;
       // e2h.print("e_2h");
       
       
        e2h = (arma::solve(A2h,r_2h.t())).t();
        e2h(0,0) = 0; e2h(0,n/2)=0;
        e2h.print("e_2h");
       
        
        //5th step
        function_prolongation(e2h,eh);
        eh(0,0) = 0; eh(0,n)=0;
        //eh.print("e_h");  
        
        //6th step
        u_1_2 = u_1_2 + eh;

        //7th step
        u_0_2 = u_1_2;

        Norm_2(0,i+1) = arma::norm(u_1_2-exact_1,"fro");
    }
    //Norm_2.print("Norm:");
    Norm_2.save("Norm_2.txt", arma::raw_ascii );
