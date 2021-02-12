#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>   


void read_file(std::vector<double> &x,std::vector<double> &y ,const char S[])
{
    std::ifstream file;
    file.open(S, std::ifstream::in | std::ifstream::binary);
    int i=0;
    while(!file.eof())
    {
        double a, b;
        file >> a >> b; // extracts 2 floating point values seperated by whitespace
        //std::cout<<a<<" "<<b<<"\n";
         x.push_back(a);
         y.push_back(b);
    }
}

void get_prediction(double * beta, double * Model_Predicted_y,std::vector<double> &x, int problem)
{
int i=0;
if(problem == 1)
{
    for(auto t=x.begin(); t!=x.end(); ++t)
    {
        Model_Predicted_y[i] = beta[0] + (*t)*beta[1];
        i++;
    }
}
else
    {
        for(auto t=x.begin(); t!=x.end(); ++t)
        {
            Model_Predicted_y[i] = beta[0] + (*t)*beta[1] + pow(*t,2)*beta[2];
            i++;
        }
    }
}


double get_SS_res(double * Model_Predicted_y,std::vector<double> &Y)
{
    int i=0;
    double res=0;
    for(auto t=Y.begin(); t!=Y.end(); ++t)
        {
            res += pow(*t - Model_Predicted_y[i],2);
            i++;
        }
    return res;
}

double get_SS_tot(std::vector<double> &Y)
{
    double i=0;
    double mean = 0;
    for(auto t=Y.begin(); t!=Y.end(); ++t)
        {
            mean += *t;
            i++;
        }
    mean = mean/i;

    double res=0;
    for(auto t=Y.begin(); t!=Y.end(); ++t)
        {
            res += pow(*t - mean,2);
        }
    return res;

}


int main()
{

std::vector<double> x1;
std::vector<double> y1;
std::vector<double> x2;
std::vector<double> y2;

char const * name = "sample1.dat";
read_file(x1,y1,name);
char const * name1 = "sample2.dat";
read_file(x2,y2,name1);

int File_Size = x1.size();
int File_Size2 = x2.size();
int nrow = File_Size-1; 
int col = 0;

// write to file
std::ofstream myfile1;
myfile1.open("Problem1_result.txt");
std::ofstream myfile2;
myfile2.open("Problem2_result.txt");

// ********************************************Problme 2, Part 1*************************************************


std::vector<double> tempX;
std::vector<double> tempY;
// When Problem ==1, solve for matrix with 2 columns, and for 3 when problem ==2
for(int problem = 1;problem<3;problem++)
{
    if(problem == 1)
    {
        col  = 2;
    }
    else
    {
        col = 3;
    }

double X[nrow][col]={};
double Xt[col][nrow]={};

// when iteration == 0, use data from sample1.dat and from sample2.dat if it is 1.
for(int iteration = 0 ; iteration < 2;iteration++)
{
   if(iteration == 0)
   {
       tempX = x1;   // intialize x, and y to some temp vector.
       tempY = y1;
   }
   else
   {
       tempX = x2;
       tempY = y2;
   }
     
   
    //Set MATRIX X
    int index = 0;
    for(auto t=tempX.begin(); t!=tempX.end()-1; ++t)
        {  
            X[index][0] = 1;
            X[index][1] = *t;
            if(problem == 2)
            {
                X[index][2] = pow(*t,2);
            }
            index ++;
        }
   
    // Set MATRIX X.T 
    index=0;
    for(auto t=tempX.begin(); t!=tempX.end()-1; ++t)
        {  
            Xt[0][index] = 1;
            Xt[1][index] = *t;
            if(problem == 2)
            {
                Xt[2][index] = pow(*t,2);
            }
            index ++;
        }


    //Calculate X*X.T
    double XXt[col][col] = {}; 
        for(int i = 0;i<col;i++)
        {
            
            for(int j=0;j<col;j++)
            {
            double temp=0;
            for(int k=0;k<nrow;k++)
                {
                     temp += Xt[i][k]*X[k][j];
                }
                XXt[i][j] = temp;
            }
        }

  

    //Inverse of XX.t
    double Inverse[col][col] = {};
    double deter =0;
    if(problem == 1)
    {
        deter = XXt[0][0]*XXt[1][1] - XXt[1][0]*XXt[0][1];
        Inverse[0][0] = XXt[1][1]/deter;
        Inverse[1][1] = XXt[0][0]/deter;
        Inverse[1][0] = -1*XXt[1][0]/deter;
        Inverse[0][1] = -1*XXt[0][1]/deter;
        std::cout<<Inverse[0][0]<<" "<<Inverse[0][1]<<" "<<Inverse[1][0]<<" "<<Inverse[1][1]<<"\n";
           
    }
    else
    {
    
    //finding determinant for 3*3
	for(int i = 0; i < col; i++)
		{
            deter = deter + (XXt[0][i] * (XXt[1][(i+1)%3] * XXt[2][(i+2)%3] - XXt[1][(i+2)%3] * XXt[2][(i+1)%3]));
        }
	
  
    // Finding Inverse for 3*3
    for(int i = 0; i < col; i++)
    {
		for(int j = 0; j < col; j++)
			{
              Inverse[i][j] = ((XXt[(j+1)%3][(i+1)%3] * XXt[(j+2)%3][(i+2)%3]) - (XXt[(j+1)%3][(i+2)%3] * XXt[(j+2)%3][(i+1)%3]))/ deter;
              std::cout<<Inverse[i][j]<<" ";
            }
        std::cout<<"\n";    
	}

    }

    //Calculate XtY
    double XtY[col] = {};
    for(int i=0;i<col;i++)
    {  double temp =0;
        for(int j=0;j<nrow;j++)
        {
            temp += Xt[i][j]*tempY[j];
        }
        XtY[i] = temp;
    }
 
    //Calculate Beta
    double beta[col] = {};
    if(problem == 1)
    {
        beta[0] = Inverse[0][0]*XtY[0] + Inverse[0][1]*XtY[1];
        beta[1] = Inverse[1][0]*XtY[0] + Inverse[1][1]*XtY[1];
    }
    else
    {
        beta[0] = Inverse[0][0]*XtY[0] + Inverse[0][1]*XtY[1] + Inverse[0][2]*XtY[2];
        beta[1] = Inverse[1][0]*XtY[0] + Inverse[1][1]*XtY[1] + Inverse[1][2]*XtY[2];
        beta[2] = Inverse[2][0]*XtY[0] + Inverse[2][1]*XtY[1] + Inverse[2][2]*XtY[2];
    }


    // For Parametres
    double Model_Predicted_y[nrow] = {};
    get_prediction(beta,Model_Predicted_y,tempX,problem);
    double SS_res = 0;
    SS_res = get_SS_res(Model_Predicted_y,tempY);
    double SS_tot = 0;
    SS_tot = get_SS_tot(tempY);
    double R2 =0;
    R2 = 1 - SS_res/SS_tot;

    // Write to File and for console window
    if(problem == 1)
    {   
        myfile1 << beta[0]<<" "<<beta[1]<<" "<<R2<<"\n";
        std::cout<<"For Problem 1.1 :"<<"\n";
        std::cout<<"\n"<<"For File "<<iteration<<", Value of beta_0 : "<<beta[0]<<" and Value of beta_1 :"<<beta[1]<<" Value of R2 is = "<<R2<<"\n";
    }
    else 
    { 
        myfile2 << beta[0]<<" "<<beta[1]<<" "<<beta[2]<<" "<<R2<<"\n";
        std::cout<<"For Problem 1.2 "<<"\n";
        std::cout<<"\n"<<"For File "<<iteration<<" Value of beta_0 : "<<beta[0]<<", Value of beta_1 :"<<beta[1]<<", Value of beta_2 :"<<beta[2]<<" Value of R2 is = "<<R2<<"\n";
    }
  
 
}

}

// Close the file.
myfile1.close();
myfile2.close();
}



/*
for(auto t=x1.begin(); t!=x1.end(); ++t)
	std::cout << *t << '\n';
for(auto t=y1.begin(); t!=y1.end(); ++t)
	std::cout << *t << '\n';
*/



/*
for(int i=0;i<nrow;i++)
{
    for(int j=0;j<col;j++)
    {
      std::cout<<"X["<<i<<"]["<<j<<"] = "<<X1[i][j]<<"  ";
    } 
    std::cout<<"\n";
}

for(int i=0;i<col;i++)
{
    for(int j=0;j<nrow;j++)
    {
      std::cout<<X1t[i][j]<<" ";
    } 
    std::cout<<"\n";
}

*/


// for Problem 2
/*
void set_X(std::vector<double> &x,double X[nrow][col])
{
    
        int index = 0;
        for(auto t = x.begin(); t!= x.end()-1; ++t)
        {  
            X[index][0] = 1;
            X[index][1] = *t;
            index ++;
        }
     
}

void set_Xt(std::vector<double> &x,double Xt[col][nrow])
{
        int index1=0;
        for(auto t=x.begin(); t!=x.end()-1; ++t)
        {  
            Xt[0][index1] = 1;
            Xt[1][index1] = *t;
            index1 ++;
        }
}

*/