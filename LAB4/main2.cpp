#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include<vector>
#include<math.h>

void read_file(std::vector<double> &x,std::vector<double> &y,std::vector<double> &z)
{
   std::ifstream file;
    file.open("data_Brandon.txt");
    while(!file.eof())
    {
        double a, b, c;
        file >> a >> b >> c; // extracts 3 floating point values seperated by whitespace
        x.push_back(a);
        y.push_back(b);
        z.push_back(c);
    }
}

double mean(std::vector<double> &x)
{

    double i=0;
    double m = 0;
    for(auto t=x.begin(); t!=x.end(); ++t)
        {
            m += *t;
            i++;
        }
    m = m/i;
    return m;
}


double variance(std::vector<double> &x)
{
    double MEAN = mean(x);
    double var = 0;
    int i=0;
   for(auto t=x.begin(); t!=x.end(); ++t)
    {
        var += pow( *t - MEAN,2);
        i++;
    }
   return var/(i-1);
}


double cov(std::vector<double> &x,std::vector<double> &y)
{
    double x_mean =0;
    double y_mean =0;
    x_mean = mean(x);
    y_mean = mean(y);
    std::vector<double> temp;
    std::vector<double> temp1;
    int index=0;
    for(auto t=x.begin(); t!=x.end(); ++t)
    {
        temp.push_back(*t - x_mean);
        index++;
    }
    for(auto t=y.begin(); t!=y.end(); ++t)
    {
        temp1.push_back(*t - y_mean);
    }
    
    double sum =0;
    for(int i=0;i<index;i++)
    {
        sum += temp[i]*temp1[i];
    }
    
    sum = sum/index;
    return sum;
}


double cor(std::vector<double> &x,std::vector<double> &y)
{
  double COV = 0; 
  COV =  cov(x,y);
  double COR  = COV/(sqrt(variance(x))*sqrt(variance(y)));
  return COR;
}

using namespace std;
int main()
{     
    std::vector<double> rain;
    std::vector<double> solar;
    std::vector<double> temprature;
    read_file(rain,solar,temprature);

    double COR_S_T = 0;
    double COR_S_R = 0;
    COR_S_T = cor(solar,temprature);
    COR_S_R = cor(solar,rain);
    std::cout<<"Correlation between Solar Radiation and rain is :"<<COR_S_R<<", and Between Solar radiation and Temprature is: "<<COR_S_T<<"\n";
}
