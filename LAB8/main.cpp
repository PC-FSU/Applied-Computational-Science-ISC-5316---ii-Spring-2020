#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<algorithm>
#include<fstream>
#include <iterator>
#include <bits/stdc++.h>


int iteration = 10000;
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


double Fun(double x, double y)
{
    double temp = 0;
    temp = 1 + pow(sin(x),2) + pow(sin(y),2) - 0.1*exp(-pow(x,2)-pow(y,2));
    return temp; 
}

void get_Cloud(std::vector< std::vector<double> > &x,double N)
{
 
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-20, 20);
    double temp2 = 0;
    for (int n = 0; n < N; ++n) 
    {
        // Use dis to transform the random unsigned int generated by gen into a 
        // double in [20, 20). Each call to dis(gen) generates a new random double
        std::vector<double> temp;
        temp.push_back(dis(gen));
        temp.push_back(dis(gen));
        temp2 = Fun(temp[0],temp[1]);
        temp.push_back(temp2);
        x.push_back(temp);
    }
}

void get_subspace(int * index_list,int n,int N)
{
  static std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution (0,N-1);
  for(int i=0;i<=n;i++)
  {
    index_list[i] = std::round(distribution(generator));
  }
}

double get_max(std::vector< std::vector<double> > &x,double N)
{

    std::vector<double> temp;
    for(int i = 0;i<N;i++)
    {
        temp.push_back(x.at(i).at(2));
    }    
    double maxElementIndex = std::max_element(temp.begin(),temp.end()) - temp.begin();
    double maxElement = *std::max_element(temp.begin(), temp.end());
    return maxElementIndex;
}

int Find_new(std::vector< std::vector<double> > &x,int * index_list,int n,int N)
{   
    double val_1 = 0,val_2=0,val_3=0,MAX,x_centroid = 0,y_centroid = 0,P_x,P_y,P_val;
    int min_index = 0,max_index = 0,status=0;
    val_1 = x.at(index_list[0]).at(2);
    val_2 = x.at(index_list[1]).at(2);
    val_3 = x.at(index_list[2]).at(2);
    MAX = std::max({val_1,val_2,val_3});
    
    //To find the index of max element
    for(int i = 0;i<=n;i++)
    {
        if(index_list[i]==MAX)
        {
            max_index = i;
        }
    }
   
   for (int i =0;i<=n;i++)
   { 
          if(i != max_index)
          {
            x_centroid += x.at(index_list[i]).at(0);
            y_centroid += x.at(index_list[i]).at(1);
          }
   }

   x_centroid = x_centroid/(n); y_centroid = y_centroid/(n);
   P_x = 2*x_centroid -  x.at(index_list[max_index]).at(0);
   P_y = 2*y_centroid -  x.at(index_list[max_index]).at(1);
   P_val = Fun(P_x,P_y);

   // Reusing val_1 to find the index of global minima 
   val_1 = get_max(x,N);
   //Resuing val_2 to find the value of global minima
   val_2 = x.at(val_1).at(2);
   if( ( P_val < val_2 ) && ( ( (P_x <= 20) && (P_x >= -20) ) && ( (P_y <= 20) && (P_y >= -20) ) ) )
   {
       x.at(val_1).at(0) = P_x;
       x.at(val_1).at(1) = P_y;
       x.at(val_1).at(2) = P_val;
       status = 1;
   }
   return status;
}

int main()
{
   
    int N = 0, n=0, flag=0,total_iterations=0;
    int index_list[n+1] = {};
    N=1000;
    n=2;
    std::vector< std::vector<double> > coords;
    get_Cloud(coords,N);
    
    while(flag < iteration)
    {   
        // To check that we found point p.
        int status = 0;   
        get_subspace(index_list,n,N);       
        status = Find_new(coords,index_list,n,N);
        
        // Updata the counter if we are able to find the point p.
        if(status == 1)
        {
            flag++;
        }

        if(flag==1000)
        {   
            const char * name = "Coords_1000";
            write_file(coords,name);
        }
        if(flag==5000)
        {   
            const char * name1 = "Coords_5000";
            write_file(coords,name1);
        }
        if(flag==7500)
        {   
            const char * name2 = "Coords_7500";
            write_file(coords,name2);
        }
        if(flag==10000)
        {   
            const char * name3 = "Coords_10000";
            write_file(coords,name3);
        }
        total_iterations++;
    } 
    std::cout<<total_iterations<<"\n";
}
