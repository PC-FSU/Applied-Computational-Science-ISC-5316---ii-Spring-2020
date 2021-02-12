#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include <algorithm>    // std::random_shuffle
#include<math.h>
#include<random>


// To write to file
void write_file(double *list,int n,const char S[])
{
    std::ofstream file;
    file.open(S,std::ios::out|std::ios::trunc);
    if(file.is_open())
        {
        for(int i=0;i<n;i++)
            {
               file<<i<<" ";
               file<<list[i]<<" ";
            }
        file.close();
        }
    else
        {
            std::cout<<"unable to open the file"<<"\n";
        }
}


// to generate random index for randomly picking chunk of array
int get_index(int File_Size,int n)
{
    int temp = 0;
    bool temp1 = true;
    while(temp1)
    {
        temp = rand()%File_Size;
        if (abs(File_Size -  temp) > n)
        {
            temp1 = false;
        } 
    }
    return temp;
}

// To create subvector from generated random vector
std::vector<int> sub_vector(std::vector<int> const &data, int m ,int n)
{
    auto first = data.cbegin() + m;
    auto last = data.cbegin() + n+1;
    std::vector<int> temp(first,last);
    return temp;
}

// To get rejection probability
void get_proportion(int ** result, double * proportion_list,int row,int columns)
{
 for(int i = 0;i<row;i++)
 { 
     double count = 0;
     double total = 0;
     for(int j =0;j<columns;j++)
     {
        if(result[i][j] == 0)
        {
            count++;
        }
        total++;      
     }
     proportion_list[i] = count/total;
     
 }
}

// To find probability of rejection
void function(std::vector<int> &data, int * list_n, int Number_of_experiment, int **result,int File_Size)
{
    

   for(int j = 0 ; j<Number_of_experiment;j++)
   { 
       for(int i = 0;i<5;i++)
        {
         double count = 0;
         double p = 0;
         double SD = 0;
         int index = get_index(File_Size,list_n[i]);    
         std::vector<int> Sub_Vector = sub_vector(data,index,list_n[i]+index);
         for (auto k = Sub_Vector.begin(); k != Sub_Vector.end(); ++k) 
         {
             if(*k == 1)
             {
                 count++;
             }
         }

         p = count/list_n[i];
         SD = sqrt((p*(1-p))/list_n[i]);
         if((p + 2*SD) >= 0.5 && (p - 2*SD)<=0.5)
         {
               result[i][j] = 1;
         }
         //std::cout<<count;
         //std::cout<<i<< " ";
        }
        //std::cout<<j<<" ";
   }

}


int main()
{  
    srand(time(0));
    //To read from binary file
    std::vector<int> data;
    std::ifstream file;
    file.open("survey1.dat", std::ifstream::in | std::ifstream::binary);
    std::string a;
    int temp = 0;
    while(!file.eof())
    {
        std::getline(file,a,',');
        temp  = std::stoi(a);
        data.push_back(temp);
    }
    int File_Size = data.size();
    std::cout<<File_Size<<"\n";
    
    // initialize the 2D dynamic list which holds the results of 10000 experiment for different N.
    int list_n[] = {50,100,250,500,1000};
    int Number_of_experiment = 10000;
    int **result = new int*[5];
    for(int i =0;i<5;i++)
        {
            result[i] = new int[Number_of_experiment]; 
        }

    for(int i =0;i<5;i++)
    {
        for(int j = 0; j <Number_of_experiment;j++)
        {
            result[i][j] = 0; 
        }
    }

    // Run the experiments
    function(data,list_n,Number_of_experiment,result,File_Size);

    // calculate the results of experiments(proportion) for different n
    int row =  5;
    int columns = Number_of_experiment;
    double proportion_list[5]  = {};
    get_proportion(result,proportion_list,row,columns);
   
   
    //Write proportions to file
    char const * name = "Result.txt";
    write_file(proportion_list,5,name);


    //Free each sub-array
    for(int i = 0; i < 5; ++i)
    {
        delete[] result[i];   
    }
    //Free the array of pointers
    delete[] result;
    
}