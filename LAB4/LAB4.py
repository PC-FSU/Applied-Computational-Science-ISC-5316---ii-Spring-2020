#!/usr/bin/env python
# coding: utf-8

# In[219]:


import numpy as np
import matplotlib.pyplot as plt


# In[ ]:





# In[198]:


#Read from file
temp = np.loadtxt('sample1.dat')
X1 = temp[:,0].reshape(-1,1)
Y1 = temp[:,1].reshape(-1,1)

temp2 = np.loadtxt('sample2.dat')
X2 = temp2[:,0].reshape(-1,1)
Y2 = temp2[:,1].reshape(-1,1)


# In[203]:


#Read beta from file
beta_1= np.loadtxt('Problem1_result.txt')
beta_2= np.loadtxt('Problem2_result.txt')


# In[223]:


def plot(x,y,beta,problem,name):
    if(problem == 1):
        plt.plot(x,y,label ='True')
        y_pred = beta[0] + beta[1]*x;
        plt.plot(x,y_pred,label ='Fitted')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.savefig(name)
    if(problem == 2):
        plt.plot(x,y,label ='True')
        y_pred = beta[0] + beta[1]*x + beta[2]*np.square(x);
        plt.plot(x,y_pred,label ='Fitted')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(name)


# #### For Problem 1,Part 1

# In[224]:


beta_temp = beta_1[0][:2]
plot(X1,Y1,beta_temp,1,"File 1")


# In[225]:


beta_temp = beta_1[1][:2]
plot(X2,Y2,beta_temp,1,"File 2")


# #### Problem 1,Part 2

# In[229]:


beta_temp = beta_2[0][:3]
plot(X1,Y1,beta_temp,2,"File 2,1")


# In[230]:


beta_temp = beta_2[1][:3]
plot(X2,Y2,beta_temp,2,"File 2,2")


# #### Problem 2
# To preprocess the file having informaton about temprature,rain and solar radiation.

# In[231]:


with open('data.txt') as f:
    lines = list((line for line in f if not line.startswith('T')))
    
    #FH = np.loadtxt(lines, delimiter=deli)
x=[]
y=[]
z=[]
for line in lines:
     first = line.split(' ')[1]
     x.append(float(first.split(',')[1]))
     y.append(float(first.split(',')[2]))
     z.append(float(first.split(',')[3].split('\n')[0]))   
x = np.array(x).reshape(-1,1)
y = np.array(y).reshape(-1,1)
z = np.array(z).reshape(-1,1)

res = np.concatenate((y,z), axis=1)
res = np.concatenate((x,res),axis =1)
np.savetxt('data1.txt',res, fmt='%10.5f')

