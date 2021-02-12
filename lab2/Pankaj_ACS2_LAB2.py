#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


# # Problem 1

# In[3]:


test_point = np.random.uniform(0,1,size = (1,1000))
np.savetxt("Number_list_1.txt",test_point)


# In[4]:


temp = np.loadtxt("Inverse_random_sample.txt")


# In[14]:


Sample = np.linspace(0,max(temp),100)
plt.hist(temp,density = True, bins =  50)
plt.plot(Sample,np.exp(-Sample),label = r'$\rho(x)$')
plt.xlabel('X')
plt.ylabel('Distribution')
plt.legend()
plt.savefig("fig1")
plt.show()


# The result of Inverse sampling produce the exponential distribution as suspected. To validatate the results one can use QQ-plot or KS-test, which can vaildate the generate distribution.

# # Problem 2

# In[6]:


temp2 = np.loadtxt("Number_list_2.txt")


# In[7]:


stats.describe(temp2)


# In[8]:


#theoretical result
mean = np.sqrt(np.pi/2)
variance  = (4 - np.pi )/ 2
mean,variance


# In[9]:


temp3 = np.loadtxt("Number_list_3.txt")


# In[10]:


stats.describe(temp3)


# From C++ code:<br>
# 
# probability of acceptence,Uniform (Problem 2,part 1): 0.167661
# 
# Theoretical mean and variance : 1.253, 0.429204
# 
# Calculated mean and variance : 1.24775, 0.420909
# 
# probability of acceptence,Piecewise(Problem 2,part 2) : 0.312725
# 
# Theoretical mean and variance : 1.253, 0.429204
# 
# Calculated mean and variance : 1.25385, 0.429029

# # Problem 3

# In[16]:


temp4 = np.loadtxt("stats_sample_100.txt")
_ = plt.hist(temp4,density = True,bins = 25)
x = np.linspace(min(temp4),max(temp4), 100)
plt.plot(x, stats.norm.pdf(x, 0, np.sqrt(variance)),label = "Normal_dist")
plt.title('Stats For N=100')
plt.xlabel('X')
plt.ylabel('Distribution')
plt.legend()
plt.savefig("fig2")
plt.show()


# In[17]:


temp5 = np.loadtxt("stats_sample_1000.txt")
_ = plt.hist(temp5,density = True,bins = 25)
x = np.linspace(min(temp5),max(temp5), 100)
plt.plot(x, stats.norm.pdf(x, 0, np.sqrt(variance)),label = "Normal_dist")
plt.title('Stats for N=1000')
plt.xlabel('X')
plt.ylabel('Distribution')
plt.legend()
plt.savefig("fig3")
plt.show()


# As the number of observation increase from N=100 to N=1000, the generated distribution become more like Normal.This can also be acheived If we Increase the number of iterartion for a given N.

# In[ ]:




