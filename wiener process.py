# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 20:11:17 2018

@author: Matthew Hamm
"""
import timeit
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import norm,linregress
from scipy.spatial.distance import cdist
n = 1000
l=36
def createwalk(N):
    ''' Creates a self avoiding random walk using the pivot algorithm and gives a 2d numpy array
    Parameters
    -----------
    N: integer, which denotes the length of the random walk
    
    Returns
    --------
    chain: A 2d numpy array, with the coords of th random walk
    '''
    chain=np.zeros(N)
    
    for i in range(1, N):
        val = random.randint(1, 2)
        step=random.gauss(0, 1)
        if val == 1:
            chain[i] = chain[i - 1] + step

        elif val == 2:
            chain[i] = chain[i - 1] - step

    
        
    return chain

x=(createwalk(n))


def lengthfind(L,X):
    '''Finds the length between points a certain number of step apart
    Parameters
    -----------
    L : an integer for the step gap
    X:A 2d numpy array for a walk
    Returns
    --------
    Length:A numpy array of lengths a certain number of steps apart'''
    length=[]

    for i in range(L,(np.size(X))):
        length.append(X[i]-X[i-L])

    
    length=np.array(length)

    return length
length=lengthfind(l,x)
np.savetxt('length.txt',length)
# save numpy array to a text doucment

def plotdist(Length):
    '''Plots a distribution as a histogram and finds a relevant gaussian
    Parameters
    ------------
    Length: A numpy array of the distribution'''
    nr_in_bin, bin_edges=np.histogram(Length)
    width = bin_edges[1:] - bin_edges[:-1]
    center = (bin_edges[:-1] + bin_edges[1:]) / 2   
                           
    s = 1/len(Length)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(center, nr_in_bin*s, align='center', width=width*0.9, 
           color='g', label='Histogram of lengths')
    mu, std = norm.fit(length)
    p=norm.pdf(center, mu, std)
    ax.plot(center,p, label='Gaussian of lengths')
    ax.set_xlabel('Lengths for step gap '+ str(l))
    ax.set_ylabel('Probablity density')
    ax.legend()
    ax.set_title('Distribution of lengths of a step gap of '+str(l) +' for a random walk')
    #plotting stuff:
    plt.savefig('Distribution of lengths of a step gap of '+str(l) +' for a random walk')
    plt.show()
    print(str(mu) + ' with error ' +str(std))

plotdist(length)
def plotvar(X):
    '''Plot the variance of the lengths between points number of steps apart
    Parameters
    -----------

    X :A 2d numpy array for a walk'''
    variance=[]
    nsteps=np.logspace(1/2,3/2,400,dtype=np.int64)
    
    for i in nsteps:

        length=lengthfind(i,X)



        variance.append(np.var(length))
    slope, intercept, r_value, p_value, std_err = linregress(np.log(nsteps), np.log(variance))
    fig=plt.figure()
    ax = fig.add_subplot(1, 1, 1) 
    ax.set_yscale('log')
    ax.set_ylabel('Variance of length')
    ax.set_xlabel('Step gap')
    ax.set_xscale('log')
    ax.set_title('Plot of the Variance Length against step gap')
    ax.plot(nsteps,variance, label='Results')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope), np.exp(intercept)), label='linear regression')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope+std_err), np.exp(intercept)),'.k',label='upper error')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope-std_err), np.exp(intercept)),'.k', label='lower error')
    ax.plot(nsteps,np.multiply(nsteps, np.exp(intercept)), label= 'General esitmate')
    ax.legend()
    plt.savefig('Plot of the Variance Length against step gap')
    plt.show()
    print(str(slope) + ' with error ' +str(std_err) +'with intercept' +str(intercept))

def findradius(X):
    '''Finds the minimum radius of a sphere to contain the random walk
    Parameters
    -----------
    X :A 2d numpy array for a walk
    Returns
    --------
    Radius: A float of the minium radius'''
    maxlength=[]
    N=np.size(X,0)

    nsteps=np.linspace(0,N-1,N,dtype=np.int64)

    for i in nsteps:

        length=lengthfind(i,X)
        maxlength.append(np.max(length))
    radius= np.max(np.array(maxlength)/2)
    return radius
plotvar(x)
def plotmassvol(n):
    '''Plots the mass against volume containing walks'''
    r=[]
    N=np.logspace(3/2,np.log10(n),400,dtype=np.int64)
    for i in N:
        x=(createwalk(i))
        r.append(findradius(x))
    
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Volume')
    ax.set_xlabel('Mass')
    ax.set_title('Plot of mass of the total 1d random walk against radius of to contain it')
    slope, intercept, r_value, p_value, std_err = linregress(np.log(N), np.log(r))
    ax.plot(N,r)
    ax.plot(N,np.multiply(np.power(N,slope),np.exp(intercept)), label='Regression')
    ax.plot(N,np.multiply(np.power(N,slope+std_err),np.exp(intercept)),'.k',label='upper error')
    ax.plot(N,np.multiply(np.power(N,1/2),np.exp(intercept)), label= 'General esitmate')
    ax.plot(N,np.multiply(np.power(N,slope-std_err),np.exp(intercept)), '.k', label='lower error')
    ax.legend()
    plt.savefig('Plot of mass of the total 1d random walk against radius of to contain it')
    plt.show()
    print(str(slope) + ' with error ' +str(std_err))
plotmassvol(n)
def spherecount(x):

    count=[]
    length=x
    N=findradius(x)
    r=np.logspace(1/2,(np.log10(N)),40)
    for i in r :
        count.append(np.size(length[length<i]))
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Volume')
    ax.set_xlabel('Mass')
    ax.set_title('Plot of mass against volume of a random walk')
    slope, intercept, r_value, p_value, std_err = linregress(np.log(count), np.log(r))
    ax.plot(count,r)
    ax.plot(count,np.multiply(np.power(count,slope),np.exp(intercept)), label='Regression')
    ax.plot(count,np.multiply(np.power(count,slope+std_err),np.exp(intercept)),'.k',label='upper error')
    ax.plot(count,np.multiply(np.power(count,1/2),np.exp(intercept)), label= 'General esitmate')
    ax.plot(count,np.multiply(np.power(count,slope-std_err),np.exp(intercept)), '.k', label='lower error')
    ax.legend()
    plt.savefig('Plot of mass against volume of a random walk')
    plt.show()
    print(str(slope) + ' with error ' +str(std_err))
    return slope
   


