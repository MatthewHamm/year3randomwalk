# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 20:11:17 2018

@author: Matthew Hamm
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import norm,linregress
from scipy.spatial.distance import cdist
import timeit
n = 200
l=25
def createwalk(N):
    ''' Creates a self avoiding random walk and gives a 2d numpy array
    Parameters
    -----------
    N: integer, which denotes the length of the random walk
    
    Returns
    --------
    chain_init: A 2d numpy array, with the coords of the random walk
    '''
    i=1

    chain=np.zeros([N,2])
    
    while i<N :
        val = random.randint(1, 4)
        if val == 1 and not (np.any(chain[chain[:,0]==chain[i-1,0]+1,1]==chain[i-1,1])):
            chain[i,0] = chain[i - 1,0] + 1
            chain[i,1] = chain[i - 1,1]
            i+=1
        elif val == 2 and not (np.any(chain[chain[:,0]==chain[i-1,0]-1,1]==chain[i-1,1])):
            chain[i,0] = chain[i - 1,0] - 1
            chain[i,1] = chain[i - 1,1]
            i+=1
        elif val == 3 and not (np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]+1)):
            chain[i,0] = chain[i - 1,0] 
            chain[i,1] = chain[i - 1,1] +1
            i+=1
        elif val==4 and not (np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]-1)):
            chain[i,0] = chain[i - 1,0]
            chain[i,1] = chain[i - 1,1]-1
            i+=1
        elif np.any(chain[chain[:,0]==chain[i-1,0]+1,1]==chain[i-1,1]) and np.any(chain[chain[:,0]==chain[i-1,0]-1,1]==chain[i-1,1]) and np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]+1) and np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]-1):
            chain=np.zeros([N,2])
            i=1
    
    return chain
def creategrid(N):
    '''Creates a 2d list of strings to print out later
    Paremters
    ----------
    n: integer denoting the size of the grid to create
    
    Returns
    --------
    grid: a 2d list of strings, which are spaces
    '''
    grid=[]


    grid=[[' ' for _ in range(2*N)] for i in range(2*N) ]

    return grid

def printgrid(X):
    '''Prints a grid displaying x 
    Parameters
    ----------
    X: A 2d numpy array for a path
'''

    N0=np.max(np.abs(X[:,0]))
    N1=np.max(np.abs(X[:,1]))
    #size of grid found from path
    N=abs(max(N0,N1))

    grid=creategrid(2*N+1)
    for i in range (np.size(X,0)):
    
        coordx=N+X[i,0]
        coordy=N+X[i,1]
        choice=x[i]-x[i-1]

    
        if choice[0]==1:
            grid[2*coordy][2*coordx]="."
            grid[2*coordy][2*coordx-1]="-"
        elif choice[0]==-1:
            grid[2*coordy][2*coordx]="."
            grid[2*coordy][2*coordx+1]="-"
        elif choice[1]==1:
            grid[2*coordy][2*coordx]="."
            grid[2*coordy-1][2*coordx]="|"
        elif choice[1]==-1:
            grid[2*coordy][2*coordx]="."
            grid[2*coordy+1][2*coordx]="|"
    
    
    grid[2*N][2*N]="."
    print(len(grid))
    for j in range(len(grid)):
        print('')
        for i in range(len(grid)):
            print(grid[len(grid)-j-1][i], end='')
    print('')
x=np.int64(createwalk(n))

printgrid(x)
def lengthfind(L,X):
    '''Finds the length between points a certain number of step apart
    Parameters
    -----------
    L : an integer for the step gap
    X:A 2d numpy array for a walk
    Returns
    --------
    Length:A numpy array of lengths a certain number of steps apart'''
    lengthy=np.array([])
    lengthx=np.array([])
    for i in range(L,(np.size(X,0))):
        lengthx=np.append(lengthx,X[i,0]-X[i-L,0])
        lengthy=np.append(lengthy,X[i,1]-X[i-L,1])
    length=np.sqrt(np.square(lengthy)+np.square(lengthx))
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
    fig=plt.figure()
    ax = fig.add_subplot(1, 1, 1)        
    ax.set_title("Random Walk ($n = " + str(n) + "$ steps)")
    ax.plot(x[:,0],x[:,1],'r' )
    plt.savefig("rand_walk"+str(n)+".png",bbox_inches="tight",dpi=600)
    plt.show()
plotdist(length)
def plotvar(X):
    '''Plot the variance of the lengths between points number of steps apart
    Parameters
    -----------

    X :A 2d numpy array for a walk'''
    variance=np.array([])
    mean=([])
    nsteps=np.logspace(1/2,3/2,40,dtype=np.int64)
    
    for i in nsteps:

        length=lengthfind(i,X)


        mean=np.append(mean,np.mean(np.square(length)))
        variance=np.append(variance,np.var(length))
    slope, intercept, r_value, p_value, std_err = linregress(np.log(nsteps), np.log(variance))
    fig=plt.figure()
    ax = fig.add_subplot(1, 1, 1) 
    ax.set_yscale('log')
    ax.set_ylabel('Variance of length')
    ax.set_xlabel('Step gap')
    ax.set_xscale('log')
    ax.set_title('Plot of the Variance Length aginst step gap')
    ax.plot(nsteps,variance, label='Results')
    ax.plot(nsteps,mean, label='Results 2')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope), np.exp(intercept)), label='linear regression')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope+std_err), np.exp(intercept)),'.k',label='upper error')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope-std_err), np.exp(intercept)),'.k', label='lower error')
    ax.plot(nsteps,np.multiply(np.power(nsteps, 3/2), np.exp(intercept)), label= 'General esitmate')
    ax.legend()
    plt.show()
    print(str(slope) + ' with error ' +str(std_err))

def findradius(X):
    '''Finds the minimum radius of a sphere to contain the random walk
    Parameters
    -----------
    X :A 2d numpy array for a walk
    Returns
    --------
    Radius: A float of the minium radius'''
    maxlength=np.array([])
    N=np.size(X,0)

    nsteps=np.linspace(0,N-1,N,dtype=np.int64)

    for i in nsteps:

        length=lengthfind(i,X)
        maxlength=np.append(maxlength,np.max(length))
    radius= np.max(maxlength)/2
    return radius
plotvar(x)
def plotmassvol():
    '''Plots the mass against volume containing walks'''
    r=np.array([])
    N=np.logspace(1,2,40,dtype=np.int64)
    for i in N:
        x=np.int64(createwalk(i))
        r=np.append(r,findradius(x))
    
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Volume')
    ax.set_xlabel('Mass')
    ax.set_title('Plot of mass against volume of a random walk')
    slope, intercept, r_value, p_value, std_err = linregress(np.log(N), np.log(r))
    ax.plot(N,r)
    ax.plot(N,np.multiply(np.power(N,slope),np.exp(intercept)), label='Results')
    ax.plot(N,np.multiply(np.power(N,slope+std_err),np.exp(intercept)),'.k',label='upper error')
    ax.plot(N,np.multiply(np.power(N,3/4),np.exp(intercept)), label= 'General esitmate')
    ax.plot(N,np.multiply(np.power(N,slope-std_err),np.exp(intercept)), '.k', label='lower error')
    ax.legend()
    plt.show()
    print(str(slope) + ' with error ' +str(std_err))
plotmassvol()
def spherecount(x):

    count=np.array([])
    length=np.sqrt(np.square(x[:,0])+np.square(x[:,1]))
    r=np.logspace(1/2,2,100)
    for i in r :
        count=np.append(count,np.size(length[length<i]))
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Volume')
    ax.set_xlabel('Mass')
    ax.set_title('Plot of mass against volume of a random walk')
    slope, intercept, r_value, p_value, std_err = linregress(np.log(count), np.log(r))
    ax.plot(count,r)
    ax.plot(count,np.multiply(np.power(count,slope),np.exp(intercept)), label='Results')
    ax.plot(count,np.multiply(np.power(count,slope+std_err),np.exp(intercept)),'.k',label='upper error')
    ax.plot(count,np.multiply(np.power(count,3/4),np.exp(intercept)), label= 'General esitmate')
    ax.plot(count,np.multiply(np.power(count,slope-std_err),np.exp(intercept)), '.k', label='lower error')
    ax.legend()
    plt.show()
    print(str(slope) + ' with error ' +str(std_err))
spherecount(x)
def timewalk(N):
   times=np.array([])
   Nsteps=np.logspace(1,np.log10(N), 50, dtype=np.int64)
   for i in Nsteps:
       print(i)
       times=np.append(times,timeit.timeit('createwalk('+ str(i) +')', setup="from __main__ import createwalk", number=5))
   fig = plt.figure()
   ax = fig.add_subplot(1, 1, 1)
   slope, intercept, r_value, p_value, std_err = linregress(np.log(Nsteps), np.log(times))
   ax.set_yscale('log')
   ax.set_xscale('log')
   ax.set_ylabel('Time taken')
   ax.set_xlabel('Number of steps')
   ax.set_title('Time taken to run exact self avoiding walk')
   ax.plot(Nsteps,times,label='Results')
   ax.plot(Nsteps,np.multiply(np.power(Nsteps,slope),np.exp(intercept)), label='Regression with slope {:.4f} with error {:.4f} and intercept {:.4f}'.format(slope,std_err,intercept))
   ax.plot(Nsteps,np.multiply(np.power(Nsteps,slope+std_err),np.exp(intercept)),'.k',label='upper error')
   ax.plot(Nsteps,np.multiply(np.power(Nsteps,slope-std_err),np.exp(intercept)), '.k', label='lower error')
   ax.legend()
   plt.savefig('Time taken to run exact self avoiding walk')   
   plt.show()
   print(str(slope) + ' with error ' +str(std_err))
   print('intercept '+str(intercept))
timewalk(500)
