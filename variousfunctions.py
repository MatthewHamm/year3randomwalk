# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 15:21:44 2019

@author: Matthew Hamm
This contains other functions ,which were subbed in for report
"""
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import norm
from scipy.spatial.distance import cdist
def pivotcreatewalk(N):
    ''' Creates a self avoiding random walk using the pivot algorithm and gives a 2d numpy array
    Parameters
    -----------
    N: integer, which denotes the length of the random walk
    
    Returns
    --------
    chain_init: A 2d numpy array, with the coords of th random walk
    '''
    t=3/4*N
    #number of rotations
    curr_rot=0
    #current number of rotations
    chain_init=np.zeros([N,2])
    chain_init[:,0]=np.arange(N)
    #creates an array, which represents a straight line
    
    while curr_rot<t:
    #loop untill number of raotation is reached

        pivot_site=random.randint(0, N-1)
        #picks a site to pivot a around
        rot=random.randint(1,3)
        #picks a rotation(90,180,270)
        chain=None
        #empties roated chain
        side=random.randint(1,2)
        #picks side to rotate
        if side==1:
            chain=chain_init[pivot_site:,]-chain_init[pivot_site]
            #selects chain to rotate
            if rot==1:
                #array for rotation
                rotation=np.array([[0,-1],[1,0]])
                for i in range(len(chain)):
                    chain[i]=np.matmul(rotation,chain[i])+chain_init[pivot_site]
            elif rot==2:
                rotation=np.array([[-1,0],[0,-1]])
                #array for rotation
                for i in range(len(chain)):
                    chain[i]=np.matmul(rotation,chain[i])+chain_init[pivot_site]
            elif rot==3:
                #array for rotation
                rotation=np.array([[0,1],[-1,0]])
                for i in range(len(chain)):
               
                    chain[i]=np.matmul(rotation,chain[i])+chain_init[pivot_site]
            overlap = cdist(chain,chain_init)
            
            overlap = overlap.flatten()
            if len(np.nonzero(overlap)[0]) == len(overlap):
                #checks if rotated chain overlaps with the rest
                chain_init[pivot_site:,]=chain
                chain_init=chain_init-chain_init[0]
                curr_rot+=1
        elif side==2:
            chain=chain_init[:pivot_site,]-chain_init[pivot_site]
        
            if rot==1:
                rotation=np.array([[0,-1],[1,0]])
                for i in range(len(chain)):
                    chain[i]=np.matmul(rotation,chain[i])+chain_init[pivot_site]
            elif rot==2:
                rotation=np.array([[-1,0],[0,-1]])
                for i in range(len(chain)):
                    chain[i]=np.matmul(rotation,chain[i])+chain_init[pivot_site]
            elif rot==3:
                rotation=np.array([[0,1],[-1,0]])
                for i in range(len(chain)):
               
                    chain[i]=np.matmul(rotation,chain[i])+chain_init[pivot_site]
    
            overlap = cdist(chain,chain_init)
            overlap = overlap.flatten()
            if len(np.nonzero(overlap)[0]) == len(overlap):
                chain_init[:pivot_site,]=chain        
                chain_init=chain_init-chain_init[0]
                curr_rot+=1
    
    return chain_init
def exactcreatewalk(N):
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

def 1Dcreatewalk(N):
    ''' Creates a random walk and gives a 1d numpy array
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
def 3Dcreatewalk(N):
    ''' Creates a random walk and gives a 3d numpy array
    Parameters
    -----------
    N: integer, which denotes the length of the random walk
    
    Returns
    --------
    chain: A 3d numpy array, with the coords of the random walk
    '''
    chain=np.zeros([N,3])
    
    for i in range(1, N):
        val = random.randint(1, 6)
        if val == 1:
            chain[i,0] = chain[i - 1,0] + 1
            chain[i,1] = chain[i - 1,1]
            chain[i,2] =chain[i-1,2]
        elif val == 2:
            chain[i,0] = chain[i - 1,0] - 1
            chain[i,1] = chain[i - 1,1]
            chain[i,2] = chain[i-1,2]
        elif val == 3:
            chain[i,0] = chain[i - 1,0] 
            chain[i,1] = chain[i - 1,1] +1
            chain[i,2] = chain[i-1,2]
        elif val==4:
            chain[i,0] = chain[i - 1,0]
            chain[i,1] = chain[i - 1,1]-1
            chain[i,2] = chain[i-1,2]
        elif val==5:
            chain[i,0] = chain[i - 1,0]
            chain[i,1] = chain[i - 1,1]
            chain[i,2] =chain[i-1,2]+1   
        else :
            chain[i,0] = chain[i - 1,0]
            chain[i,1] = chain[i - 1,1]
            chain[i,2] = chain[i-1,2]-1   
    
        
    return chain
def 3Dlengthfind(L,X):
    '''Finds the length between points a certain number of step apart
    Parameters
    -----------
    L : an integer for the step gap
    X:A 2d numpy array for a walk
    Returns
    --------
    Length:A numpy array of lengths a certain number of steps apart'''
    lengthy=[]
    lengthx=[]
    lengthz=[]
    for i in range(L,(np.size(X,0))):
        lengthx.append(X[i,0]-X[i-L,0])
        lengthy.append(X[i,1]-X[i-L,1])
        lengthz.append(X[i,2]-X[i-L,2])
    
    length=np.sqrt(np.square(np.array(lengthy))+np.square(np.array(lengthx))+np.square(np.array(lengthz)))
    np.savetxt('length.txt',length)
    # save numpy array to a text doucment
    return length
def 3Dspherecount(x):
    '''Plots the radius ofrom the origin sand the number of stepss contained within
    Parameters
    ----------
    x is a 3d numpy array representing the walk to be checked'''
    count=[]
    length=np.sqrt(np.square(x[:,0])+np.square(x[:,1])+np.square(x[:,2]))
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
def 1Dlengthfind(L,X):
    '''Finds the length between points a certain number of step apart
    Parameters
    -----------
    L : an integer for the step gap
    X:A 1d numpy array for a walk
    Returns
    --------
    Length:A numpy array of lengths a certain number of steps apart'''
    length=[]

    for i in range(L,(np.size(X))):
        length.append(X[i]-X[i-L])

    
    length=np.array(length)
    np.savetxt('length.txt',length)
    # save numpy array to a text doucment
    return length
def 1Dspherecount(x):
    '''Plots the radius ofrom the origin sand the number of stepss contained within
    Parameters
    ----------
    x is a 1d numpy array representing the walk to be checked'''
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
