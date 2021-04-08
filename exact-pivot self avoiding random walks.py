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

def createwalk(N):
    ''' Creates a self avoiding random walk using a mix of the pivot algorithm and an exact method and gives a 2d numpy array
    Parameters
    -----------
    N: integer, which denotes the length of the random walk
    
    Returns
    --------
    chain: A 2d numpy array, with the coords of the random walk
    '''
    i=1

    chain=np.array([[0,0]])

    while i<N :
    

        val = random.randint(1, 4)
        if val == 1 and not (np.any(chain[chain[:,0]==chain[i-1,0]+1,1]==chain[i-1,1])):
            chain=np.append(chain,np.array([[chain[i - 1,0] + 1,chain[i - 1,1]]]),axis=0) 
            
            i+=1
        elif val == 2 and not (np.any(chain[chain[:,0]==chain[i-1,0]-1,1]==chain[i-1,1])):
            chain=np.append(chain,np.array([[chain[i - 1,0] - 1,chain[i - 1,1]]]),axis=0)
            i+=1
        elif val == 3 and not (np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]+1)):
            chain=np.append(chain,np.array([[chain[i - 1,0] ,chain[i - 1,1]+1]]),axis=0)
            i+=1
        elif val==4 and not (np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]-1)):
            chain=np.append(chain,np.array([[chain[i - 1,0] ,chain[i - 1,1]-1]]),axis=0)
            i+=1
        elif np.any(chain[chain[:,0]==chain[i-1,0]+1,1]==chain[i-1,1]) and np.any(chain[chain[:,0]==chain[i-1,0]-1,1]==chain[i-1,1]) and np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]+1) and np.any(chain[chain[:,0]==chain[i-1,0],1]==chain[i-1,1]-1):
            m=len(chain)
            side=random.randint(1,2)
            #picks side to rotate
            pivot_site=random.randint(0, m-1)
            #picks a site to pivot a around
            rot=random.randint(1,3)
            #picks a rotation(90,180,270)
            if side==1:
                Nchain=chain[pivot_site:,]-chain[pivot_site]
        
                #selects chain to rotate
                if rot==1:
                    #array for rotation
                    rotation=np.array([[0,-1],[1,0]])
                    for j in range(len(Nchain)):
                        Nchain[j]=np.matmul(rotation,Nchain[j])+chain[pivot_site]
                elif rot==2:
                    rotation=np.array([[-1,0],[0,-1]])
                    #array for rotation
                    for j in range(len(Nchain)):
                        Nchain[j]=np.matmul(rotation,Nchain[j])+chain[pivot_site]
                elif rot==3:
                    #array for rotation
                    rotation=np.array([[0,1],[-1,0]])
                    for j in range(len(Nchain)):
                       
                        Nchain[j]=np.matmul(rotation,Nchain[j])+chain[pivot_site]
                overlap = cdist(Nchain,chain)
                    
                overlap = overlap.flatten()
                if len(np.nonzero(overlap)[0]) == len(overlap):
                    #checks if rotated chain overlaps with the rest
                    print(Nchain)
                    chain[pivot_site:,]=Nchain
        
                    chain=chain-chain[0]
            elif side==2:
                Nchain=chain[:pivot_site,]-chain[pivot_site]
        
                #selects chain to rotate
                if rot==1:
                    #array for rotation
                    rotation=np.array([[0,-1],[1,0]])
                    for j in range(len(Nchain)):
                        Nchain[j]=np.matmul(rotation,Nchain[j])+chain[pivot_site]
                elif rot==2:
                    rotation=np.array([[-1,0],[0,-1]])
                    #array for rotation
                    for j in range(len(Nchain)):
                        Nchain[j]=np.matmul(rotation,Nchain[j])+chain[pivot_site]
                elif rot==3:
                    #array for rotation
                    rotation=np.array([[0,1],[-1,0]])
                    for j in range(len(Nchain)):
                       
                        Nchain[j]=np.matmul(rotation,Nchain[j])+chain[pivot_site]
                overlap = cdist(Nchain,chain)
                    
                overlap = overlap.flatten()
                if len(np.nonzero(overlap)[0]) == len(overlap):
                    #checks if rotated chain overlaps with the rest
    
                    chain[:pivot_site,]=Nchain
        
                    chain=chain-chain[0]
       
    
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
        choice=X[i]-X[i-1]

    
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

def lengthfind(L,X):
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
    for i in range(L,(np.size(X,0))):
        lengthx.append(X[i,0]-X[i-L,0])
        lengthy.append(X[i,1]-X[i-L,1])
    
    length=np.sqrt(np.square(np.array(lengthy))+np.square(np.array(lengthx)))
    np.savetxt('length.txt',length)
    # save numpy array to a text doucment
    return length


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
    ax.set_title('Distribution of lengths of a step gap of '+str(l) +' for a SAW walk')
    #plotting stuff:
    plt.savefig('Distribution of lengths of a step gap of '+str(l) +' for a SAW walk')
    plt.show()
    print(str(mu) + ' with error ' +str(std))

def plotvar(X):
    '''Plot the variance of the lengths between points number of steps apart
    Parameters
    -----------

    X :A 2d numpy array for a walk'''
    variance=[]
    nsteps=np.logspace(1/2,3/2,40,dtype=np.int64)
    
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
    ax.set_title('Plot of the Variance Length aginst step gap for a SAW')
    ax.plot(nsteps,variance, label='Results')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope), np.exp(intercept)), label='linear regression')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope+std_err), np.exp(intercept)),'.k',label='upper error')
    ax.plot(nsteps,np.multiply(np.power(nsteps, slope-std_err), np.exp(intercept)),'.k', label='lower error')
    ax.plot(nsteps,np.multiply(np.power(nsteps,3/2), np.exp(intercept)), label= 'General esitmate')
    ax.legend()
    plt.savefig('Plot of the Variance Length against step gap for SAW')
    plt.show()
    print(str(slope) + ' with error ' +str(std_err)+' intercept ' +str(intercept))

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

def plotmassvol(n):
    '''Plots the number of steps against the radius need to contain the walks
    Parameters
    ------------
    n is an integer, which denotes the upper limit of the stepsize to be of walk to be checked'''
    r=[]
    N=np.logspace(1,np.log10(n),40,dtype=np.int64)
    for i in N:
        x=np.int64(createwalk(i))
        r.append(findradius(x))
    
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Volume')
    ax.set_xlabel('Mass')
    ax.set_title('Plot of mass of the total SAW walk against radius of to contain it')
    slope, intercept, r_value, p_value, std_err = linregress(np.log(N), np.log(r))
    ax.plot(N,r)
    ax.plot(N,np.multiply(np.power(N,slope),np.exp(intercept)), label='Regression')
    ax.plot(N,np.multiply(np.power(N,slope+std_err),np.exp(intercept)),'.k',label='upper error')
    ax.plot(N,np.multiply(np.power(N,3/4),np.exp(intercept)), label= 'General esitmate')
    ax.plot(N,np.multiply(np.power(N,slope-std_err),np.exp(intercept)), '.k', label='lower error')
    ax.legend()
    plt.savefig('Plot of mass of the total SAW walk against radius of to contain it')
    plt.show()
    print(str(slope) + ' with error ' +str(std_err))

def spherecount(x):
    '''Plots the radius from the origin against the number of steps contained within
    Parameters
    ----------
    x is a 2d numpy array representing the walk to be checked'''
    count=[]
    length=np.sqrt(np.square(x[:,0])+np.square(x[:,1]))
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
    ax.set_title('Plot of mass against volume of a SAW walk')
    slope, intercept, r_value, p_value, std_err = linregress(np.log(count), np.log(r))

   
def findslope(x):
    '''Finds the slope, but doesn't plot the graph of spherecount
    Parameters
    ----------
    x:a numpy array representing the walk
    Returns
    --------
    slope: the slope of the logs of the radius and the number of steps conatined within'''
    count=np.array([])
    length=np.sqrt(np.square(x[:,0])+np.square(x[:,1]))
    r=np.logspace(1/2,2,100)
    for i in r :
        count=np.append(count,np.size(length[length<i]))
    

    slope, intercept, r_value, p_value, std_err = linregress(np.log(count), np.log(r))

    return slope
def meanslope(Samples,N):
    '''Plots the distrubution of a slopes fromfrom the function spherecount
    Parameters
    ------------
    Samples is an integer representing the number of samples to be taken form the function
    
    N is an integer for the size of the create walk being tested
    '''
    grad=[]
    for i in range(Samples):
        X=createwalk(N)
        grad.append(findslope(X))
    nr_in_bin, bin_edges=np.histogram(grad)
    width = bin_edges[1:] - bin_edges[:-1]
    center = (bin_edges[:-1] + bin_edges[1:]) / 2   
                           
    s = 1/len(grad)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(center, nr_in_bin*s, align='center', width=width*0.9, 
           color='g', label='Histogram of slopes')
    mu, std = norm.fit(grad)
    p=norm.pdf(center, mu, std)
    ax.plot(center,p, label='Gaussian of slopes')
    ax.set_xlabel('slopes')
    ax.set_ylabel('Probablity density')
    ax.legend()
    ax.set_title('Distribution of the slopes upto' +str(n)+' steps for a SAW')
    #plotting stuff:
    plt.savefig('Distribution of the slopes upto'+str(n)+' steps for a SAW')
    plt.show()
    print(str(mu) + ' with error ' +str(std))

def timewalk(N):
    '''This times the function create walk for different steps up to N and plots them
    Parameters 
    -----------
    N is an integer of the limit of the steps.'''
    times=np.array([])
    Nsteps=np.logspace(1,np.log10(N), 50, dtype=np.int64)
    for i in Nsteps:
       
        times=np.append(times,timeit.timeit('createwalk('+ str(i) +')', setup="from __main__ import createwalk", number=5))
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    slope, intercept, r_value, p_value, std_err = linregress(np.log(Nsteps), np.log(times))
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Time taken')
    ax.set_xlabel('Number of steps')
    ax.set_title('Time taken to run mixed self avoiding walk')
    ax.plot(Nsteps,times,label='Results')
    ax.plot(Nsteps,np.multiply(np.power(Nsteps,slope),np.exp(intercept)), label='Regression with slope {:.4f} with error {:.4f} and intercept {:.4f}'.format(slope,std_err,intercept))
    ax.plot(Nsteps,np.multiply(np.power(Nsteps,slope+std_err),np.exp(intercept)),'.k',label='upper error')
    ax.plot(Nsteps,np.multiply(np.power(Nsteps,slope-std_err),np.exp(intercept)), '.k', label='lower error')
    ax.legend()
    plt.savefig('Time taken to run mixed self avoiding walk')   
    plt.show()
    
    print(str(slope) + ' with error ' +str(std_err))
    print('intercept '+str(intercept))
run=True
while run==True:
    n=int(input('What size of self avoiding walk do you want?'))
    x=createwalk(n)
    fig=plt.figure()
    ax = fig.add_subplot(1, 1, 1)        
    ax.set_title("Random Walk ($n = " + str(n) + "$ steps)")
    ax.plot(x[:,0],x[:,1],'r' )
    plt.savefig("rand_walk"+str(n)+".png",bbox_inches="tight",dpi=600)
    plt.show()
    print('Menu')
    print('-------------')
    print('Print walk[1]')
    print('Plot distribution of lengths for a certain step gap[2]')
    print('Plot the variance of the lengths between points number of steps[3]')
    print('Plots the number of steps against the radius needed to contain the walk[4]')
    print('Plots the radius from the origin against the number of steps contained within[5]')
    print('plots the distribution of the slope of the radius from the origin against the number of steps contained within[6]')
    print('Plot the time to create a walk against the number of steps[7]')
    print('Exit[8]')
    choice=int(input('Choose:'))
    if choice==1:
        printgrid(x)
    elif choice==2:
        l=int(input('What step gap do you want:'))
        length=lengthfind(l,x)
        plotdist(length)
    elif choice==3:
        plotvar(x) 
    elif choice==4:
        plotmassvol(n)
    elif choice==5:
        spherecount(x)
    elif choice==6:
        samples=int(input('How many samples do you want:'))
        meanslope(samples,n)
    elif choice==7:
        timewalk(n)
    elif choice==8:
        run=False

