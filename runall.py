# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 22:33:21 2021

@author: Mein PC
"""

import numpy as np
import scipy as sp
import math

import matrizen as ma
    
    
def initex2_40():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 80
    faces = int(n/2)
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 1.5
    
    arad=0.5
    brad=0.5
    
    filename='ex2_40.txt'
    
    m=int(2*3*faces)
    
    vi = np.zeros((n,2), dtype=np.float)
    for i in range(0,n):
        ti = 2*sp.pi*i/n
        
        #Circle and Ellipse
        vi[i,0] = arad*math.cos(ti)
        vi[i,1] = brad*math.sin(ti)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    

def initex3_40():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 80
    faces = int(n/2)
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 2.1
    
    arad=0.5
    brad=0.5
    
    filename='ex3_40.txt'
    
    m=int(2*3*faces)
    
    vi = np.zeros((n,2), dtype=np.float)
    for i in range(0,n):
        ti = 2*sp.pi*i/n
        
        #Circle and Ellipse
        vi[i,0] = arad*math.cos(ti)
        vi[i,1] = brad*math.sin(ti)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def initex2_40e():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 80
    faces = int(n/2)
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 1.4
    
    arad=1
    brad=0.5
    
    filename='ex2_40e.txt'
    
    m=int(2*3*faces)
    
    vi = np.zeros((n,2), dtype=np.float)
    for i in range(0,n):
        ti = 2*sp.pi*i/n
        
        #Circle and Ellipse
        vi[i,0] = arad*math.cos(ti)
        vi[i,1] = brad*math.sin(ti)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def initex2_40k():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 80
    faces = int(n/2)
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 0.9
    
    filename='ex2_40k.txt'
    
    m=int(2*3*faces)
    
    vi = np.zeros((n,2), dtype=np.float)
    for i in range(0,n):
        ti = 2*sp.pi*i/n
        
        #Kite
        vi[i,0] = 0.75*math.cos(ti) + 0.3*math.cos(2*ti)
        vi[i,1] = math.sin(ti)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def initex3_40k():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 80
    faces = int(n/2)
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 1.1
    
    filename='ex3_40k.txt'
    
    m=int(2*3*faces)
    
    vi = np.zeros((n,2), dtype=np.float)
    for i in range(0,n):
        ti = 2*sp.pi*i/n
        
        #Kite
        vi[i,0] = 0.75*math.cos(ti) + 0.3*math.cos(2*ti)
        vi[i,1] = math.sin(ti)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    

def initex2_46q():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 24
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 1.5
    
    filename='ex2_46q.txt'
        
    si=np.linspace(0,1,n)
    v1=np.array([[0], [0]])
    v2=np.array([[1], [0]])
    v3=np.array([[1], [1]])
    v4=np.array([[0], [1]])
    vi1=v1*(1-si)+v2*si
    vi2=v2*(1-si)+v3*si
    vi3=v3*(1-si)+v4*si
    vi4=v4*(1-si)+v1*si
    vi=np.concatenate((vi1[:,:-1], vi2[:,:-1], vi3[:,:-1], vi4), axis=1).T
    n=vi.shape[0]-1
    
    faces = int(n/2)
    
    m=int(2*3*faces)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def initex3_46q():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 24
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 1.8
    
    filename='ex3_46q.txt'
        
    si=np.linspace(0,1,n)
    v1=np.array([[0], [0]])
    v2=np.array([[1], [0]])
    v3=np.array([[1], [1]])
    v4=np.array([[0], [1]])
    vi1=v1*(1-si)+v2*si
    vi2=v2*(1-si)+v3*si
    vi3=v3*(1-si)+v4*si
    vi4=v4*(1-si)+v1*si
    vi=np.concatenate((vi1[:,:-1], vi2[:,:-1], vi3[:,:-1], vi4), axis=1).T
    n=vi.shape[0]-1
    
    faces = int(n/2)
    
    m=int(2*3*faces)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def initex2_40C():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 80
    faces = int(n/2)
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 2+0.5j
    
    arad=0.5
    brad=0.5
    
    filename='ex2_20C.txt'
    
    m=int(2*3*faces)
    
    vi = np.zeros((n,2), dtype=np.float)
    for i in range(0,n):
        ti = 2*sp.pi*i/n
        
        #Circle and Ellipse
        vi[i,0] = arad*math.cos(ti)
        vi[i,1] = brad*math.sin(ti)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def initex2_46qC():
    
    mu    = 1/16
    lam   = 1/4
    sigma1= 1
    sigma2= 4
    alpha = (1-math.sqrt(3/5))/2
    n     = 24
    
    N = 24
    l=20
    tol_rank = 10**(-2)
    R  = 0.25
    mp = 2+0.5j
    
    filename='ex2_46qC.txt'
        
    si=np.linspace(0,1,n)
    v1=np.array([[0], [0]])
    v2=np.array([[1], [0]])
    v3=np.array([[1], [1]])
    v4=np.array([[0], [1]])
    vi1=v1*(1-si)+v2*si
    vi2=v2*(1-si)+v3*si
    vi3=v3*(1-si)+v4*si
    vi4=v4*(1-si)+v1*si
    vi=np.concatenate((vi1[:,:-1], vi2[:,:-1], vi3[:,:-1], vi4), axis=1).T
    n=vi.shape[0]-1
    
    faces = int(n/2)
    
    m=int(2*3*faces)
        
    # create collocation nodes vh_i
    vhi = np.zeros((3*faces,2), dtype=np.float)
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = ma.index(k,n-1)
        ind3 = ma.index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+1,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = ma.lagrange1(s)*x1+ma.lagrange2(s)*x2+ma.lagrange3(s)*x3
        vhi[L+2,1] = ma.lagrange1(s)*y1+ma.lagrange2(s)*y2+ma.lagrange3(s)*y3
        k = k+2
        L = L+3
    
    ma.initBeyn(mu,lam,sigma1,sigma2,alpha,faces,n,vi,vhi,N,l,m,tol_rank,R,mp,filename)
    ma.beyn()
    
    
def runAllPaper():
    
    #Example 1
    error = np.array([ma.testD(10, 1),ma.testD(20, 1),ma.testD(40, 1),ma.testD(80, 1)])
    eoc = np.array([0.0,0.0,0.0])
    eoc[0] = (1.0 / math.log(2)) * math.log(error[0] / error[1])
    eoc[1] = (1.0 / math.log(2)) * math.log(error[1] / error[2])
    eoc[2] = (1.0 / math.log(2)) * math.log(error[2] / error[3])
    #with open("testDL_1.txt", "w") as external_file:
    print("test DL, omega=1")
    print(10, error[0])
    print(20, error[1],eoc[0])
    print(40, error[2],eoc[1])
    print(80, error[3],eoc[2])
    #    external_file.close()
    
    error = np.array([ma.testD(10, 1j),ma.testD(20, 1j),ma.testD(40, 1j),ma.testD(80, 1j)])
    eoc = np.array([0.0,0.0,0.0])
    eoc[0] = (1.0 / math.log(2)) * math.log(error[0] / error[1])
    eoc[1] = (1.0 / math.log(2)) * math.log(error[1] / error[2])
    eoc[2] = (1.0 / math.log(2)) * math.log(error[2] / error[3])
    #with open("testDL_1i.txt", "w") as external_file:
    print("test DL, omega=i")
    print(10, error[0])
    print(20, error[1],eoc[0])
    print(40, error[2],eoc[1])
    print(80, error[3],eoc[2])
    #    external_file.close()
    
    #Example 2
    error = np.array([ma.testS(10, 1),ma.testS(20, 1),ma.testS(40, 1),ma.testS(80, 1)])
    eoc = np.array([0.0,0.0,0.0])
    eoc[0] = (1.0 / math.log(2)) * math.log(error[0] / error[1])
    eoc[1] = (1.0 / math.log(2)) * math.log(error[1] / error[2])
    eoc[2] = (1.0 / math.log(2)) * math.log(error[2] / error[3])
    #with open("testSL_1.txt", "w") as external_file:
    print("test SL, omega=1")
    print(10, error[0])
    print(20, error[1],eoc[0])
    print(40, error[2],eoc[1])
    print(80, error[3],eoc[2])
    #    external_file.close()
    
    error = np.array([ma.testS(10, 1j),ma.testS(20, 1j),ma.testS(40, 1j),ma.testS(80, 1j)])
    eoc = np.array([0.0,0.0,0.0])
    eoc[0] = (1.0 / math.log(2)) * math.log(error[0] / error[1])
    eoc[1] = (1.0 / math.log(2)) * math.log(error[1] / error[2])
    eoc[2] = (1.0 / math.log(2)) * math.log(error[2] / error[3])
    #with open("testSL_1i.txt", "w") as external_file:
    print("test SL, omega=i")
    print(10, error[0])
    print(20, error[1],eoc[0])
    print(40, error[2],eoc[1])
    print(80, error[3],eoc[2])
    #    external_file.close()
        
    #D1: Circle
    print("Circle")
    initex2_40()
    initex3_40()
    
    #D2: Ellipse
    print("Ellipse")
    initex2_40e()
    
    #D3: Kite 
    print("Kite")
    initex2_40k()
    initex3_40k()
    
    #D4: Square
    print("Square")
    initex2_46q()
    initex3_46q()
    
    #Complex
    print("Complex")
    initex2_40C()
    initex2_46qC()
    
    
if __name__ == "__main__":

    runAllPaper()
