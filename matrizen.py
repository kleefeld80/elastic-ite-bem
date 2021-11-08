# -*- coding: utf-8 -*-
"""

"""

import sys
import scipy as sp
import numpy as np
import numpy.matlib
import cmath
import math
from scipy.integrate import quad
from scipy.special import hankel1
from concurrent import futures

import warnings
warnings.filterwarnings('ignore')

def beyn():
    global l
    
    print("faces: ",faces)
    
    #Kontur diskretisieren
    phi = np.zeros(N, dtype=np.complex)
    dphi = np.zeros(N, dtype=np.complex)
    for i in range(0,N):
        ti = 2*sp.pi*i/N
        phi[i] = mp + R*np.exp(1j*ti)
        dphi[i] = R*1j*np.exp(1j*ti)
    
    while (l<=m):
        print(l)
    
        # 2. zufaellige Matrix Vdach erzeugen
        Vdach=np.matlib.rand(m,l) + np.matlib.rand(m,l) * 1j
        #print(Vdach)
        
        # 3. A_{0,N} und A-{1,N} berechnen
        A0 = np.zeros((m, l), dtype=np.complex)
        A1 = np.zeros((m, l), dtype=np.complex)
                   
        if __name__ == "__main__":
            it=list(range(0,N))
            with futures.ProcessPoolExecutor(max_workers=4) as e:
                fs = {e.submit(matr, phi[i]): i for i in it}
                for f in futures.as_completed(fs):
                    i=fs[f]
                    T=f.result()
                    A0 += np.linalg.solve(T,(Vdach*dphi[i]))
                    A1 += phi[i]*np.linalg.solve(T,(Vdach*dphi[i]))
                    
        for i in range(0,N):
            T = matr(phi[i])
            A0 += np.linalg.solve(T,(Vdach*dphi[i]))
            A1 += phi[i]*np.linalg.solve(T,(Vdach*dphi[i]))
         
         
        A0 = A0/(1j*N)
        A1 = A1/(1j*N)
    
        # 4. SVD von A0 berechnen
        [V, s, WH] = np.linalg.svd(A0, full_matrices=True)
        #print(V)
        W = WH.conj().T
    
        #5. Rang Test
        # for ind in range(0,len(s)):
        #     if s[ind] >= tol:
        #         k=ind+1
        k = np.sum(s > tol_rank)
        if l == k and l != m:
            l += 5
        else:
            break
  
    # 5. k < l => bestimme V0, W0 und S0
    V0 = V[0:m,0:k]
    W0 = W[0:l,0:k]
    S0inv = np.diag(1/s)[0:k,0:k]
    V0H = V0.conj().T
  
    # 6. B berechnen
    B = np.matmul(np.matmul(np.matmul(V0H,A1),W0), S0inv)
        
    # 7. Eigenwerte s von B, Eigenvektoren V
    sb, V = np.linalg.eig(B)
    
    printeigs(sb)
    return 

def matr(omegas):
    global omega
    
    omega=omegas*math.sqrt(sigma2)
    
    Sk = createSmatrix()
    A  = createDmatrix()
    Dk = 0.5*np.identity(A.shape[0])+A
    
    omega = omegas*math.sqrt(sigma1)
    
    Skn = createSmatrix()
    A   = createDmatrix()
    Dkn = 0.5*np.identity(A.shape[0])+A
    
    A = np.linalg.solve(Skn,Dkn) - np.linalg.solve(Sk,Dk)
    return A    

def createDmatrix():
    
    A = np.zeros((2*3*faces,2*3*faces), dtype=np.complex)
    M0 = np.zeros((2*3*faces,2*3*faces), dtype=np.complex)
    
    for i in range(0,3*faces):
        k = 0
        L = 0
        xh = vhi[i,0]
        yh = vhi[i,1]
        for j in range(0,faces):
            ind1 = k
            ind2 = index(k,n-1)
            ind3 = index(k+1,n-1)
            x1 = vi[ind1,0]
            x2 = vi[ind2,0]
            x3 = vi[ind3,0]
            y1 = vi[ind1,1]
            y2 = vi[ind2,1]
            y3 = vi[ind3,1]
            k=k+2
            
            if i==L:
                E1 = np.zeros((2,2), dtype=np.complex)
                E1c = np.zeros((2,2), dtype=np.complex)
            else:
                E1 = integrateD(x1, x2, x3, y1, y2, y3, xh, yh, 1)
                E1c = integrateD0(x1, x2, x3, y1, y2, y3, xh, yh, 1)
                
            if i==L+1:
                E2 = np.zeros((2,2), dtype=np.complex)
                E2c = np.zeros((2,2), dtype=np.complex)
            else:
                E2 = integrateD(x1, x2, x3, y1, y2, y3, xh, yh, 2)
                E2c = integrateD0(x1, x2, x3, y1, y2, y3, xh, yh, 2)
                
            if i==L+2:
                E3 = np.zeros((2,2), dtype=np.complex)
                E3c = np.zeros((2,2), dtype=np.complex)
            else:
                E3 = integrateD(x1, x2, x3, y1, y2, y3, xh, yh, 3)
                E3c = integrateD0(x1, x2, x3, y1, y2, y3, xh, yh, 3)
                
            M0[2*i:2*i+2,2*(L  ):2*L+2] = E1c
            M0[2*i:2*i+2,2*(L+1):2*(L+1)+2] = E2c
            M0[2*i:2*i+2,2*(L+2):2*(L+2)+2] = E3c
              
            A[2*i:2*i+2,2*(L  ):2*L+2] = E1
            A[2*i:2*i+2,2*(L+1):2*(L+1)+2] = E2
            A[2*i:2*i+2,2*(L+2):2*(L+2)+2] = E3
            
            L=L+3;
    
    #Trick -1/2*I
    #omega-0 seems to be zero on the diagonal
    for i in range(0,3*faces):
        A[2*i:2*i+2,2*i] = [-0.5,0] - np.sum(M0[2*i:2*i+2,0:2*3*faces+1:2],axis=1)
        A[2*i:2*i+2,2*i+1] = [0,-0.5] - np.sum(M0[2*i:2*i+2,1:2*3*faces+1:2],axis=1);
    return A


def createSmatrix():
    
    A = np.zeros((2*3*faces,2*3*faces), dtype=np.complex)
    
    for i in range(0,3*faces):
        k = 0
        L = 0
        xh = vhi[i,0]
        yh = vhi[i,1]
        for j in range(0,faces):
            ind1 = k
            ind2 = index(k,n-1)
            ind3 = index(k+1,n-1)
            x1 = vi[ind1,0]
            x2 = vi[ind2,0]
            x3 = vi[ind3,0]
            y1 = vi[ind1,1]
            y2 = vi[ind2,1]
            y3 = vi[ind3,1]
            k=k+2
            
            if i==L:
                E1 = integrateS(x1, x2, x3, y1, y2, y3, xh, yh, 1, 1)
            else:
                E1 = integrateS(x1, x2, x3, y1, y2, y3, xh, yh, 1, 0)
                
            if i==L+1:
                E2 = integrateS(x1, x2, x3, y1, y2, y3, xh, yh, 2, 2)
            else:
                E2 = integrateS(x1, x2, x3, y1, y2, y3, xh, yh, 2, 0)
                
            if i==L+2:
                E3 = integrateS(x1, x2, x3, y1, y2, y3, xh, yh, 3, 3)
            else:
                E3 = integrateS(x1, x2, x3, y1, y2, y3, xh, yh, 3, 0)
            
            A[2*i:2*i+2,2*(L  ):2*L+2] = E1
            A[2*i:2*i+2,2*(L+1):2*(L+1)+2] = E2
            A[2*i:2*i+2,2*(L+2):2*(L+2)+2] = E3
            
            L=L+3;
    
    return A

            
    
def integrateD(x1,x2,x3,y1,y2,y3,xh,yh,basis):
    eps=1e-6
    E = np.zeros((2,2), dtype=np.complex)
    E[0,0] = quad(D11, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D11, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    E[0,1] = quad(D12, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D12, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    E[1,0] = quad(D21, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D21, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    E[1,1] = quad(D22, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D22, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    
    return E


def integrateD0(x1,x2,x3,y1,y2,y3,xh,yh,basis):
    eps=1e-6
    E = np.zeros((2,2), dtype=np.complex)
    E[0,0] = quad(D11_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D11_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    E[1,0] = quad(D12_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D12_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    E[0,1] = quad(D21_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D21_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    E[1,1] = quad(D22_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(D22_0, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    
    return E


def integrateS(x1,x2,x3,y1,y2,y3,xh,yh,basis,singular):
    E = np.zeros((2,2), dtype=np.complex)
    eps=1e-6
    if singular==1:
        E[0,0] = quad(S11, 0, alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, 0, alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] = quad(S12, 0, alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, 0, alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,1] = quad(S22, 0, alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, 0, alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,0] += quad(S11, alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] += quad(S12, alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,1] += quad(S22, alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,0] = E[0,1]
    elif singular==2:
        E[0,0] = quad(S11, 0, 0.5-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, 0, 0.5-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] = quad(S12, 0, 0.5-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, 0, 0.5-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,1] = quad(S22, 0, 0.5-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, 0, 0.5-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,0] += quad(S11, 0.5+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, 0.5+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] += quad(S12, 0.5+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, 0.5+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,1] += quad(S22, 0.5+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, 0.5+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,0] = E[0,1]
    elif singular==3:
        E[0,0] = quad(S11, 0, 1-alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, 0, 1-alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] = quad(S12, 0, 1-alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, 0, 1-alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,1] = quad(S22, 0, 1-alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, 0, 1-alpha-eps, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,0] += quad(S11, 1-alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, 1-alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] += quad(S12, 1-alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, 1-alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,1] += quad(S22, 1-alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, 1-alpha+eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,0] = E[0,1]
    else:
        E[0,0] = quad(S11, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S11, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[0,1] = quad(S12, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S12, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
        E[1,0] = E[0,1]
        E[1,1] = quad(S22, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,0))[0] + 1j* quad(S22, eps, 1, args=(x1,x2,x3,y1,y2,y3,xh,yh,basis,1))[0]
    
    return E


def D11(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a=xh-((1-3*s+2*s**2)*x1+(4*s-4*s**2)*x2+(2*s**2-s)*x3)
    b=yh-((1-3*s+2*s**2)*y1+(4*s-4*s**2)*y2+(2*s**2-s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    phiS1d = -1j/(4*mu)*ks*H1ks-1j/(2*omega**2*r)*((kp*H1kp-ks*H1ks)/r+(ks**2*H0ks-kp**2*H0kp)/2)
    phiS2 = 1j/(2*omega**2*r)*(ks*H1ks-kp*H1kp)+1j/(4*omega**2)*(kp**2*H0kp-ks**2*H0ks)
    phiS2d =  1j/(2*omega**2)*(ks**2*H0ks/r-kp**2*H0kp/r-2*(ks*H1ks-kp*H1kp)/(r*r))+1j/(4*omega**2)*(ks**3*H1ks-kp**3*H1kp)
    U1 = (-lam-2*mu)*ny1*a-mu*ny2*b
    U1J = ((-lam-2*mu)*ny1*a**3-lam*ny1*a*b**2-2*mu*ny2*a**2*b)/(r**2)
    U2 = (-lam-4*mu)*ny1*a-mu*ny2*b+4*mu*(ny1*a**3+ny2*a**2*b)/(r**2)
    kernel = phiS1d/r*U1+phiS2d/r*U1J+phiS2/(r**2)*U2
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def D12(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a=xh-((1-3*s+2*s**2)*x1+(4*s-4*s**2)*x2+(2*s**2-s)*x3)
    b=yh-((1-3*s+2*s**2)*y1+(4*s-4*s**2)*y2+(2*s**2-s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    phiS1d = -1j/(4*mu)*ks*H1ks-1j/(2*omega**2*r)*((kp*H1kp-ks*H1ks)/r+(ks**2*H0ks-kp**2*H0kp)/2)
    phiS2 = 1j/(2*omega**2*r)*(ks*H1ks-kp*H1kp)+1j/(4*omega**2)*(kp**2*H0kp-ks**2*H0ks)
    phiS2d =  1j/(2*omega**2)*(ks**2*H0ks/r-kp**2*H0kp/r-2*(ks*H1ks-kp*H1kp)/(r*r))+1j/(4*omega**2)*(ks**3*H1ks-kp**3*H1kp)
    U1 = -1*lam*ny2*a-mu*ny1*b
    U1J = (-lam*ny2*a**3-2*mu*ny1*a**2*b+(-lam-2*mu)*ny2*a*b**2)/(r**2)
    U2 = (-lam-2*mu)*ny2*a-mu*ny1*b+4*mu*(ny1*a**2*b+ny2*a*b**2)/(r**2)
    kernel = phiS1d/r*U1+phiS2d/r*U1J+phiS2/(r**2)*U2
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag
    

def D21(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a=xh-((1-3*s+2*s**2)*x1+(4*s-4*s**2)*x2+(2*s**2-s)*x3)
    b=yh-((1-3*s+2*s**2)*y1+(4*s-4*s**2)*y2+(2*s**2-s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    phiS1d = -1j/(4*mu)*ks*H1ks-1j/(2*omega**2*r)*((kp*H1kp-ks*H1ks)/r+(ks**2*H0ks-kp**2*H0kp)/2)
    phiS2 = 1j/(2*omega**2*r)*(ks*H1ks-kp*H1kp)+1j/(4*omega**2)*(kp**2*H0kp-ks**2*H0ks)
    phiS2d =  1j/(2*omega**2)*(ks**2*H0ks/r-kp**2*H0kp/r-2*(ks*H1ks-kp*H1kp)/(r*r))+1j/(4*omega**2)*(ks**3*H1ks-kp**3*H1kp)
    U1 = -1*lam*ny1*b-mu*ny2*a
    U1J = ((-lam-2*mu)*ny1*a**2*b-lam*ny1*b**3-2*mu*ny2*a*b**2)/(r**2)
    U2 = (-lam-2*mu)*ny1*b-mu*ny2*a+4*mu*(ny1*a**2*b+ny2*a*b**2)/(r**2)
    kernel = phiS1d/r*U1+phiS2d/r*U1J+phiS2/(r**2)*U2
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def D22(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a=xh-((1-3*s+2*s**2)*x1+(4*s-4*s**2)*x2+(2*s**2-s)*x3)
    b=yh-((1-3*s+2*s**2)*y1+(4*s-4*s**2)*y2+(2*s**2-s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    phiS1d = -1j/(4*mu)*ks*H1ks-1j/(2*omega**2*r)*((kp*H1kp-ks*H1ks)/r+(ks**2*H0ks-kp**2*H0kp)/2)
    phiS2 = 1j/(2*omega**2*r)*(ks*H1ks-kp*H1kp)+1j/(4*omega**2)*(kp**2*H0kp-ks**2*H0ks)
    phiS2d =  1j/(2*omega**2)*(ks**2*H0ks/r-kp**2*H0kp/r-2*(ks*H1ks-kp*H1kp)/(r*r))+1j/(4*omega**2)*(ks**3*H1ks-kp**3*H1kp)
    U1 = (-lam-2*mu)*ny2*b-mu*ny1*a
    U1J = (-lam*ny2*a**2*b-2*mu*ny1*a*b**2+(-lam-2*mu)*ny2*b**3 )/(r**2)
    U2 = (-lam-4*mu)*ny2*b-mu*ny1*a+4*mu*(ny1*a*b**2+ny2*b**3)/(r**2)
    kernel = phiS1d/r*U1+phiS2d/r*U1J+phiS2/(r**2)*U2;
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def D11_0(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = math.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = math.sqrt(a**2 + b**2)
    faktor = 1/(2*sp.pi*(lam+2*mu))
    normAbl = (a*ny1+b*ny2)/r**2
    summand = (mu+2*(lam+mu)*a**2/r**2)*normAbl
    kernel = faktor*summand
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
    
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def D12_0(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = math.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    ty1 = dxds/J
    ty2 = dyds/J
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = math.sqrt(a**2 + b**2)
    faktor = 1/(2*sp.pi*(lam+2*mu))
    normAbl = (a*ny1+b*ny2)/r**2
    tangAbl = (a*ty1+b*ty2)/r**2
    summand = 2*(lam+mu)*a*b*normAbl/r**2 - mu*tangAbl
    kernel = faktor*summand
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
    
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def D21_0(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = math.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    ty1 = dxds/J
    ty2 = dyds/J
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = math.sqrt(a**2 + b**2)
    faktor = 1/(2*sp.pi*(lam+2*mu))
    normAbl = (a*ny1+b*ny2)/r**2
    tangAbl = (a*ty1+b*ty2)/r**2
    summand = 2*(lam+mu)*a*b*normAbl/r**2 + mu*tangAbl
    kernel = faktor*summand
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
    
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def D22_0(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = math.sqrt(dxds**2 + dyds**2)
    
    ny1 = dyds/J
    ny2 = -dxds/J
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = math.sqrt(a**2 + b**2)
    faktor = 1/(2*sp.pi*(lam+2*mu))
    normAbl = (a*ny1+b*ny2)/r**2
    summand = (mu+2*(lam+mu)*b**2/r**2)*normAbl
    kernel = faktor*summand
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
    
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def S11(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    s1 = 1j/(4*mu)*H0ks-1j/(4*omega**2*r)*(ks*H1ks-kp*H1kp)
    s2 = 1j/(4*omega**2)*(2/r*(ks*H1ks-kp*H1kp)-(ks**2*H0ks-kp**2*H0kp))
    kernel = s1 + s2*a**2 /(r**2)
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag


def S12(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    s2 = 1j/(4*omega**2)*(2/r*(ks*H1ks-kp*H1kp)-(ks**2*H0ks-kp**2*H0kp))
    kernel = s2*a*b /(r**2)
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag
    
    
def S22(s,x1,x2,x3,y1,y2,y3,xh,yh,basis,part):
    # Jacobian
    dxds = (4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3
    dyds = (4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3
    J = cmath.sqrt(dxds**2 + dyds**2)
    
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    
    a = xh-(lagrange1(s)*x1 + lagrange2(s)*x2 + lagrange3(s)*x3)
    b = yh-(lagrange1(s)*y1 + lagrange2(s)*y2 + lagrange3(s)*y3)
    r = cmath.sqrt(a**2 + b**2)
    H0ks = hankel1(0,ks*r)
    H1ks = hankel1(1,ks*r)
    H0kp = hankel1(0,kp*r)
    H1kp = hankel1(1,kp*r)
    s1 = 1j/(4*mu)*H0ks-1j/(4*omega**2*r)*(ks*H1ks-kp*H1kp)
    s2 = 1j/(4*omega**2)*(2/r*(ks*H1ks-kp*H1kp)-(ks**2*H0ks-kp**2*H0kp))
    kernel = s1 + s2*b**2 /(r**2)
    
    if basis==1:
        f = lagrange1a(s)
    elif basis==2:
        f = lagrange2a(s)
    elif basis==3:
        f = lagrange3a(s)
        
    if part==0:
        return (kernel*J*f).real
    elif part==1:
        return (kernel*J*f).imag
    

def lagrange1(s):
    u=1-s
    return u*(2*u-1)


def lagrange1a(s):
    u=1-s
    return (u-alpha)/(1-2*alpha) * (1-2*s)/(1-2*alpha)


def lagrange2(s):
    u=1-s
    return 4*s*u


def lagrange2a(s):
    u=1-s
    return 4*(s-alpha)/(1-2*alpha) * (u-alpha)/(1-2*alpha)


def lagrange3(s):
    return s*(2*s-1)


def lagrange3a(s):
    return (s-alpha)/(1-2*alpha) * (2*s-1)/(1-2*alpha)


def index(j,n):
    if j==n:
        return 0
    else:
        return j+1
    
    
def testD(ns,omegas):
    if ns%2:
        print("n has to be even\n")
        return sys.float_info.max
    
    alpha = (1-math.sqrt(3/5))/2
    n = ns
    faces = int(n/2)
    mu = 1
    lam = 1
    omega = omegas
    arad = 2
    brad = 1
    
    # create triangulation mit knodes v_i
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
        ind2 = index(k,n-1)
        ind3 = index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = lagrange1(s)*x1+lagrange2(s)*x2+lagrange3(s)*x3
        vhi[L,1] = lagrange1(s)*y1+lagrange2(s)*y2+lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = lagrange1(s)*x1+lagrange2(s)*x2+lagrange3(s)*x3
        vhi[L+1,1] = lagrange1(s)*y1+lagrange2(s)*y2+lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = lagrange1(s)*x1+lagrange2(s)*x2+lagrange3(s)*x3
        vhi[L+2,1] = lagrange1(s)*y1+lagrange2(s)*y2+lagrange3(s)*y3
        k = k+2
        L = L+3
    
    initTest(mu,lam,omega,alpha,faces,n,vi,vhi)
      
    A = createDmatrix()
    
    rhs = np.zeros(2*3*faces, dtype=np.complex)
    for i in range(0,3*faces):
        rhs[2*i:2*i+2] = truesolution(vhi[i,0], vhi[i,1])
     
    lsg = np.linalg.solve(0.5*np.identity(2*3*faces)+A, rhs)
    
    P = np.array([3,3], dtype=np.float)
    sol = 0
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = index(k,n-1)
        ind3 = index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        k = k+2
        E1 = integrateD(x1, x2, x3, y1, y2, y3, P[0], P[1], 1)
        E2 = integrateD(x1, x2, x3, y1, y2, y3, P[0], P[1], 2)
        E3 = integrateD(x1, x2, x3, y1, y2, y3, P[0], P[1], 3)
        sol = sol + np.matmul(E1,lsg[2*L:2*L+2]) + np.matmul(E2,lsg[2*(L+1):2*(L+1)+2]) + np.matmul(E3,lsg[2*(L+2):2*(L+2)+2])
        L = L+3
        
    truesol = truesolution(P[0], P[1])
    return sp.linalg.norm(abs(sol-truesol), 2)


def testS(ns,omegas):
    if ns%2:
        print("n has to be even\n")
        return sys.float_info.max
    
    alpha = (1-math.sqrt(3/5))/2
    n = ns
    faces = int(n/2)
    mu = 1
    lam = 1
    omega = omegas
    arad = 2
    brad = 1
    
    # create triangulation mit knodes v_i
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
        ind2 = index(k,n-1)
        ind3 = index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        s = alpha
        vhi[L,0] = lagrange1(s)*x1+lagrange2(s)*x2+lagrange3(s)*x3
        vhi[L,1] = lagrange1(s)*y1+lagrange2(s)*y2+lagrange3(s)*y3
        s = 0.5
        vhi[L+1,0] = lagrange1(s)*x1+lagrange2(s)*x2+lagrange3(s)*x3
        vhi[L+1,1] = lagrange1(s)*y1+lagrange2(s)*y2+lagrange3(s)*y3
        s = 1-alpha
        vhi[L+2,0] = lagrange1(s)*x1+lagrange2(s)*x2+lagrange3(s)*x3
        vhi[L+2,1] = lagrange1(s)*y1+lagrange2(s)*y2+lagrange3(s)*y3
        k = k+2
        L = L+3
    
    initTest(mu,lam,omega,alpha,faces,n,vi,vhi)
    
    A = createSmatrix()
    
    rhs = np.zeros(2*3*faces, dtype=np.complex)
    for i in range(0,3*faces):
        rhs[2*i:2*i+2] = truesolution(vhi[i,0], vhi[i,1])
    
     
    lsg = np.linalg.solve(A, rhs)
    
    P = np.array([3,3], dtype=np.float)
    sol = 0
    k = 0
    L = 0
    for j in range(0,faces):
        ind1 = k
        ind2 = index(k,n-1)
        ind3 = index(k+1,n-1)
        x1 = vi[ind1,0]
        x2 = vi[ind2,0]
        x3 = vi[ind3,0]
        y1 = vi[ind1,1]
        y2 = vi[ind2,1]
        y3 = vi[ind3,1]
        k = k+2
        E1 = integrateS(x1, x2, x3, y1, y2, y3, P[0], P[1], 1, 0)
        E2 = integrateS(x1, x2, x3, y1, y2, y3, P[0], P[1], 2, 0)
        E3 = integrateS(x1, x2, x3, y1, y2, y3, P[0], P[1], 3, 0)
        sol = sol + np.matmul(E1,lsg[2*L:2*L+2]) + np.matmul(E2,lsg[2*(L+1):2*(L+1)+2]) + np.matmul(E3,lsg[2*(L+2):2*(L+2)+2])
        L = L+3
        
    truesol = truesolution(P[0], P[1])
    return sp.linalg.norm(abs(sol-truesol), 2)
    
    
def truesolution(x,y):
    r = cmath.sqrt(x**2 + y**2)
    kp = omega/ cmath.sqrt(lam+2*mu)
    ks = omega/ cmath.sqrt(mu)
    s1 = 1j/(4*mu)*hankel1(0,ks*r)-1j/(4*omega**2*r)*(ks*hankel1(1,ks*r)-kp*hankel1(1,kp*r))
    s2 = 1j/(4*omega**2)*(2/r*(ks*hankel1(1,ks*r)-kp*hankel1(1,kp*r))-(ks**2*hankel1(0,ks*r)-kp**2*hankel1(0,kp*r)))
    J = np.array([[x**2, x*y], [x*y, y**2]], dtype=np.complex) / (x**2 + y**2)
    tmp = s1 * np.identity(2,dtype=np.complex) + s2 * J
    return tmp[:,0]

# initialisiere Parameter fuer testD und testS
def initTest(imu,ilam,iomega,ialpha,ifaces,inn,ivi,ivhi):
    global mu
    global lam
    global omega
    global alpha
    global faces
    global n
    global vi
    global vhi
    
    mu = imu
    lam = ilam
    omega = iomega
    alpha = ialpha
    faces = ifaces
    n = inn
    vi = ivi
    vhi = ivhi
    
#initialisiere Parameter fuer beyn
def initBeyn(imu,ilam,isigma1,isigma2,ialpha,ifaces,inn,ivi,ivhi,iN,il,im,itol,iR,imp, ifile):
    global mu
    global lam
    global sigma1
    global sigma2
    global alpha
    global faces
    global n
    global vi
    global vhi
    global N
    global l
    global m
    global tol_rank
    global R
    global mp
    global filename
    
    mu = imu
    lam = ilam
    sigma1 = isigma1
    sigma2 = isigma2
    alpha = ialpha
    faces = ifaces
    n = inn
    vi = ivi
    vhi = ivhi
    N = iN
    l = il
    m = im
    tol_rank = itol
    R = iR
    mp = imp
    filename = ifile
    
    
def printeigs(sb):
    np.set_printoptions(precision=12)
    #with open(filename, "w") as external_file:
    print("mu=" ,mu)
    print("lambda=" ,lam)
    print("alpha=" ,alpha)
    print("faces=" ,faces)
    print("n=" ,n)
    print("sigma1=" ,sigma1)
    print("sigma2=" ,sigma2)
    print("l=" ,l)
    print("m=" ,m)
    print("N=" ,N)
    print("tol_rank=" ,tol_rank)
    print("mp=" ,mp)
    print("R=" ,R)
    print()
    for i in sb:
        print(i)
    #external_file.close()
    
