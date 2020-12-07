import numpy as np, random
from scipy import optimize, integrate, interpolate
from math import *


########## Square ##########
class Square:
    def __init__(self, build=False):
        self.id = "S"
        self.numJ = 2
    def build(self):
        num=100000
        table_k=np.linspace(-80, 80, num)
        table_Tc=np.zeros(num)
        for i in range(num):
            table_Tc[i]=self.Tc_binary(table_k[i])
            if i%(num//10)==num//10-1:
                print("*\n" if i==num-1 else ".", end='')
        self.k_to_Tc = interpolate.interp1d(table_k, table_Tc, fill_value=(table_Tc[0], table_Tc[-1]), bounds_error=False)

    def J0(self,Js):
        return 2*np.sum(Js)
    def gen_Js(self, k):
        Js = np.array([1.0, 10**k])
        Js/=self.J0(Js)
        return Js

    def z(self, Js, T):
        assert(len(Js)==self.numJ)
        return 2*Js/T
    def Tc_binary(self, k):
        Js=self.gen_Js(k)
        f=lambda T: np.prod(np.sinh(self.z(Js,T)))-1
        return optimize.brentq(f, 5e-3, 0.57)
    def Tc(self, k):
        Js=self.gen_Js(k)
        f=lambda T: np.prod(np.sinh(self.z(Js,T)))-1
        return optimize.newton(f, self.k_to_Tc(k))
    def Tc_Js(self, Js):
        f=lambda T: np.prod(np.sinh(self.z(Js,T)))-1
        return optimize.brentq(f, 5e-3, 0.57)

    def Lmd(self, Js, q1, q2):
        return 1.0/np.sum(Js) * ( Js[0]*cos(q1) + Js[1]*cos(q2) )
    def P(self, X, Js):
        return integrate.dblquad(lambda q1, q2: self.integrand(X,Js,q1,q2), -pi, pi, lambda q2: -pi, lambda q2: pi)
    def integrand(self, X, Js, q1, q2):
        return 1/(2*pi)**2/ ( 1-X*self.Lmd(Js, q1, q2) )
    def PolyCoeff(self, n, Js):
        return integrate.dblquad(lambda q1, q2: self.Lmd(Js,q1,q2)**n/(2*pi)**2 , -pi, pi, lambda q2: -pi, lambda q2: pi)
    

########## Triangle ##########
class Triangle:
    def __init__(self):
        self.id = "T"
        self.numJ = 3
    
    def J0(self,Js):
        return 2*np.sum(Js)
    def gen_Js(self,k1, k2):
        Js = np.array([1.0, 10**k1, 10**k2])
        Js/=self.J0(Js)
        return Js
    
    def z(self, Js, T):
        assert(len(Js)==self.numJ)
        return np.exp(-2.0*Js/T)
    def Tc(self, k1, k2):
        Js=self.gen_Js(k1, k2)
        f = lambda T: np.dot(self.z(Js, T), np.roll(self.z(Js, T), 1))-1 
        return optimize.brentq(f, 1e-3, 0.61)
    def Tc_Js(self,Js):
        f = lambda T: np.dot(self.z(Js, T), np.roll(self.z(Js, T), 1))-1 
        return optimize.brentq(f, 1e-3, 0.61)
    
    def Lmd(self, Js, q1, q2):
        return 1.0/np.sum(Js) * ( Js[0]*cos(q1) + Js[1]*cos(q2) + Js[2]* cos(q1+q2) ) 
    
    def P(self, X, Js):
        assert(len(Js)==self.numJ)
        return integrate.dblquad(lambda q1, q2: self.integrand(X,Js,q1,q2), -pi, pi, lambda q2: -pi, lambda q2: pi)
    def integrand(self, X, Js, q1, q2):
        return 1/(2*pi)**2/ ( 1-X*self.Lmd(Js, q1, q2) )
    
    def PolyCoeff(self, n, Js):
        assert(len(Js)==self.numJ)
        return integrate.dblquad(lambda q1, q2: self.Lmd(Js,q1,q2)**n/(2*pi)**2 , -pi, pi, lambda q2: -pi, lambda q2: pi)


########## HoneyComb ##########
class HoneyComb:
    def __init__(self):
        self.id = "H"
        self.numJ = 3
    
    def J0(self,Js):
        return np.sum(Js)
    def gen_Js(self,k1, k2):
        Js = np.array([1.0, 10**k1, 10**k2])
        Js/=self.J0(Js)
        return Js
    
    def z(self, Js, T):
        return np.tanh(Js/T)
    def Tc(self, k1, k2):
        Js=self.gen_Js(k1, k2)
        f = lambda T: np.dot(self.z(Js, T), np.roll(self.z(Js, T), 1))-1 
        return optimize.brentq(f, 1e-10, 0.51)
    def Tc_Js(self, Js):
        f = lambda T: np.dot(self.z(Js, T), np.roll(self.z(Js, T), 1))-1 
        return optimize.brentq(f, 1e-10, 0.51)
    
    def Lmd(self, Js, q1, q2):
        return 1.0/(np.sum(Js*Js)+2*np.sum(Js*np.roll(Js, 1))) * ( np.sum(Js*Js) + 2*Js[0]*Js[1]*cos(q1)+ 2*Js[1]*Js[2]*cos(q2) +  2*Js[2]*Js[0]*cos(q1+q2) )
        
    def P(self, X, Js):
        assert(len(Js)==self.numJ)
        return integrate.dblquad(lambda q1, q2: self.integrand(X,Js,q1,q2), -pi, pi, lambda q2: -pi, lambda q2: pi)
    def integrand(self, X, Js, q1, q2):
        return 1/(2*pi)**2/ ( 1 -  X**2*self.Lmd(Js, q1, q2) )
    
    def PolyCoeff(self, n, Js):
        assert(len(Js)==self.numJ)
        if n%2==1:
          return [0.0, 0.0]
        else:
          n//=2
          return integrate.dblquad(lambda q1, q2: self.Lmd(Js,q1,q2)**n/(2*pi)**2 , -pi, pi, lambda q2: -pi, lambda q2: pi)





########## Cubic ##########
class Cubic:
  sample=10
  ratio=np.array([1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001])
  Tc_coeff=np.array([4.5117, 2.9297, 1.8139, 1.34307, 1.03984, 0.78297, 0.65251, 0.55589, 0.46200, 0.40832])
  
  def __init__(self):
    self.id = "C"
    self.numJ = 3
    
  def J0(self,Js):
    return 2*np.sum(Js)
    
  def gen_Js(self, i):
    assert(i<Cubic.sample)
    Js = np.array([1, Cubic.ratio[i], Cubic.ratio[i]])
    Js/=self.J0(Js)
    return Js
    
  def Tc(self, Js, i):
    return Js[0]*Cubic.Tc_coeff[i]
    
  def Lmd(self, Js, q1, q2, q3):
    return 1.0/np.sum(Js) * ( Js[0]*cos(q1) + Js[1]*cos(q2) + Js[2]*cos(q3) )

  def P(self, X, Js):
    assert(len(Js)==self.numJ)
    return integrate.tplquad(lambda q1, q2, q3: self.integrand(X,Js,q1,q2,q3), -pi, pi, lambda q3: -pi, lambda q3: pi, lambda q3,q2: -pi, lambda q3,q2: pi)
  def integrand(self, X, Js, q1, q2, q3):
    return 1/(2*pi)**3/ ( 1-X*self.Lmd(Js, q1, q2, q3) )

  def PolyCoeff(self, n, Js):
    assert(len(Js)==self.numJ)
    return integrate.tplquad(lambda q1, q2, q3: self.Lmd(Js,q1,q2,q3)**n/(2*pi)**3 , -pi, pi, lambda q3: -pi, lambda q3: pi, lambda q3,q2: -pi, lambda q3,q2: pi)
    
    
    
########## Chain ##########
class Chain:
  
  def __init__(self):
    self.id = "Ch"
    self.numJ = 1
    
  def gen_Js(self):
    return 0.5
    
  def Lmd(self, Js, q1):
    return cos(q1)

  def P(self, X, Js):
    assert(len(Js)==self.numJ)
    return integrate.quad(lambda q1: self.integrand(X,Js,q1), -pi, pi)
  def integrand(self, X, Js, q1):
    return 1/(2*pi)/ ( 1-X*self.Lmd(Js, q1) )

  def PolyCoeff(self, n, Js):
    assert(len(Js)==self.numJ)
    return integrate.dblquad(lambda q1: self.Lmd(Js,q1)**n/(2*pi) , -pi, pi)
    
H=HoneyComb()
T=Triangle()
S=Square()
C=Cubic()
Ch=Chain()