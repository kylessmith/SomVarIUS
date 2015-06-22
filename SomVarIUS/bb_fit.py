from scipy.special import gammaln, beta
from scipy.misc import factorial, comb
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt


class bb_model:
    
    def __init__(self, xi=None, ni=None, alpha=None, beta=None):
        
        if alpha != None and beta != None:
            self.alpha = float(alpha)
            self.beta = float(beta)
            
        elif alpha != None and xi != None:
            print 'Error: class can only accept xi and ni or alpha and beta'
            exit()
        
        else:
        
            self.xi = np.array(xi, dtype=float)
            self.ni = np.array(ni, dtype=float)
            self.pi=self.xi/self.ni

            # initializing weights and gamma. The weights will be updated recursively later on.
            wi = np.repeat(1, len(self.xi))
            self.gamma=1

            for i in xrange(1,1000):
                w = sum(wi)
                p = sum(wi*self.pi)/w
                q = (1-p)
    
                S = sum(wi*(self.pi-p)**2)
    
                mu = p
                self.old_gamma = self.gamma
                self.gamma = (S-p*q*(sum((wi/self.ni)*(1-wi/w))))/(p*q*(sum(wi*(1-wi/w))-sum((wi/self.ni)*(1-wi/w))))
                wi=self.ni/(1+self.gamma*(self.ni-1))
                if(abs(self.old_gamma-self.gamma)<0.001):
                    break
        
            self.rho=self.theta=self.gamma/(1-self.gamma)
            self.beta=(1-mu)*(1/self.gamma-1)
            self.alpha=mu*self.beta/(1-mu)

        self.mean = self.alpha/(self.alpha+self.beta)
        var_num = (self.alpha*self.beta)*(self.alpha+self.beta+1)
        var_den = ((self.alpha+self.beta)**2)*(self.alpha+self.beta+1)
        self.var = var_num/var_den
        skew_1 = (self.alpha+self.beta+2)*(self.beta-self.alpha)
        skew_2 = (self.alpha+self.beta+2)
        skew_3 = sqrt((1+self.alpha+self.beta)/((self.alpha*self.beta)*(1+self.alpha+self.beta)))
        self.skew = (skew_1/skew_2)*skew_3
        
    def update_mean(self, n):
        self.mean = self.mean*n
    
    def update_var(self, n):
        var_num = (self.alpha*self.beta*n)*(self.alpha+self.beta+n)
        var_den = ((self.alpha+self.beta)**2)*(self.alpha+self.beta+1)
        self.var = var_num/var_den
        
    def update_skew(self, n):
        skew_1 = (self.alpha+self.beta+2*n)*(self.beta-self.alpha)
        skew_2 = (self.alpha+self.beta+2)
        skew_3 = sqrt((1+self.alpha+self.beta)/((n*self.alpha*self.beta)*(n+self.alpha+self.beta)))
        self.skew = (skew_1/skew_2)*skew_3
        
    def pdf(self, k, n):
        p = comb(n,k) * beta(k+self.alpha, n-k+self.beta) / beta(self.alpha,self.beta)
        return p
    
    def cdf(self, k, n):
        
        def log_dcm(counts,alpha):
            N = sum(counts)
            A = sum(alpha)
            return gammaln(N+1) - sum(gammaln(counts+1)) + gammaln(A) - gammaln(N + A) + sum(gammaln(counts+alpha) - gammaln(alpha))

        def log_betabin(a,b,k,N):
            return log_dcm(np.array([k,N-k]),np.array([a,b]))
        
        if k > 0.5 * n:
            p = 1. - sum([np.exp(log_betabin(self.alpha, self.beta, x, n)) for x in range(k+1,n)])
        else:
            p = sum([np.exp(log_betabin(self.alpha, self.beta, x, n)) for x in range(k+1)])
            
        return p
        
    def plot_pdf(self, n,out_name=None):
        x = np.arange(n)
        funct = np.vectorize(self.pdf)
        y = funct(x, x[-1])
        line = plt.plot(x,y)
        if out_name == None:
            plt.show()
        else:
            plt.savefig(out_name)
            
    def plot_cdf(self, n,out_name=None):
        x = np.arange(n)
        funct = np.vectorize(self.cdf)
        y = funct(x, x[-1])
        line = plt.plot(x,y)
        if out_name == None:
            plt.show()
        else:
            plt.savefig(out_name)