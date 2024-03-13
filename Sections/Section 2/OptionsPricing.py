from scipy.stats import norm
import numpy as np

N = norm.cdf

class OptionsPricing():
    
    def __init__(self, S=100, K=100, T=1, r=0.05, q=0, sigma=0.2):
        
        self.S = S
        self.K = K
        self.T = T
        self.r = r
        self.q = q
        self.sigma = sigma
        
    def PRICE_CALL(self):

        d1 = (np.log(self.S / self.K) + (self.r - self.q + self.sigma**2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
        d2 = d1 - self.sigma * np.sqrt(self.T)
        
        return self.S * np.exp(-self.q * self.T) * N(d1) - self.K * np.exp(-self.r * self.T) * N(d2)
    
    def PRICE_PUT(self):
        
        d1 = (np.log(self.S / self.K) + (self.r - self.q + self.sigma**2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
        d2 = d1 - self.sigma * np.sqrt(self.T)
        
        return self.K * np.exp(-self.r * self.T) * N(-d2) - self.S * np.exp(-self.q * self.T) * N(-d1)

    def DELTA_CALL(self):
        
        d1 = (np.log(self.S / self.K) + (self.r - self.q + self.sigma**2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
        delta = np.exp(-self.q * self.T) * N(d1)
        
        return delta
    
    def DELTA_PUT(self):
        
        d1 = (np.log(self.S / self.K) + (self.r - self.q + self.sigma**2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
        delta = np.exp(-self.q * self.T) * (N(d1) - 1)
        
        return delta
    
    def GAMMA_CALL(self):
        
        d1 = (np.log(self.S / self.K) + (self.r - self.q + self.sigma**2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
        gamma = np.exp(-(d1**2 / 2 + self.q * self.T)) / (self.S * self.sigma * np.sqrt(2 * np.pi * self.T))
        
        return gamma
    
    def GAMMA_CALL(self):
        
        d1 = (np.log(self.S / self.K) + (self.r - self.q + self.sigma**2 / 2) * self.T) / (self.sigma * np.sqrt(self.T))
        gamma = np.exp(-(d1**2 / 2 + self.q * self.T)) / (self.S * self.sigma * np.sqrt(2 * np.pi * self.T))
        
        return gamma

