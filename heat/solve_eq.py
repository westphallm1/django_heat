#for use in eval'ing
from numpy import *
import numpy as np



def valfstr(eqnstr,x):
    """Parses a string in the form of f(x) into a value at the specified x """
    #exit if any functions from "suspect" modules are being called
    if 'os' in eqnstr or 'sys' in eqnstr:
        return None
        
    return eval(eqnstr)

def valsfstr(eqnstr, xmin, xmax, num = 50):
    """Parses a string in the form of f(x) into a numpy array over the 
    specified range
    """
    #exit if any functions from "suspect" modules are being called
    if 'os' in eqnstr or 'sys' in eqnstr:
        return None
        
    x = np.linspace(xmin,xmax,num)
    f_x = eval(eqnstr)
    #turn constants into arrays
    if not isinstance(f_x, np.ndarray):
        f_x = np.ones(num)*f_x
    return x,f_x
    
    
def find_0(eqnstr,xmin,xmax, num_places = 4):
    """approximate the zero of a function"""
    closeness = 10.**-num_places
    old_0 = 0.
    new_0 = 2.*closeness
    num_iters = 0
    while np.abs(old_0 - new_0) > closeness:
        num_iters +=1
        #give up if it's not converging
        if num_iters > 1000: return None
            
        old_0 = new_0
        
        x,f_x = valsfstr(eqnstr, xmin, xmax)
        
        #find where f_x crosses from positive to negative
        for i in range(f_x.size):
            try:
                if np.sign(f_x[i]) != np.sign(f_x[i+1]):
                    break
            except IndexError:
                #got to end of xrange without finding a zero
                return None
        
        xmin = x[i]
        xmax = x[i+1]
        new_0 = xmin
        
    return np.round(new_0,num_places) 
        
        