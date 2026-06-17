import os
import re
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, Data, RealData
from scipy.constants import Boltzmann, Avogadro, c, atomic_mass, hbar, electron_volt, gravitational_constant, elementary_charge, electron_mass, proton_mass, mu_0, epsilon_0, pi, h, Wien
plt.rcParams['figure.dpi'] = 400
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   purple = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'
   
def exp_nobg(p, x):
    """
    Exponential without background:
        y = p[0] * exp(p[1] * x)
    Parameters:
        p[0] = A
        p[1] = b
    """
    return p[0] * np.exp(p[1] * x)

def power(p, x):
    x = np.clip(x, 1e-6, None)  # Avoid zero or negative values
    return p[0] * x**p[1]

def gaussianfunc(p,x):
    return p[0]/(p[2]*np.sqrt(2*np.pi))*np.exp(-(x-p[1])**2/(2*p[2]**2))

def linearfunc(p,x):
    return p[0]*x + p[1]

def gaussianlinear(p,x):
    return gaussianfunc(p[0:3],x) + linearfunc(p[3:5],x)

def expfunc(p, x):
    return (p[0]*np.exp(-x*p[1])) + p[2]

def expfunc_bg(p, x):
    print("Parameters (p):", p)
    print("Length of parameters (p):", len(p))
    R0, lambda_, B = p  # This line will throw an error if p has fewer than 3 elements
    return R0 * np.exp(-lambda_ * x) + B

def residual(p, func, xvar, yvar, err):
    print("Length of parameters (p):", len(p))
    print("Parameters (p):", p)
    return (func(p, xvar) - yvar) / err

def format_with_uncertainty(value, uncertainty):
    """Formats a value with its uncertainty, keeping the correct significant figures and two decimal places."""
    if uncertainty == 0:
        return f"{value:.2f}"  # Always display with two decimal places if uncertainty is zero

    # Determine the number of significant figures from the uncertainty
    sig_figs = -int(np.floor(np.log10(abs(uncertainty)))) + 1
    
    # Round the value and uncertainty to two decimal places
    rounded_value = round(value, 2)
    rounded_uncertainty = round(uncertainty, 2)

    return f"{rounded_value:.2f} \\pm {rounded_uncertainty:.2f}"

def data_fit(p0,func,xvar, yvar, err,tmi=0):
    try:
        fit = optimize.least_squares(residual, p0, args=(func,xvar, yvar, err),verbose=tmi)
    except Exception as error:
        print("Something has gone wrong:",error)
        return p0,np.zeros_like(p0),np.nan,np.nan
    pf = fit['x']

    print()

    try:
        cov = np.linalg.inv(fit['jac'].T.dot(fit['jac']))          
        # This computes a covariance matrix by finding the inverse of the Jacobian times its transpose
        # We need this to find the uncertainty in our fit parameters
    except:
        # If the fit failed, print the reason
        print('Fit did not converge')
        print('Result is likely a local minimum')
        print('Try changing initial values')
        print('Status code:', fit['status'])
        print(fit['message'])
        return pf,np.zeros_like(pf),np.nan,np.nan
            #You'll be able to plot with this, but it will not be a good fit.

    chisq = sum(residual(pf,func,xvar, yvar, err) **2)
    dof = len(xvar) - len(pf)
    red_chisq = chisq/dof
    pferr = np.sqrt(np.diagonal(cov)) # finds the uncertainty in fit parameters by squaring diagonal elements of the covariance matrix
    print('Converged with chi-squared {:.2f}'.format(chisq))
    print('Number of degrees of freedom, dof = {:.2f}'.format(dof))
    print('Reduced chi-squared {:.2f}'.format(red_chisq))
    print()
    Columns = ["Parameter #","Initial guess values:", "Best fit values:", "Uncertainties in the best fit values:"]
    print('{:<11}'.format(Columns[0]),'|','{:<24}'.format(Columns[1]),"|",'{:<24}'.format(Columns[2]),"|",'{:<24}'.format(Columns[3]))
    for num in range(len(pf)):
        print('{:<11}'.format(num),'|','{:<24.3e}'.format(p0[num]),'|','{:<24.3e}'.format(pf[num]),'|','{:<24.3e}'.format(pferr[num]))
    return pf, pferr, chisq,dof