#!/usr/bin/env python
# coding: utf-8
#python version 3.8.16
#bilby version 1.1.1
#numpy version 1.19.5
#astropy version 5.0.0

import numpy as np
from scipy.special import erf
from scipy.stats import beta as beta_dist
from scipy.stats import truncnorm
from scipy.interpolate import interp1d
import bilby
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
from bilby.core.prior import Uniform
from bilby.core.sampler import run_sampler
from bilby.core.prior import LogUniform
from bilby.hyper.likelihood import HyperparameterLikelihood
from scipy import interpolate
from scipy import integrate
import random
import argparse#
import os
import glob
import warnings
warnings.filterwarnings("ignore")

try:
    import cupy as xp
    from cupyx.scipy.special import erf, gammaln  # noqa

    CUPY_LOADED = True
except ImportError:
    import numpy as xp
    from scipy.special import erf, gammaln  # noqa

    CUPY_LOADED = False


#summary:
#1#####fundmental models
###the structure of fundmental models
#+uniform model
#+log-uniform
#+pow
#+log-normal
#+SST
#+Gamma
############
#2####Gaussian distribution series models
###the structure of Gaussian distribution series models
#+Gaussian model 
#+Gaussian model with fixed min
#+Gaussian model with fixed max
############
#3#### 2-Gaussian distributions series
###the structure of 2-Gaussian distribution series models
#+2G with free paramaters maxmum and minimum mass 
#+2G fixed max
#+2G fixed min
#+2G fixed min and max
############
#4####3Gaussian distribution series
###the structure of 3-Gaussian distribution series models
#+3G fixed min and max
############
#5######turn on power-law model series
###the structure of turn on power-law model series
#+top
#+top with fixed m_max
#+top + G 
#+top + G with fixed m_max
###############

#begin uniform model
def hyper_prior_U(dataset,mlo,mup):
    return (( dataset['mu'] >= mlo) & (dataset['mu'] <= mup) ) / abs(mup-mlo)

hp_priors_U = dict(mlo=Uniform(0.9, 1.5, 'mlo',r'$\rm{m^l}$'),
                mup=Uniform(1.5, 2.9, 'mup',r'$\rm{m^u}$'))
#end uniform model

#begin log-uniform
def hyper_prior_logu(dataset,mlo,mup):
    return (( dataset['mu'] >= mlo) & (dataset['mu'] <= mup))/( dataset['mu'] * np.log(mup / mlo) )

hp_priors_logu = dict(mlo=Uniform(0.9, 1.5, 'mlo',r'$\rm{m^l}$'),
                mup=Uniform(1.5, 2.9, 'mup',r'$\rm{m^u}$'))
#end log-uniform

#begin pow
def hyper_prior_pow(dataset,mlo,mup,beta):
    beta=-1*beta
    return (( dataset['mu'] >= mlo) & (dataset['mu'] <= mup))*((1+beta)/(mup**(1+beta)-mlo**(1+beta)))*dataset['mu']**beta

hp_priors_pow = dict(mlo=Uniform(0.9, 1.5, 'mlo',r'$\rm{m^l}$'),
                mup=Uniform(1.5, 2.9, 'mup',r'$\rm{m^u}$'),
                beta=Uniform(-5, 25, 'beta','$\\beta$'))
#end pow

#begin log-normal
def hyper_prior_lognorm(dataset, s_mu, s_sigma):
    return np.exp(- (np.log(dataset['mu']) - s_mu)**2 / (2 * s_sigma**2)) /\
        (2 * np.pi * s_sigma**2)**0.5/(dataset['mu'])
hp_priors_lognorm = dict(s_mu=Uniform(0.01, 1, 's_mu', '$\mu$'),
                 s_sigma=Uniform(0.01, 0.5, 's_sigma', '$\sigma$') )
#end log-normal

#begin SST
from scipy.special import beta
def hyper_prior_sst(dataset, mu,sigma,nu,tau):
        c = 2 * nu * ((1 + nu ** 2) *
                                beta(0.5, tau / 2) *
                                tau ** 0.5) ** -1
        m = ((2 * tau ** 0.5) * (nu - nu ** -1)) / (
                (tau - 1) * beta(0.5, 0.5 * tau))
        s2 = ((tau / (tau - 2)) * (
                nu ** 2 + nu ** -2 - 1) - m ** 2)
        mu_0 = mu - (sigma * m / np.sqrt(s2))
        sigma_0 = sigma / np.sqrt(s2)
        z = (dataset['mu'] - mu_0) / sigma_0
        p = np.where(dataset['mu'] < mu_0,
                     (c / sigma_0) * (1 + ((nu ** 2) * (z ** 2)) / tau) ** (
                             -(tau + 1) / 2),
                     (c / sigma_0) * (1 + (z ** 2) / ((nu ** 2) * tau)) ** (
                             -(tau + 1) / 2))
        return p

hp_priors_sst = dict(mu=Uniform(0.9, 2.9, 'mlo',r'$\rm{m^l}$'),
                sigma=Uniform(0.01, 2, 'sigma',r'$\rm{m^u}$'),
                nu=Uniform(0,8,'nu'),
                   tau=Uniform(2.001,20,'tau') )
#end SST

#begin gamma distribution
from scipy.special import beta
from scipy.special import gamma
def hyper_prior_gamma(dataset, k,theta):
    return (1 / (gamma(k)*theta**k)) * dataset['mu']**(k-1) *np.exp(-dataset['mu']/theta)

hp_priors_gamma = dict(k=Uniform(0, 80, 'k',r'$k$'),
                theta=Uniform(0.01, 0.1, 'theta',r'$\theta$') )
#end gamma distribution

#####Gaussian distribution series models
###the structure of Gaussian distribution series models
#+Gaussian model 
#+Gaussian model with fixed maxmum and minimum mass
#+Gaussian model with fixed min
#+Gaussian model with fixed max
############

#begin Gaussian model 
def hyper_prior_G(dataset, mu, sigma,mlo,mup):
    normalisingTerm = 0.5 * ( erf((mu-mlo)/(np.sqrt(2) * sigma)) -  erf((mu-mup)/(np.sqrt(2) * sigma)) )
    return ( ( dataset['mu'] >= mlo) & (dataset['mu'] <= mup))*((mu>mlo)&(mu<mup)) * (np.exp(- (dataset['mu'] - mu)**2 / (2 * sigma**2)) /\
        (2 * np.pi * sigma**2)**0.5) / normalisingTerm 
hp_priors_G = dict(mu=Uniform(0.9, 2.9, 's_mu', '$\mu$'),
                 sigma=Uniform(0.01, 2, 's_sigma', '$\sigma$'),
                 mlo=Uniform(0.9, 1.5, 'mlow', '$mlow$'),
                 mup=Uniform(1.5, 2.9, 'mup', '$mup$')
                              )
#end Gaussian model 

#begin Gaussian model with fixed max
def hyper_prior_G_fixed_max(dataset, mu, sigma,mlo):
    mup=2.9
    normalisingTerm = 0.5 * ( erf((mu-mlo)/(np.sqrt(2) * sigma)) -  erf((mu-mup)/(np.sqrt(2) * sigma)) )
    return ( ( dataset['mu'] >= mlo) & (dataset['mu'] <= mup))*((mu>mlo)&(mu<mup)) * (np.exp(- (dataset['mu'] - mu)**2 / (2 * sigma**2)) /\
        (2 * np.pi * sigma**2)**0.5) / normalisingTerm 
hp_priors_G_fixed_max= dict(mu=Uniform(0.9, 2.9, 's_mu', '$\mu$'),
                 sigma=Uniform(0.01, 2, 's_sigma', '$\sigma$'),
                      mlo=Uniform(0.9, 1.5, 'mlow', '$mlow$')     )
#end Gaussian with fixed upper and lower mass

##### 2-Gaussian distributions series
###the structure of 2-Gaussian distribution series models
#+2G with free paramaters maxmum and minimum mass 
#+2G fixed max
#+2G fixed min
#+2G fixed min and max
############

#2G with free paramaters maxmum and minimum mass 
def hyper_prior_2G(dataset, mu1, sigma1,mu2,sigma2,alpha,mup,mlo):
    normalisingTerm1 = 0.5 * ( erf((mu1-mlo)/(np.sqrt(2) * sigma1)) -  erf((mu1-mup)/(np.sqrt(2) * sigma1)) )
    normalisingTerm2 = 0.5 * ( erf((mu2-mlo)/(np.sqrt(2) * sigma2)) -  erf((mu2-mup)/(np.sqrt(2) * sigma2)) )
    return ((mu2 < mup ) & (mu1 > mlo) & (mu1 < mu2)  & ( dataset['mu'] >= mlo) & (dataset['mu'] <= mup)) *\
        ( (( alpha*(np.exp(- (dataset['mu'] - mu1)**2 / (2 * sigma1**2)) /(2 * np.pi * sigma1**2)**0.5)) /normalisingTerm1) +\
        (1-alpha)*( ((np.exp(- (dataset['mu'] - mu2)**2 / (2 * sigma2**2)) /(2 * np.pi * sigma2**2)**0.5) ) / normalisingTerm2) )
hp_priors_2G = dict(mu1=Uniform(0.9, 2.9, 'mu1', '$\mu_1$'),
                 sigma1=Uniform(0.01, 2, 'sigma1', '$\sigma_1$'),
                mu2=Uniform(0.9, 2.9, 'mu2', '$\mu_2$'),
                sigma2=Uniform(0.01, 2, 'sigma2', '$\sigma_2$'),
                alpha=Uniform(0.01, 1, 'alpha', '$\\alpha$'),
                mup=Uniform(1.5, 2.9, 'mup',r'$\rm{m^u}$'),
                mlo=Uniform(0.9, 1.5, 'mlo',r'$\rm{m^l}$') )
#end 2G with free paramaters maxmum and minimum mass 

#2G fixed max
def hyper_prior_2G_fixed_max(dataset, mu1, sigma1,mu2,sigma2,alpha,mlo):
    mup=2.9
    normalisingTerm1 = 0.5 * ( erf((mu1-mlo)/(np.sqrt(2) * sigma1)) -  erf((mu1-mup)/(np.sqrt(2) * sigma1)) )
    normalisingTerm2 = 0.5 * ( erf((mu2-mlo)/(np.sqrt(2) * sigma2)) -  erf((mu2-mup)/(np.sqrt(2) * sigma2)) )
    return ( (mu2 < mup ) & (mu1 > mlo) & (mu1 < mu2)  & ( dataset['mu'] >= mlo) & (dataset['mu'] <= mup)) *\
        ( (( alpha*(np.exp(- (dataset['mu'] - mu1)**2 / (2 * sigma1**2)) /(2 * np.pi * sigma1**2)**0.5)) /normalisingTerm1) +\
        (1-alpha)*( ((np.exp(- (dataset['mu'] - mu2)**2 / (2 * sigma2**2)) /(2 * np.pi * sigma2**2)**0.5) ) / normalisingTerm2) )
hp_priors_2G_fixed_max = dict(mu1=Uniform(0.9, 2.9, 'mu1', '$\mu_1$'),
                 sigma1=Uniform(0.01, 2, 'sigma1', '$\sigma_1$'),
                mu2=Uniform(0.9, 2.9, 'mu2', '$\mu_2$'),
                sigma2=Uniform(0.01, 2, 'sigma2', '$\sigma_2$'),
                alpha=Uniform(0.01, 1, 'alpha', '$\\alpha$'),
                mlo=Uniform(0.9, 1.5, 'mlo',r'$\rm{m^l}$') )
#end 2G fixed max

#begin 2G fixed min
def hyper_prior_2G_fixed_min(dataset, mu1, sigma1,mu2,sigma2,alpha,mup):
    mlo=0.9
    normalisingTerm1 = 0.5 * ( erf((mu1-mlo)/(np.sqrt(2) * sigma1)) -  erf((mu1-mup)/(np.sqrt(2) * sigma1)) )
    normalisingTerm2 = 0.5 * ( erf((mu2-mlo)/(np.sqrt(2) * sigma2)) -  erf((mu2-mup)/(np.sqrt(2) * sigma2)) )
    return ((mu2 < mup ) & (mu1 > mlo) & (mu1 < mu2)  & ( dataset['mu'] >= mlo) & (dataset['mu'] <= mup)) *\
        ( (( alpha*(np.exp(- (dataset['mu'] - mu1)**2 / (2 * sigma1**2)) /(2 * np.pi * sigma1**2)**0.5)) /normalisingTerm1) +\
        (1-alpha)*( ((np.exp(- (dataset['mu'] - mu2)**2 / (2 * sigma2**2)) /(2 * np.pi * sigma2**2)**0.5) ) / normalisingTerm2) )
hp_priors_2G_fixed_min = dict(mu1=Uniform(0.9, 2.9, 'mu1', '$\mu_1$'),
                 sigma1=Uniform(0.01, 2, 'sigma1', '$\sigma_1$'),
                mu2=Uniform(0.9, 2.9, 'mu2', '$\mu_2$'),
                sigma2=Uniform(0.01, 2, 'sigma2', '$\sigma_2$'),
                alpha=Uniform(0.01, 1, 'alpha', '$\\alpha$'),
                mup=Uniform(1.5, 2.9, 'mup',r'$\rm{m^u}$') )
#end 2G fixed min


#begin two-Gausssian model with fixed max and min mass
def hyper_prior_2G_fixed_max_min(dataset, mu1, sigma1,mu2,sigma2,alpha):
    mup=2.9
    mlo=0.9
    normalisingTerm1 = 0.5 * ( erf((mu1-mlo)/(np.sqrt(2) * sigma1)) -  erf((mu1-mup)/(np.sqrt(2) * sigma1)) )
    normalisingTerm2 = 0.5 * ( erf((mu2-mlo)/(np.sqrt(2) * sigma2)) -  erf((mu2-mup)/(np.sqrt(2) * sigma2)) )
    return ((mu2 < mup ) & (mu1 > mlo) & (mu1 < mu2)  & ( dataset['mu'] >= mlo) & (dataset['mu'] <= mup)) *\
        ( (( alpha*(np.exp(- (dataset['mu'] - mu1)**2 / (2 * sigma1**2)) /(2 * np.pi * sigma1**2)**0.5)) /normalisingTerm1) +\
        (1-alpha)*( ((np.exp(- (dataset['mu'] - mu2)**2 / (2 * sigma2**2)) /(2 * np.pi * sigma2**2)**0.5) ) / normalisingTerm2) )
hp_priors_2G_fixed_max_min = dict(mu1=Uniform(0.9, 2.9, 'mu1', '$\mu_1$'),
                 sigma1=Uniform(0.01, 2, 'sigma1', '$\sigma_1$'),
                mu2=Uniform(0.9, 2.9, 'mu2', '$\mu_2$'),
                sigma2=Uniform(0.01, 2, 'sigma2', '$\sigma_2$'),
                alpha=Uniform(0.01, 1, 'alpha', '$\\alpha$'))
#end two-Gausssian model with fixed lower and upper mass

#####3Gaussian distribution series
###the structure of 3-Gaussian distribution series models
#+3G fixed min and max
############

#begin three Gausssian model with fixed lower and upper mass
def hyper_prior_3G_fixed_max_min(dataset, mu1, sigma1,mu2,sigma2,alpha,mu3,sigma3,beta):
    mup=2.9
    mlo=0.9
    normalisingTerm1 = 0.5 * ( erf((mu1-mlo)/(np.sqrt(2) * sigma1)) -  erf((mu1-mup)/(np.sqrt(2) * sigma1)) )
    normalisingTerm2 = 0.5 * ( erf((mu2-mlo)/(np.sqrt(2) * sigma2)) -  erf((mu2-mup)/(np.sqrt(2) * sigma2)) )
    normalisingTerm3 = 0.5 * ( erf((mu3-mlo)/(np.sqrt(2) * sigma3)) -  erf((mu3-mup)/(np.sqrt(2) * sigma3)) )
    if mu1 < mu2 and mu3>mu2  and alpha+beta<=1:
        return ((alpha*(np.exp(- (dataset['mu'] - mu1)**2 / (2 * sigma1**2)) /(2 * np.pi * sigma1**2)**0.5))/normalisingTerm1)\
        +((beta*(np.exp(- (dataset['mu'] - mu2)**2 / (2 * sigma2**2)) /(2 * np.pi * sigma2**2)**0.5))/normalisingTerm2)\
        +(((1-alpha-beta)*(np.exp(- (dataset['mu'] - mu3)**2 / (2 * sigma3**2)) /(2 * np.pi * sigma3**2)**0.5))/normalisingTerm3)
    else:
        return 0
hp_priors_3G_fixed_max_min = dict(mu1=Uniform(0.9, 2.9, 'mu1', '$\mu_1$'),
                 sigma1=Uniform(0.01, 2, 'sigma1', '$\sigma_1$'),
                mu2=Uniform(0.9, 2.9, 'mu2', '$\mu_2$'),
                sigma2=Uniform(0.01, 2, 'sigma2', '$\sigma_2$'),
                alpha=Uniform(0.01, 1, 'alpha', '$\\alpha$'),
                mu3=Uniform(0.9, 2.9, 'mu3', '$\mu_3$'),
                sigma3=Uniform(0.01, 2, 'sigma3', '$\sigma_3$'),
                beta=Uniform(0.01, 1, 'beta', '$\\beta$'))
#end three Gausssian model with fixed lower and upper mass
                

#######turn on power-law model series
###the structure of turn on power-law model series
#+top
#+top with fixed m_max
#+top + G 
#+top + G with fixed m_max
###############

#begin turn_on_pow 
def window(masses, mmin, mmax, delta_m):

    """
    Apply a one sided window between mmin and mmin + delta_m to the
    mass pdf.
    The upper cut off is a step function,
    the lower cutoff is a logistic rise over delta_m solar masses.
    See T&T18 Eqs 7-8
    Note that there is a sign error in that paper.
    S = (f(m - mmin, delta_m) + 1)^{-1}
    f(m') = delta_m / m' + delta_m / (m' - delta_m)
    See also, https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
    """
    window = xp.ones_like(masses)
    if delta_m > 0.0:
        smoothing_region = (masses >= mmin) & (masses < (mmin + delta_m))
        shifted_mass = masses[smoothing_region] - mmin
        if shifted_mass.size:
            exponent = xp.nan_to_num(
                delta_m / shifted_mass + delta_m / (shifted_mass - delta_m)
            )
            window[smoothing_region] = 1 / (xp.exp(exponent) + 1)
    window[(masses < mmin) | (masses > mmax)] = 0
    return window    

def extract_mass_parameters(parameters):
    """extract the parameters of the mass distribution hyperparameters used in
    T&T18 from either a list or dictionary."""
    if isinstance(parameters, list):
        return parameters
    elif isinstance(parameters, dict):
        keys = ['alpha', 'mmin', 'mmax', 'delta_m']
        return [parameters[key] for key in keys]

def ppow(masses, parameters):
    """1d unnormalised powerlaw mass probability with smoothed low-mass end"""
    alpha, mmin, mmax, delta_m = extract_mass_parameters(parameters)
    return masses**(-alpha) * window(masses, mmin, mmax, delta_m) 

def norm_ppow(parameters):
    """normalise ppow, requires m1s, an array of m values, and dm, the spacing of
    that array"""
    m1s = np.linspace(0.9, 2.9, 500)
    return np.trapz(ppow(m1s, parameters), m1s)

def turn_on_pow(masses, parameters, pow_norm):
    alpha, mmin, mmax, delta_m = extract_mass_parameters(parameters)
    p_pow = ppow(masses, parameters) / pow_norm
    return  p_pow 

def hyper_prior_turn_on_pow(dataset, alpha, mmin, mmax, delta_m):
    parameters = dict(
        alpha=alpha, mmin=mmin, mmax=mmax, delta_m=delta_m)
    pow_norm = norm_ppow(parameters)
    probability = turn_on_pow(dataset['mu'], parameters, pow_norm)
    return probability

hp_priors_turn_on_pow= dict(alpha=Uniform(-5, 25, 'alpha', '$\\alpha$'),
                 mmin=Uniform(0.9, 1.5, 'mmin', '$mmin$'),
                mmax=Uniform(1.5, 2.9, 'mmax', '$mmax$'),
                delta_m=Uniform(0.01, 1, 'delta', '$\\delta$'))
#end turn_on_pow 

#begin turn_on_pow with fixed m_max
def window_fix(masses, mmin,  delta_m):
    mmaxs_fix=2.9
    mmax=mmaxs_fix

    """
    Apply a one sided window between mmin and mmin + delta_m to the
    mass pdf.
    The upper cut off is a step function,
    the lower cutoff is a logistic rise over delta_m solar masses.
    See T&T18 Eqs 7-8
    Note that there is a sign error in that paper.
    S = (f(m - mmin, delta_m) + 1)^{-1}
    f(m') = delta_m / m' + delta_m / (m' - delta_m)
    See also, https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
    """
    window_fix = xp.ones_like(masses)
    if delta_m > 0.0:
        smoothing_region = (masses >= mmin) & (masses < (mmin + delta_m))
        shifted_mass = masses[smoothing_region] - mmin
        if shifted_mass.size:
            exponent = xp.nan_to_num(
                delta_m / shifted_mass + delta_m / (shifted_mass - delta_m)
            )
            window_fix[smoothing_region] = 1 / (xp.exp(exponent) + 1)
    window_fix[(masses < mmin) | (masses > mmax)] = 0
    return window_fix    

def extract_mass_parameters_fix(parameters):
    """extract the parameters of the mass distribution hyperparameters used in
    T&T18 from either a list or dictionary."""
    if isinstance(parameters, list):
        return parameters
    elif isinstance(parameters, dict):
        keys = ['alpha', 'mmin',  'delta_m']
        return [parameters[key] for key in keys]

def ppow_fix(masses, parameters):
    """1d unnormalised powerlaw mass probability with smoothed low-mass end"""
    alpha, mmin,  delta_m = extract_mass_parameters_fix(parameters)
    return masses**(-alpha) * window_fix(masses, mmin, delta_m) 

def norm_ppow_fix(parameters):
    """normalise ppow, requires m1s, an array of m values, and dm, the spacing of
    that array"""
    m1s = np.linspace(0.9, 2.9, 500)
    return np.trapz(ppow_fix(m1s, parameters), m1s)

def turn_on_pow_fix(masses, parameters, pow_norm_fix):
    alpha, mmin, delta_m = extract_mass_parameters_fix(parameters)
    p_pow_fix = ppow_fix(masses, parameters) / pow_norm_fix
    return  p_pow_fix 

def hyper_prior_turn_on_pow_fix(dataset, alpha, mmin, delta_m):
    parameters = dict(alpha=alpha, mmin=mmin, delta_m=delta_m)
    pow_norm_fix = norm_ppow_fix(parameters)
    probability_fix = turn_on_pow_fix(dataset['mu'], parameters, pow_norm_fix)
    return probability_fix

hp_priors_turn_on_pow_fix= dict(alpha=Uniform(-5, 25, 'alpha', '$\\alpha$'),
                 mmin=Uniform(0.9, 1.5, 'mmin', '$mmin$'),
                delta_m=Uniform(0.01, 1, 'delta', '$\\delta$'))
#end turn_on_pow fixed m_max

#begin top + Gaussian   
def pow_add_G(masses, alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam):
    """1d unnormalised powerlaw mass probability with smoothed low-mass end"""
    return ( lam*masses**(-alpha) + (1-lam)* np.exp(-(masses - mu_g)**2 / (2 * sigma_g**2))*(mu_g > (mmin + delta_m) )  )\
             *window(masses, mmin, mmax, delta_m)

def sum_pow_add_G(alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam):
    mii=np.linspace(0.9, 2.9, 500)
    return np.trapz(pow_add_G(mii, alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam), mii)

def norm_pow_add_G(masses,alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam):
    return pow_add_G(masses, alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam)\
           /sum_pow_add_G(alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam)

def hyper_prior_turn_on_pow_G(dataset,alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam):
    return norm_pow_add_G(dataset['mu'],alpha, mmin, mmax, delta_m, mu_g,sigma_g,lam)

hp_priors_turn_on_pow_G = dict(alpha=Uniform(-5, 25, 'alpha', '$\\alpha$'),
                 mmin=Uniform(0.9, 1.5, 'mmin', '$mmin$'),
                 mmax=Uniform(1.5, 2.9, 'mmax', '$mmax$'),
                lam=Uniform(0.0, 1, 'lam', '$\\lambda$'),
                 mu_g=Uniform(1.5, 2.9, 'mpp', '$mpp$'),
                 sigma_g=Uniform(0.01, 2, 'sigpp', '$\\sigma$'),
                delta_m=Uniform(0.01, 1, 'delta', '$\\delta$')) 
#end top + Gaussian 

#begin top + Gaussian with fixed max 
def window_fixed_max(masses, mmin, delta_m):

    """
    Apply a one sided window between mmin and mmin + delta_m to the
    mass pdf.
    The upper cut off is a step function,
    the lower cutoff is a logistic rise over delta_m solar masses.
    See T&T18 Eqs 7-8
    Note that there is a sign error in that paper.
    S = (f(m - mmin, delta_m) + 1)^{-1}
    f(m') = delta_m / m' + delta_m / (m' - delta_m)
    See also, https://en.wikipedia.org/wiki/Window_function#Planck-taper_window
    """
    mmax=2.9
    window = xp.ones_like(masses)
    if delta_m > 0.0:
        smoothing_region = (masses >= mmin) & (masses < (mmin + delta_m))
        shifted_mass = masses[smoothing_region] - mmin
        if shifted_mass.size:
            exponent = xp.nan_to_num(
                delta_m / shifted_mass + delta_m / (shifted_mass - delta_m)
            )
            window[smoothing_region] = 1 / (xp.exp(exponent) + 1)
    window[(masses < mmin) | (masses > mmax)] = 0
    return window    

def pow_add_G_fixed_max(masses, alpha, mmin, delta_m, mu_g,sigma_g,lam):
    mmax=2.9
    """1d unnormalised powerlaw mass probability with smoothed low-mass end"""
    return ( lam*masses**(-alpha) + (1-lam)* np.exp(-(masses - mu_g)**2 / (2 * sigma_g**2))*\
             (mu_g > (mmin + delta_m) )  )*window_fixed_max(masses, mmin, delta_m)

def sum_pow_add_G_fixed_max(alpha, mmin, delta_m, mu_g,sigma_g,lam):
    mmax=2.9
    mii=np.linspace(0.9, 2.9, 500)
    return np.trapz(pow_add_G_fixed_max(mii, alpha, mmin, delta_m, mu_g,sigma_g,lam), mii)

def norm_pow_add_G_fixed_max(masses,alpha, mmin, delta_m, mu_g,sigma_g,lam):
    mmax=2.9
    return pow_add_G_fixed_max(masses, alpha, mmin,delta_m, mu_g,sigma_g,lam)\
           /sum_pow_add_G_fixed_max(alpha, mmin, delta_m, mu_g,sigma_g,lam)

def hyper_prior_turn_on_pow_G_fixed_max(dataset,alpha, mmin, delta_m, mu_g,sigma_g,lam):
    mmax=2.9
    return norm_pow_add_G_fixed_max(dataset['mu'],alpha, mmin, delta_m, mu_g,sigma_g,lam)

hp_priors_turn_on_pow_G_fixed_max = dict(alpha=Uniform(-5, 25, 'alpha', '$\\alpha$'),
                 mmin=Uniform(0.9, 1.5, 'mmin', '$mmin$'),
                lam=Uniform(0.0, 1, 'lam', '$\\lambda$'),
                 mu_g=Uniform(1.5, 2.9, 'mpp', '$mpp$'),
                 sigma_g=Uniform(0.01, 2, 'sigpp', '$\\sigma$'),
                delta_m=Uniform(0.01, 1, 'delta', '$\\delta$'))   
#end top + Gaussian with fixed max