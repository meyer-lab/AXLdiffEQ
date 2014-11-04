from ctypes import cdll, c_double
import numpy as np
import collections
from math import log10
import cPickle as pickle
import string
import random
import copy
from scipy import stats
import time

Parameter = collections.namedtuple('Parameter', 'name value')

global lib
lib = cdll.LoadLibrary('./libOptimize.dylib')
lib.pyEntry.restype = c_double;
lib.pyDiffTPS_Activation.restype = c_double;

def likelihood(mcmc, position):    
    if position[10] > 0:
        return 100000000.0;    
    
    pIn = np.power(10,np.array(position))

    value = lib.pyEntry(pIn.ctypes)
    lkh = -stats.chi2.logpdf(value, 24);
    
    # Spatial start
    dIn = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    spatAct = lib.pyDiffTPS_Activation(c_double(16.0), pIn.ctypes, c_double(30), dIn.ctypes)
    
    if spatAct < 0.01:
        spatAct = 0.01;

    return lkh + (100 / spatAct)

def step(mcmc):
    if mcmc.iter % 40 == 0 and mcmc.iter > 100:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)
        
    if time.time() - mcmc.options.model.saveTime > 60*60:
        mcmcSaver = copy.deepcopy(mcmc);
        mcmcSaver.prune(0, 10);
    
        print("Pickle!");
        pickle.dump(mcmcSaver, open(mcmc.options.model.uniqueName,'wb'))
        mcmc.options.model.saveTime = time.time()
        

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def prior(mcmc, position):
    """TODO ...mean and variance from calculate_prior_stats"""
    
    priorVal = (position[0] - log10(0.06)) ** 2 # kfb2
    priorVal += (position[1] - log10(6)) ** 2 # krb2
    
    for x in xrange(2, 5):
        priorVal += (position[x] + 1) ** 2 / 10;
     
    #priorVal += max(position[2] - position[4], 0) ** 10;
    #priorVal += max(position[5] - position[3], 0) ** 10;
    
    priorVal += (position[6] - log10(3E-2)) ** 2 # kint1
    priorVal += (position[7] - log10(3E-1)) ** 2 # kint2
    priorVal += (position[8] - log10(5.8E-2)) ** 2 # kRec
    priorVal += (position[9] - log10(2.2E-3)) ** 2 # kDeg
    
    priorVal += max(position[10], 0) ** 10;    #felse positive
    
    priorVal += (position[10] - log10(0.1)) ** 2 # fElse
    priorVal += (position[11] - log10(100)) ** 2 # AXL    
    priorVal += (position[12] - log10(0.001)) ** 2 # Gas6
    
    return priorVal

class Model(object):
    def __init__(self):
        
        self.parameters = [None] * 13
        self.parameters[0] = Parameter('kfb2', 0.6)
        self.parameters[1] = Parameter('krb2', 6.0)
        self.parameters[2] = Parameter('f1', 0.00081)
        self.parameters[3] = Parameter('r1', 0.34571)
        self.parameters[4] = Parameter('f3', 0.0010493)
        self.parameters[5] = Parameter('r3', 0.017322)
        
        self.parameters[6] = Parameter('internalize', 0.03)
        self.parameters[7] = Parameter('pYint', 0.3)
        self.parameters[8] = Parameter('kRec', 5.8E-2)
        self.parameters[9] = Parameter('kDeg', 2.2E-3)
        self.parameters[10] = Parameter('fElse', 0.01)
        
        self.parameters[11] = Parameter('AXL', 1000)
        self.parameters[12] = Parameter('Gas', 0.001)

        self.uniqueName = 'SPAT'.join(id_generator());
        self.saveTime = time.time()
        
        
    def params(self):
        out = []
        for p in self.parameters:
            out.append(p)

        return out
