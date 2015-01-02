from ctypes import cdll, c_double
from math import pow, log10
import pymc;
from numpy import empty;

global lib
lib = cdll.LoadLibrary('./libOptimize.dylib')
lib.pyEntry.restype = c_double;
inn = empty(14, dtype=c_double);

minValue = 1E10;

x = empty(13, dtype=object);
x[0] = pymc.Normal('kfb2', mu=log10(0.06), tau=1)
x[1] = pymc.Normal('krb2', mu=log10(6), tau=1)
x[2] = pymc.Normal('f1', mu=-4.0, tau=1.0)
x[3] = pymc.Normal('r1', mu=1.00, tau=1.0)
x[4] = pymc.Normal('f3', mu=-4.0, tau=1.0)
x[5] = pymc.Normal('r3', mu=0.00, tau=1.0)
x[6] = pymc.Normal('internalize', mu=-2.0, tau=1)
x[7] = pymc.Normal('pYInt', mu=-1.0, tau=1)
x[8] = pymc.Normal('kRec', mu=log10(5.8E-2), tau=1)
x[9] = pymc.Normal('kDeg', mu=log10(2.2E-3), tau=1)
x[10] = pymc.TruncatedNormal('fElse', mu=-1.0,        tau=1, a=-6.0, b=0.0)
x[11] = pymc.Normal('AXL',            mu=2.0,         tau=1.0)
x[12] = pymc.Normal('Gas',   mu=-3.0, tau=10)
pD1 = pymc.Bernoulli('pD1', p=0.5)

@pymc.deterministic(plot=True)
def errMeas(x_in = x):
    global inn;
    global minValue;
    for i in range(0,13):
        inn[i] = pow(10, x_in[i]);
        
    inn[13] = pD1;
    value = lib.pyEntry(inn.ctypes);
    
    if value < minValue:
    	minValue = value;
    	print(value);
        
    return value
    
@pymc.potential
def psi_i(x_in = errMeas):
    return -x_in;

if __name__ == '__main__':
    # Instantiate model
    M = pymc.MCMC([x, errMeas, psi_i, pD1],db='pickle')
    # Sample
    # M.use_step_method(pymc.AdaptiveMetropolis, x);
    M.sample(iter=int(2E6), burn=int(1E6), thin=1000)

    pymc.Matplot.plot(M)
    
    
    pymc.Matplot.geweke_plot(pymc.geweke(M), 'png');
