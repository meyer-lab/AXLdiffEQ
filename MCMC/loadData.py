# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 12:17:14 2014

@author: asm
"""
import cPickle as pickle
import numpy
import csv
from ctypes import cdll, c_double
from math import pow
from numpy import empty;
import pymc;

global lib
lib = cdll.LoadLibrary('./libOptimize.dylib')
lib.pyEntry.restype = c_double;


def errMeas(x_in):
    inn = empty(14, dtype=c_double);
    
    for i in range(0,13):
        inn[i] = pow(10, x_in[i]);
        
    inn[13] = x_in[13];
    value = lib.pyEntry(inn.ctypes);
    return value

fp = open('MCMC.pickle');
data = pickle.load(fp)



x = numpy.ndarray(shape=(15,len(data['f1'][0])))

for i in range(0, len(data['f1'][0])):
    x[0][i] = data['kfb2'][0][i];
    x[1][i] = data['krb2'][0][i];
    x[2][i] = data['f1'][0][i];
    x[3][i] = data['r1'][0][i];
    x[4][i] = data['f3'][0][i];
    x[5][i] = data['r3'][0][i];
    x[6][i] = data['internalize'][0][i];
    x[7][i] = data['pYInt'][0][i];
    x[8][i] = data['kRec'][0][i];
    x[9][i] = data['kDeg'][0][i];
    x[10][i] = data['fElse'][0][i];
    x[11][i] = data['AXL'][0][i];
    x[12][i] = data['Gas'][0][i];
    x[13][i] = data['pD1'][0][i];
    
    y = list(x[...,i]);    
    x[14][i] = 0.0; #errMeas(y[0:14])
    #print(i)
    
    
    
x = numpy.transpose(x);



with open('out.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(x)
