# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 09:39:38 2021

@author: AR
"""
import numpy as np
def LinearInterpolation(x, y, xval):
#def LinearInterpolation(double[] x, double[] y, double xval):
        #print(x)
        #print(y)
        #print(xval)
        yval = 0.0;
        for i in range (0, len(x) - 1,1):
                if (xval >= x[i]) and (xval < round(x[i + 1],15)):
                    #print(x[i],xval,round(x[i + 1],4))
                    yval = y[i] + (xval - x[i]) * (y[i + 1] - y[i]) / ((x[i + 1] - x[i]) * 1.0);
                    #print(yval)
        return yval;


def LinearInterpolationVector( x, y, xvals):
#def LinearInterpolationVector(double[] x, double[] y, double[] xvals)
    maxval=0.0;
    yvals = [0] * len(xvals)
    for i in range (0, len(xvals),1):
            temp = LinearInterpolation(x, y, xvals[i])
            if (yvals[i] > maxval):     
                maxval = yvals[i]
            yvals[i] = temp;
    return yvals;
        


def BMUniformlyLoad(xvals, BM, L):
#def BMUniformlyLoad(double[] xvals, double BM, double L)
    q = 8 * BM / (L * L);
    BMval = [0] * len(xvals)
    for i in range (0, len(xvals),1):
        BMval[i] = (q * L * 0.5 * xvals[i] - q * xvals[i] * xvals[i] * 0.5)
    return BMval;

def BMDoubleForceLoad(xvals, BM, L, b):
#def BMUniformlyLoad(double[] xvals, double BM, double L)
    BMval = [0] * len(xvals)
    for i in range (0, len(xvals),1):
        if (xvals[i] < b ):
            BMval[i] = xvals[i] * (BM/b)
        elif(xvals[i] > b  and xvals[i]< L-b):
            BMval[i] = ((BM))
        else:
            BMval[i] = ((BM) - (xvals[i]-(L-b))* (BM/b))
    return (BMval);
      

def DeflectionIntegral(xvals, curvals, L, k):
#def DeflectionIntegral(double[] xvals, double[] curvals, double L, double k)
    #k = 0.5;
    Def = 0.0;
    Deflect = [0] * len(xvals)
    yj = 0.0;
    taj = 0.0;
    for i in range (0, len(xvals),1):   
        Ai = (curvals[i] + curvals[i - 1]) * 0.5 * (xvals[i] - xvals[i - 1])
        x = xvals[i] - 0.5 * (xvals[i] - xvals[i - 1])
        yj = yj + 1.08*(L - x) * Ai
        if (k * L > x):
                taj = taj + Ai * (k * L - x)        
    Def = yj * (k * L) / L - taj
    return Def
