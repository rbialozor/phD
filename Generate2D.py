# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 09:39:38 2021

@author: AR
"""
import numpy as np

def Generate2D(h, b, ConcStripNo, ztopfi, zbottomfi, zcfrp, IsCFRP):
    delta_h = h / (ConcStripNo)
    z=[]
    
    for i in range(ConcStripNo):
        z.append((delta_h * i + delta_h * 0.5 - h * 0.5))
    z.append(ztopfi);
    z.append(zbottomfi);
    z.append(ztopfi);
    z.append(zbottomfi);
    if (IsCFRP == True):
        z.append(zcfrp)
        
    AllStripNo = len(z)
    xi_count=[]
    xi_count=[1] * AllStripNo
    lodx_all=np.array(z).reshape(AllStripNo,1)
    xi_count=np.array(xi_count).reshape(AllStripNo,1)   
    lodx_all=np.append(lodx_all,xi_count,axis=1)  
    Zvector=lodx_all
    
    
    return Zvector, AllStripNo