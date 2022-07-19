# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 09:39:38 2021

@author: AR
"""
import numpy as np
from materialstress5 import concretestress,CFRPstress, UserSteelStress, steelstress, steelstress,ParabolicConcreteStress,UserParabolicConcreteStress

def SectionStressVector(efibvec, Nfib, ConcStripNo, ConcreteClassIndex, ModelIndex, gammaC,gammaS, IsCFRP,fcki,Ecfrp):
            sigmaodepsilon=[]
            for m in range(0,Nfib,1):
            
                if (m < ConcStripNo):
                    if (ModelIndex == 0):  sigmaodepsilon.append(concretestress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.BilinearConcreteStress(efibvec[m, 0], gammaC)[0])
                    elif (ModelIndex == 1): sigmaodepsilon.append(ParabolicConcreteStress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.ParabolicConcreteStress(efibvec[m, 0], gammaC)[0])
                    elif (ModelIndex == 2): sigmaodepsilon.append(UserParabolicConcreteStress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.NonlinearEC2(efibvec[m, 0], gammaC)[0]) 
                    else:                   sigmaodepsilon.append(concretestress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.ParabolicConcreteStress(efibvec[m, 0], gammaC)[0])
                
                elif (m >= ConcStripNo and m < ConcStripNo +2):
                    if (ModelIndex == 2):
                        sigmaodepsilon.append(UserSteelStress(efibvec[m])[0]) 
                    else:
                        sigmaodepsilon.append(steelstress(efibvec[m])[0]) 
                    #sigmaodepsilon.Add(Steel.SteelStress(efibvec[m, 0], gammaS)[0])
                
                elif (m >= ConcStripNo + 2 and m < ConcStripNo + 4):
                
                    if (ModelIndex == 0):  sigmaodepsilon.append(concretestress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.BilinearConcreteStress(efibvec[m, 0], gammaC)[0])
                    elif (ModelIndex == 1): sigmaodepsilon.append(ParabolicConcreteStress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.ParabolicConcreteStress(efibvec[m, 0], gammaC)[0])
                    elif (ModelIndex == 2): sigmaodepsilon.append(UserParabolicConcreteStress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.NonlinearEC2(efibvec[m, 0], gammaC)[0]) 
                    else:                   sigmaodepsilon.append(concretestress(efibvec[m],fcki)[0])
                        #sigmaodepsilon.Add(Concrete.ParabolicConcreteStress(efibvec[m, 0], gammaC)[0])
                
                else:
                    if (IsCFRP == True): sigmaodepsilon.append(CFRPstress(efibvec[m],Ecfrp)[0]) 
                        #sigmaodepsilon.Add(CFRP.CFRPStress(efibvec[m, 0])[0])
                
            return sigmaodepsilon; 
       
        
       
