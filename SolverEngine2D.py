# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 09:39:38 2021

@author: AR
"""
import numpy as np
from MyMath import LinearInterpolation,LinearInterpolationVector,BMUniformlyLoad,DeflectionIntegral,BMDoubleForceLoad
from StressVector import SectionStressVector
def Solve(BendingMomentValue, AxialForceValue, steps, AllStripNo, ConcStripNo, lodx, h,  AfibMatrix, IsCFRP, ConcreteClassIndex, ConcreteMaterialIndex, gammaC, gammaS,C,Ecfrp,L):

    iterahistory=[]
    dDjxhistory1=[]
    dDjxhistory2=[]
    sigma_total=[]
    strain_total=[]
    cause=''


    itera=0
    delta_b= 0.011*0.001
    delta_a=delta_b/h 
    NA=0
    MA=0
    a1=0 
    b1=0
    deltaab=[]
    sigmaA=[]
    sigmaB=[]
    dNhistory=[]
    dMhistory=[]

    #print()



    krzywizna=[]
    deltaMhistory=[]
    krzywiznalist=[]
    deltaMhistorylist=[]
    uhistory=[]
    uhistorylist=[]
    eshistory=[]
    eshistorylist=[]
    ethistory=[]
    ethistorylist=[]
    echistory=[]
    echistorylist=[]


    for fcki in C:
        for deltaM in np.arange(0,steps,1):
            

            if (BendingMomentValue != 0):
                    delta_b = 0.0011 * 0.001;
                    delta_a = delta_b / h;
            elif (AxialForceValue == 0):
                    delta_b = 0;
                    delta_a = delta_b / h;
            else:
                    delta_b = 0.0011 * 0.001;
                    delta_a = 0;              
            itera=0

            NA=0
            MA=0
            a1=0 
            b1=0
            deltaab=[]
            sigmaA=[]
            sigmaB=[]
            dNhistory=[]
            dMhistory=[]
            iterahistory=[]
            
            BendingMomentIteration = deltaM * BendingMomentValue / (steps)
            #M1=deltaM*M123/steps
            
            for a in np.arange(1,100,1):

                ez1=np.array([[(a1+delta_a),b1]])
                ez2=np.array([[(a1),b1+delta_b]])
                
                #print('a1', a1)
                #print('b1', b1)
                #print('delta_a', delta_a)
                #print('delta_b', delta_b)
        
                #print('ez1', ez1)
                efib1=np.dot(lodx,np.transpose(ez1))
                #print(efib1)
                efib2=np.dot(lodx,np.transpose(ez2))
        
                sigmaA=SectionStressVector(efib1, AllStripNo, ConcStripNo, ConcreteClassIndex, ConcreteMaterialIndex, gammaC, gammaS, IsCFRP,fcki,Ecfrp);
                sigmaB=SectionStressVector(efib2, AllStripNo, ConcStripNo, ConcreteClassIndex, ConcreteMaterialIndex, gammaC, gammaS, IsCFRP,fcki,Ecfrp);
                #sigmaA=naprezenia(efib1,AllStripNo)
                #sigmaB=naprezenia(efib2,AllStripNo)
                #print('odksztalcenia efib1')
                #print(efib1)
        
                #print( )
                #print('ez2', ez2)
                #print('odksztalcenia efib2')
                #print(efib2)
         
                #print( )
                
                

                #print('naprezenia sigmaA')
                #print(sigmaA)
                for i in range(0,len(sigmaA)):
                    sigmaA[i]=sigmaA[i]*10**6
                MNDA=np.dot(np.transpose(lodx),AfibMatrix)
                #print('MNDA',MNDA)   
                MNDA=np.dot(MNDA,(sigmaA))  
                #print('MNDA',MNDA)   
                #print('naprezenia sigmaA')
                #print(sigmaA)
                for i in range(0,len(sigmaA)):
                    sigmaA[i]=sigmaA[i]*0.1**6
        
        
                    
                    
                for i in range(0,len(sigmaB)):
                    sigmaB[i]=sigmaB[i]*10**6
                MNDB=np.dot(np.transpose(lodx),AfibMatrix)
                MNDB=np.dot(MNDB,(sigmaB))
                #print('MNDB',MNDB)
                #print('naprezenia sigmaB')
                #print(sigmaB)        
                for i in range(0,len(sigmaB)):
                    sigmaB[i]=sigmaB[i]*0.1**6
                    
                #print('NA',NA)  
                #print('MA',MA)  
                #print('P1',P1)  
                #print('M1',M1)  
        
                WW=(MNDA[1]-NA)*(MNDB[0]-MA)-(MNDA[0]-MA)*(MNDB[1]-NA)
                #print('WW',WW)
                WA=(AxialForceValue-NA)*(MNDB[0]-MA)-(BendingMomentIteration-MA)*(MNDB[1]-NA)
                #print('WA',WA)
                WB=(MNDA[1]-NA)*(BendingMomentIteration-MA)-(MNDA[0]-MA)*(AxialForceValue-NA)
                #print('WB',WB)
                if WW != 0 :
                    ZA = WA / WW
                    ZB = WB / WW
                else:
                    ZA = 0
                    ZB = 0
                #print('ZA i ZB')
                #print (ZA)
                #print (ZB)
                a1=a1+ZA*delta_a
                b1=b1+ZB*delta_b
                
                #print('a1 i b1')
                #print (a1)
                #print (b1)
                #print('delta_a i delta_b')
                #print (delta_a)
                #print (delta_b)
                
                
                ezAB=np.array([[a1,b1]])
               # print ('ezAB',ezAB)
                efibAB=np.dot(lodx,np.transpose(ezAB))
                #print('odksztalcenia efibAB')
                #print(efibAB)
                sigmaAB = SectionStressVector(efibAB, AllStripNo, ConcStripNo, ConcreteClassIndex, ConcreteMaterialIndex, gammaC, gammaS, IsCFRP,fcki,Ecfrp)
                #sigmaAB=naprezenia(efibAB,Nfiber)
                for i in range(0,len(sigmaAB)):
                    sigmaAB[i]=sigmaAB[i]*10**6
                #print('naprezenia sigmaAB')
                #print(sigmaAB) 
                MNDAB=np.dot(np.transpose(lodx),AfibMatrix)
                #print(np.transpose(lodx))
                #print(MNDAB)
                mm=0
                kk=0
                mj=0
                for mm in MNDAB[0]:
                    #print(mm,'\t',sigmaAB[kk])
                    mj=mj+mm*(sigmaAB[kk])
                    kk=kk+1
                #print(mj)
                MNDAB=np.dot(MNDAB,(sigmaAB))
                #print(MNDAB)
                for i in range(0,len(sigmaAB)):
                    sigmaAB[i]=sigmaAB[i]*0.1**6
                MA=MNDAB[0] 
                NA=MNDAB[1]   
                
                MN=[]
                MN=[BendingMomentIteration]+[AxialForceValue]
                MN=(np.array((MN)))
                
                DIFF=MN-MNDAB
                dNhistory.append(DIFF[1]*0.001)
                dMhistory.append(DIFF[0]*0.001)
        
                #print(MN,DIFF,itera)
                if abs(DIFF[0])<1 and abs(DIFF[1])<1:
        
                    #print('JEST!!!!!!!!!!!!!!!!!!!',DIFF[0],DIFF[1])
                    #print((DIFF[0])<1)
                    #print((DIFF[1])<1)
                    itera=itera+1
                    iterahistory.append(itera)
                    break
                deltaab.append(b1)
                
                #print('ITERACJA############################################################', itera)
                

                dDjxhistory1.append(DIFF[0])
                dDjxhistory2.append(DIFF[1])
        

                iterahistory.append(itera)
                
                itera=itera+1
        
            strain_total.append(efibAB[0:ConcStripNo])
            sigma_total.append(sigmaAB[0:ConcStripNo])     
            deltaMhistory.append(1*MA*0.001)
            krzywizna.append(a1)
            eshistory.append(efibAB[ConcStripNo+1])
            ethistory.append(efibAB[-1])
            echistory.append(efibAB[ConcStripNo-3])
            
            t=500
            iksy = np.arange(0,L, L / t)   
            BMList = BMUniformlyLoad(iksy, max(deltaMhistory)-0.1, L);
            #BMList = BMDoubleForceLoad(iksy, max(deltaMhistory)-0.1, L,1);
            Curvatures = LinearInterpolationVector(deltaMhistory, krzywizna, BMList)
            Deflect = DeflectionIntegral(iksy, Curvatures, L, 0.5);
            uhistory.append(Deflect)
            
            if float(max((efibAB[0:ConcStripNo])))>0.0035:
                print('Concrete ',float(max((efibAB[0:ConcStripNo]))))
                print('Steel ',float(max(abs(efibAB[ConcStripNo:AllStripNo-3]))))
                if IsCFRP == True: 
                    print('CFRP ',max(abs(efibAB[AllStripNo-1])))
                krzywiznalist.append(krzywizna)
                deltaMhistorylist.append(deltaMhistory)
                uhistorylist.append(uhistory)
                eshistorylist.append(eshistory)
                ethistorylist.append(ethistory)
                echistorylist.append(echistory)
                krzywizna=[]
                deltaMhistory=[]
                uhistory=[]
                eshistory=[]
                ethistory=[]
                echistory=[]
                cause = 'concrete'
                break
            if float(max(abs(efibAB[ConcStripNo:AllStripNo-3])))>0.005:
                print('Steel ',float(max(abs(efibAB[ConcStripNo:AllStripNo-3]))))
                if IsCFRP == True: 
                    print('CFRP ',max(abs(efibAB[AllStripNo-1])))
                print('Concrete ',float(max((efibAB[0:ConcStripNo]))))
                krzywiznalist.append(krzywizna)
                deltaMhistorylist.append(deltaMhistory)
                uhistorylist.append(uhistory)
                eshistorylist.append(eshistory)
                ethistorylist.append(ethistory)
                echistorylist.append(echistory)
                krzywizna=[]
                deltaMhistory=[]
                uhistory=[]
                eshistory=[]
                ethistory=[]
                echistory=[]
                cause = 'steel'
                break
            
            if IsCFRP == True:
                if max(abs(efibAB[AllStripNo-1]))>0.0058:
                    print('CFRP ',max(abs(efibAB[AllStripNo-1])))
                    print('Steel ',float(max(abs(efibAB[ConcStripNo:AllStripNo-3]))))
                    print('Concrete ',float(max((efibAB[0:ConcStripNo]))))
                    krzywiznalist.append(krzywizna)
                    deltaMhistorylist.append(deltaMhistory)
                    uhistorylist.append(uhistory)
                    eshistorylist.append(eshistory)
                    ethistorylist.append(ethistory)
                    echistorylist.append(echistory)
                    krzywizna=[]
                    deltaMhistory=[]
                    uhistory=[]
                    eshistory=[]
                    ethistory=[]
                    echistory=[]
                    cause = 'CFRP'
                    break   
            print('deltaMhistory solver: ',len(deltaMhistory),'cause:',cause )
        print()
        print('history - aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa') 
        print(itera)   
        print(dNhistory[itera-1])     
        print(dMhistory[itera-1])  
        print()
        print('###################################################################')
        print()
        print('\u03B5 CFRP [‰] \t \u03C3 CFRP [MPa]') 
        print(round(float(efibAB[AllStripNo-1]*1000),3),'\t \t',round(sigmaAB[AllStripNo-1],3))
        print()
        print('\u03B5 BT STEEL[‰] \t \u03C3 BT STEEL [MPa]') 
        print(round(float(efibAB[AllStripNo-4]*1000),3),'\t \t',round(sigmaAB[AllStripNo-4],3))
        print()
        print('\u03B5 TOP STEEL[‰] \t \u03C3 TOP STEEL [MPa]') 
        print(round(float(efibAB[AllStripNo-5]*1000),3),'\t \t',round(sigmaAB[AllStripNo-5],3))
        print()
        print('\u03B5 TP CONC [‰] \t \u03C3 TP CONC [MPa]') 
        print(round(float(efibAB[ConcStripNo-1]*1000),3),'\t \t',round(sigmaAB[ConcStripNo-1],3))
        print()
        print('N [kN] \t \t M [kNm]') 
        print(round(float(NA*0.001),3),'\t \t',round(MA*0.001,3),BendingMomentIteration)
        print()
        print('###################################################################')
        #print(np.dot(kodx,(np.transpose(eps))))
        
    return krzywiznalist,deltaMhistorylist,strain_total,sigma_total,itera,efibAB,sigmaAB,dNhistory,dMhistory,iterahistory,fcki,uhistorylist,eshistorylist,ethistorylist,echistorylist