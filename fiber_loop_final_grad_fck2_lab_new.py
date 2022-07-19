# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 11:47:37 2017

@author: Rafal
"""

import numpy as np
from materialstress5 import concretestress,CFRPstress, steelstress,ParabolicConcreteStress,UserParabolicConcreteStress
from Generate2D import Generate2D
from GenerateAreaMatrix import GenerateAreaMatrix
from SolverEngine2D import Solve
from MyMath import LinearInterpolation,LinearInterpolationVector,BMUniformlyLoad,DeflectionIntegral,BMDoubleForceLoad

############################## INPUT DATA

h=0.294
b=0.2
Es=   208333333333.
Econc= (10000*38**(1/3)*1000*1000)
Ecfrp=171900 #MPa
#Ecfrp=230000 #MPa


A=b*h
L=2.8

M1=120001.
P1=-0.001


ndiv=30
deltah=h/(ndiv)
deltahlist=[]
deltahlist=np.arange(0,h,deltah)
zi=deltahlist+deltah*0.5-h*0.5

#odwrócenie listy:
zi=zi[::-1]


Aconcrete=zi*0+b*deltah
Aconcrete=Aconcrete.tolist()

#STRZEMIONA
fis=0.006

#ZBROJENIE GORNE
#####################
topcover=0.025      #
topfi=0.008         #
ntopfi=2            #
#####################

a1topfi=0.5*h-(topcover+fis+0.5*topfi)
Atopfi=ntopfi*3.1415*0.25*topfi**2

#ZBROJENIE DOLNE
bottomcover=0.025
#####################
bottomfi=0.012      #
nbottomfi=4         #
#####################

a1bottomfi=-0.5*h+(bottomcover+fis+0.5*bottomfi)
Abottomfi=nbottomfi*3.1415*0.25*bottomfi**2

#CFRP
#####################
tcfrp=0.0012        #
bcfrp=0.12          #
#####################

Acfrp=tcfrp*bcfrp
zcfrp=-0.5*h-0.5*tcfrp

#print (Atopfi, a1topfi,ntopfi)
#print (Abottomfi, a1bottomfi,nbottomfi)
#print (Acfrp, zcfrp)



zbottomfi=a1bottomfi
ztopfi=a1topfi
ConcStripNo=ndiv
IsCFRP=True
ConcreteMaterialIndex=1

C=[0]
steps=50


############################## MATERIAL STRESS CALCULATION ###################
eslim=0.01

'''
Calculation was moved to materialstress module
'''




def naprezenia(efibvec,Nfib):
    sigmaodepsilon=[]
    #print(Nfib)
    #print(ndiv)
    for m in range(0,Nfib,1):
        if m<ndiv:
            sigmaodepsilon.append(concretestress(efibvec[m],fcki)[0])
            #print('conc',efibvec[m],sigmaodepsilon[m])
        elif (m>=ndiv and m<Nfib-3):
            sigmaodepsilon.append(steelstress(efibvec[m])[0])   
            #print('steel',efibvec[m],sigmaodepsilon[m])
        elif (m>=Nfib-3 and m<Nfib-1):
            sigmaodepsilon.append(concretestress(efibvec[m],fcki)[0])   
            #print('steel',efibvec[m],sigmaodepsilon[m])    
        else :
            sigmaodepsilon.append(CFRPstress(efibvec[m],Ecfrp)[0]) 
            #print('cfrp',efibvec[m],sigmaodepsilon[m])
    return sigmaodepsilon



############################## Stifness Matrix [E]  ##########################
Afib=Aconcrete
Afiball=Afib+[Atopfi]+[Abottomfi]+[-Atopfi]+[-Abottomfi]+[Acfrp]
#print(Afiball)
#print(len(Afiball))
Nfiber=len(Afiball)
s=(Nfiber,Nfiber)
E0=np.zeros(s)
AfibM=np.zeros(s)
#print(E0)


for a in range(0,ndiv,1):
    E0[a,a]=Econc
    AfibM[a,a]=Afib[a]
for a in range(ndiv,ndiv+2,1):
    E0[a,a]=Es
    AfibM[a,a]=Afiball[a]
for a in range(ndiv+2,ndiv+4,1):
    E0[a,a]=Econc
    AfibM[a,a]=Afiball[a]
for a in range(ndiv+4,Nfiber):
    E0[a,a]=Ecfrp
    AfibM[a,a]=Afiball[a]
Emod=E0
#print('Emod')  
#print(Emod)

############################ FIBER GEOMTRY MATRIX l(x) ########################

xi_count=[]
for i in range(Nfiber):
        xi_count.append(1)


zi=deltahlist+deltah*0.5-h*0.5

zi_all=zi[::-1]
zi_all=zi*0+deltahlist+deltah*0.5-h*0.5


yi=np.zeros(ndiv).reshape((ndiv, 1))
zi=zi.reshape((ndiv, 1))
xi=np.zeros(ndiv).reshape((ndiv, 1))+1
#print()


#print (zi)
#print (zi_all)
zi_all=zi_all.tolist()
#print (zi_all)

zi_all=zi_all+[a1topfi]+[a1bottomfi]+[a1topfi]+[a1bottomfi]+[zcfrp]
#print (zi_all)
lodx_all=np.array(zi_all).reshape(Nfiber,1)
xi_count=np.array(xi_count).reshape(Nfiber,1)
#print (xi_count)
lodx_all=np.append(lodx_all,xi_count,axis=1)



G=Generate2D(h, b, ConcStripNo, ztopfi, zbottomfi, zcfrp, IsCFRP)

lodxa=G[0]
#print(lodxa)
lodx=lodx_all



print(np.array_equal(lodxa, lodx, equal_nan=True))
AfibMatrix = GenerateAreaMatrix(h,  b,  ConcStripNo,  Econc,  Es,  Ecfrp,  Atopfi,  Abottomfi,  Acfrp, IsCFRP)[1]
print(np.array_equal(AfibM, AfibMatrix, equal_nan=True))

AfibM=AfibMatrix


Nfiber = G[1]
AllStripNo=G[1]
lodx=G[0]
#lodx=lodx_all



     
    
    
K = Solve(M1, P1, steps, AllStripNo, ConcStripNo, lodx, h,  AfibM, IsCFRP, 0, ConcreteMaterialIndex, 1, 1,C,Ecfrp,L)   
krzywiznalist = K[0]
deltaMhistorylist= K[1]
strain_total= K[2]
sigma_total= K[3]
itera= K[4]
efibAB= K[5]
sigmaAB= K[6]
dNhistory= K[7]
dMhistory  = K[8]
iterahistory= K[9]
fcki= K[10]
uhistorylist= K[11]
eshistorylist= K[12]
ethistorylist= K[13]
echistorylist= K[14]
multiplied_list = [[2 * j for j in x] for i,x in enumerate(deltaMhistorylist)]

'''
 int t = 50;
 IEnumerable<double> iksy = Enumerable.Range(0, t + 1).Select(x => x * GlobalVarClass.HeightGlobal * 0.01 / (t));
 double[] BMList = MyMath.BMUniformlyLoad(iksy.ToArray(), (double)numUD_medqp.Value, GlobalVarClass.HeightGlobal * 0.01);
 double[] Curvatures = MyMath.LinearInterpolationVector(bmomenthistory.ToArray(), curvaturelist.ToArray(), BMList);
 double Deflect = MyMath.DeflectionIntegral(iksy.ToArray(), Curvatures, GlobalVarClass.HeightGlobal * 0.01, 0.5);
   
 '''
 
t=200
iksy = np.arange(0,L, L / t)   
BMList = BMDoubleForceLoad(iksy, max(deltaMhistorylist[len(deltaMhistorylist)-1][:])-0.1, L,1);
#BMList = BMUniformlyLoad(iksy, max(deltaMhistorylist[len(deltaMhistorylist)-1][:])-0.1, L);
Curvatures = LinearInterpolationVector(deltaMhistorylist[len(deltaMhistorylist)-1][:], krzywiznalist[len(krzywiznalist)-1][:], BMList)
Deflect = DeflectionIntegral(iksy, Curvatures, L, 0.5);
print(len(deltaMhistorylist)-1)
print (Deflect)
#print (Curvatures)
print('asfasfsfasfa')   























###############################   PLOTS   #####################################
#print(zi_all)
#print(AfibM)

###############################################################################
#########################   CONCRETE / STEEL / CFRP   #########################
###############################################################################
    
    
#print()
#print()
eee=np.arange(-0.00015,0.0035,0.00002)
if (ConcreteMaterialIndex == 0):  sig=concretestress(eee,fcki) 
    #sigmaodepsilon.Add(Concrete.BilinearConcreteStress(efibvec[m, 0], gammaC)[0])
elif (ConcreteMaterialIndex == 1): sig=ParabolicConcreteStress(eee,fcki) 
    #sigmaodepsilon.Add(Concrete.ParabolicConcreteStress(efibvec[m, 0], gammaC)[0])
elif (ConcreteMaterialIndex == 2): sig=UserParabolicConcreteStress(eee,fcki) 
    #sigmaodepsilon.Add(Concrete.NonlinearEC2(efibvec[m, 0], gammaC)[0]) 
else:                   sig=concretestress(eee,fcki) 
    #sigmaodepsilon.Add(Concrete.ParabolicConcreteStress(efibvec[m, 0], gammaC)[0])


sig=np.array(sig)
#print()
#print()  
#print('concrete behaviour')
#print(eee)
#print(sig)

import matplotlib.pyplot as p
fig, axs = p.subplots(1, 3,figsize=(10,2))
ax = axs[0]

ax.plot(eee*1000,sig,'-',linewidth=2.0, c='blue')

ax.set_xlabel('\u03B5 [‰]',family='Isocpeur',size=14)
ax.set_ylabel('\u03C3 [MPa]',family='Isocpeur',size=14)
ax.set_title('Concrete',family='Isocpeur',size=15)
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)
ax.set_xlim(-0.2, 3.6)


ax = axs[1]
eee=np.arange(-eslim,eslim,0.001)
sig=steelstress(eee) 
ax.plot(eee*1000,sig,'-',linewidth=2.0, c='red')

ax.set_xlabel('\u03B5 [‰]',family='Isocpeur',size=14)
ax.set_ylabel('\u03C3 [MPa]',family='Isocpeur',size=14)
ax.set_title('Steel',family='Isocpeur',size=15)
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)
ax.set_xlim(-eslim*1000, eslim*1000)
ax.set_ylim(float(min(sig)), float(max(sig)))


ax = axs[2]
eee=np.arange(0,0.0055,0.0005)
sig=CFRPstress(eee,Ecfrp) 
ax.plot(eee*1000,sig,'-',linewidth=2.0, c='black')

ax.set_xlabel('\u03B5 [‰]',family='Isocpeur',size=14)
ax.set_ylabel('\u03C3 [MPa]',family='Isocpeur',size=14)
ax.set_title('CFRP',family='Isocpeur',size=15)
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)



p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.5)





###########################################################################
############################## MATERIAL PLOT ##############################
###########################################################################

fig, axs = p.subplots(1, 2,figsize=(12, 3))

for i in range(len(sigma_total)):
    #p.plot(zi_all[0:ndiv],sigma_total[i])
    ax = axs[0]
    ax.plot(sigma_total[i],zi_all[0:ndiv],label=str(i))
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_ylabel('h [m]')
ax.set_xlabel('\u03C3 [MPa]')
ax.set_title('Concrete \u03C3 iterations')

legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})


for i in range(len(sigma_total)):
    #p.plot(zi_all[0:ndiv],sigma_total[i])
    ax = axs[1]
    ax.plot(strain_total[i],zi_all[0:ndiv],label=str(i))
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_ylabel('h [m]')
ax.set_xlabel('\u03B5 [‰]')
ax.set_title('Concrete \u03B5 iterations')
legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),framealpha=0.5, shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})


p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.5)

###########################################################################
###########################################################################
###########################################################################

###########################################################################
###########################################################################

fig, axs = p.subplots(2, 2,figsize=(7,7))


ax = axs[0,1]
ax.plot(sigmaAB[0:ndiv],zi_all[0:ndiv],'-',linewidth=2.0, c='red')
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)
ax.set_ylabel('h [m]')
ax.set_xlabel('\u03C3 [MPa]')
ax.set_title('\u03C3 (z)')




ax = axs[1,1]
ax.plot(efibAB[0:ndiv],zi_all[0:ndiv],'-',linewidth=2.0, c='blue')
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)
ax.set_ylabel('h [m]')
ax.set_xlabel('\u03B5 [‰]')
ax.set_title('\u03B5 (z)')




# Data to plot.
epka=0.30/ndiv
wys=np.arange(0,0.30,epka)
x, y = np.meshgrid(np.arange(0,0.4,0.2), wys)
s = (2,ndiv)

yticks = np.arange(-50, 30, 10)


ax=axs[0,0]
z=np.zeros(s)
for i in np.arange(0,2,1):
    for j in np.arange(0,len(wys),1):
        z[i,j]=(sigmaAB[j])
z=z.transpose()
#print(z)
cs = ax.contourf(x, y, z,cmap="seismic")
ax.contour(cs, colors='k')
ax.grid(c='k', ls='-', alpha=0.3)


ax=axs[1,0]
z=np.zeros(s)
for i in np.arange(0,2,1):
    for j in np.arange(0,len(wys),1):
        z[i,j]=(efibAB[j])
z=z.transpose()
#print(z)
cs = ax.contourf(x, y, z,cmap="seismic")
ax.contour(cs, colors='k')
ax.grid(c='k', ls='-', alpha=0.3)
    

p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.4,
                    wspace=0.5)

###########################################################################
###########################################################################

###########################################################################


fig, axs = p.subplots(1, 2,figsize=(12, 3))
ax = axs[0]
ax.plot(iterahistory,dNhistory)
ax.grid(True, linestyle='--')
ax.set_xlabel('Iter [-]',family='Isocpeur',size=14)
ax.set_ylabel('N [kN]',family='Isocpeur',size=14)
ax.set_title('N convergence',family='oswald',size=14)

for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)

ax = axs[1]



ax.plot(iterahistory,dMhistory,label='Convergence in '+str(itera)+' iterations')


ax.set_xlabel('Iter [-]',family='Isocpeur',size=14)
ax.set_ylabel('M [kNm]',family='Isocpeur',size=14)
ax.grid(True, linestyle='--')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.5)
  
  
ax.set_title('M convergence',family='oswald',size=14)
p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25,
                    wspace=0.25)
legend = ax.legend(loc='upper right',framealpha=0.5, shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '12'})

###########################################################################



#p.figure(figsize=(16,4))
#p.plot(sigmanonlinear[0:5],zi_all[0:5],'-',linewidth=2.0, c='blue')
#p.grid()
rjxhistory1=[]
rjxhistory2=[]
iterahistory=[]
DjRxhistory1=[]
DjRxhistory2=[]
dDjxhistory1=[]
dDjxhistory2=[]




from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from Diagrams import Laboratory
MLAB=Laboratory.A1M
KLAB=Laboratory.A1K
FLAB=Laboratory.A1F
ULAB=Laboratory.A1U

MLAB2=Laboratory.A2M
KLAB2=Laboratory.A2K
FLAB2=Laboratory.A2F
ULAB2=Laboratory.A2U

#print(len(C))
#print(krzywiznalist)
#print(deltaMhistorylist)
fig, axb = p.subplots(1, 2,figsize=(12, 4))
ax=axb[0]
for m in range(0,len(C)):
    ax.plot(krzywiznalist[m],deltaMhistorylist[m],label='C'+str(C[m]),linewidth=2.0)
    #ax.scatter(krzywiznalist[m],deltaMhistorylist[m], s=30,  alpha=0.5)
    print('k',krzywiznalist[m])
    print('m',deltaMhistorylist[m])
ax.grid(True, linestyle='--')

ax.plot(KLAB,MLAB,label='LAB',linewidth=2.0,color='blue')
ax.plot(KLAB2,MLAB2,label='A2LAB',linewidth=2.0,color='black')

for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_xlabel('\u03BA [1/m]',family='Isocpeur',size=14)
ax.set_ylabel('M [kNm]',family='Isocpeur',size=14)
ax.set_title('M - k relationship',size=14,family='oswald')
ax.set_xlim(0, 0.04)
ax.set_ylim(0, float(max(max(deltaMhistorylist[len(C)-1]),MLAB2[len(MLAB2)-1])))

ax.xaxis.set_major_locator(MultipleLocator(0.005))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25,
                    wspace=0.25)
legend = ax.legend(loc="lower right", shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})









ax=axb[1]
for m in range(0,len(C)):
    ax.plot(uhistorylist[m],multiplied_list[m],label='C'+str(C[m]),linewidth=2.0)
    #ax.scatter(krzywiznalist[m],deltaMhistorylist[m], s=30,  alpha=0.5)
ax.grid(True, linestyle='--')

ax.plot(ULAB,FLAB,label='LAB',linewidth=2.0,color='blue')
ax.plot(ULAB2,FLAB2,label='LAB2',linewidth=2.0,color='black')
for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_xlabel('u [m]',family='Isocpeur',size=14)
ax.set_ylabel('F [kN]',family='Isocpeur',size=14)
ax.set_title('F - u relationship',size=14,family='oswald')
ax.set_xlim(0, 0.03)
ax.set_ylim(0, float(max(max(multiplied_list[len(C)-1]),max(FLAB2))))


ax.xaxis.set_major_locator(MultipleLocator(0.01))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25,
                    wspace=0.25)
legend = ax.legend(loc="lower right", shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})



#legend = ax.legend(loc="center left", bbox_to_anchor=(1, 0, 0.5, 1), shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})
#legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})












import xlsxwriter 
# Create an new Excel file and add a worksheet.
workbook = xlsxwriter.Workbook('file.xlsx')
worksheet = workbook.add_worksheet()
for m in range(len(krzywiznalist[0])):
    worksheet.write(m, 3, krzywiznalist[0][m])
    worksheet.write(m, 4, deltaMhistorylist[0][m])    
    worksheet.write(m, 6, uhistorylist[0][m])
    worksheet.write(m, 7, multiplied_list[0][m])    
    worksheet.write(m, 8, eshistorylist[0][m]) 
    worksheet.write(m, 10, ethistorylist[0][m]) 
    worksheet.write(m, 11, echistorylist[0][m]) 
workbook.close()


######################  K(X)  ######################################

fig, axc = p.subplots(1, 2,figsize=(12, 4))
ax=axc[0]
ax.grid(True, linestyle='--')

ax.plot(iksy,Curvatures,label='\u03BA(x)',linewidth=2.0,color='blue')

for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_xlabel('x [m]',family='Isocpeur',size=14)
ax.set_ylabel('\u03BA [1/m]',family='Isocpeur',size=14)
ax.set_title('\u03BA(x) distribution',size=14,family='oswald')
ax.set_xlim(0, L)

ax.xaxis.set_major_locator(MultipleLocator(0.4))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

#ax.set_ylim(0, float(max(Curvatures)))
p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25, wspace=0.25)
legend = ax.legend(loc="lower right", shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})

############################ M(X)######################################


ax=axc[1]
ax.grid(True, linestyle='--')

ax.plot(iksy,BMList,label='M(x)',linewidth=2.0,color='blue')

for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_xlabel('x [m]',family='Isocpeur',size=14)
ax.set_ylabel('M [kNm]',family='Isocpeur',size=14)
ax.set_title('M(x) distribution',size=14,family='oswald')
ax.set_xlim(0, L)


ax.xaxis.set_major_locator(MultipleLocator(0.4))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

#ax.set_ylim(0, float(max(BMList)))
p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25, wspace=0.25)
legend = ax.legend(loc="lower right", shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})

##################################################################





######################  F-es  ######################################

fig, axn = p.subplots(1, 2,figsize=(12, 4))
ax=axn[0]
for m in range(0,len(C)):
    ax.plot(eshistorylist[m],multiplied_list[m],label='C'+str(C[m]),linewidth=2.0,color='r')
ax.grid(True, linestyle='--')
if IsCFRP!=False:
    for m in range(0,len(C)):
        ax.plot(ethistorylist[m],multiplied_list[m],label='CFRP'+str(C[m]),linewidth=2.0,linestyle='-',color='k')

for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_xlabel('\u03B5 [‰]',family='Isocpeur',size=14)
ax.set_ylabel('F [kN]',family='Isocpeur',size=14)
ax.set_title('\u03B5s distribution',size=14,family='oswald')
#ax.set_xlim(0, L)
ax.set_xlim( -0.006,0)
ax.set_ylim(0, 200)

ax.xaxis.set_major_locator(MultipleLocator(0.001))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25, wspace=0.25)
legend = ax.legend(loc="lower right", shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})

##################################################################

######################  F-ec  ######################################

ax=axn[1]
for m in range(0,len(C)):
    ax.plot(echistorylist[m],multiplied_list[m],label='C'+str(C[m]),linewidth=2.0,color='b')
ax.grid(True, linestyle='--')


for axis in ['top','right']:
  ax.spines[axis].set_linewidth(1)
  ax.spines[axis].set_edgecolor('silver')
  ax.spines[axis].set_linestyle('--')
for axis in ['bottom','left']:
  ax.spines[axis].set_linewidth(1.0)
ax.set_xlabel('\u03B5c [‰]',family='Isocpeur',size=14)
ax.set_ylabel('F [kN]',family='Isocpeur',size=14)
ax.set_title('\u03B5c distribution',size=14,family='oswald')
#ax.set_xlim(0, L)
ax.set_xlim( 0,0.002)
ax.set_ylim(0, 200)

ax.xaxis.set_major_locator(MultipleLocator(0.0005))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

p.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.8, hspace=0.25, wspace=0.25)
legend = ax.legend(loc="lower right", shadow=False, fontsize='6',prop={'family': 'Isocpeur','size': '10'})

##################################################################