
import numpy as np

def ParabolicConcreteStress(x,y):
    sigma=[]
    fck=y
    C=[0,12,16,20,25,30,35,40,45,50]
    ec1list=[-0.0019,-0.0019,-0.0020,-0.0021,-0.0022,-0.0023,-0.0023,-0.0024,-0.0025,-0.0026]
    indexfck=C.index(fck)
    if y==0:
        fck=38.1
    ec1=ec1list[indexfck]
    
    k=0.0
    
    ec2 = ec1*(k+1)
    
    eclim= -0.0035
    
    fctm=2.37
    Ecm=24100
    n = 2
    for i in range(0,len(x)):
        ain = -x[i]
        if (ain <= 0):
            if (ain < eclim):
                sigma.append(0)
            elif (ain > ec2):
                fiberstress = (fck / 1) *(1 - ((1 - (ain / (ec2)))**n));
                sigma.append(fiberstress)
    
            else:
                fiberstress = (fck / 1)
                sigma.append(fiberstress)
        else:     
                ''' 
                if ain<=0.00015:
                    if Ecm*ain<0.9*fctm:
                        sigma.append(-Ecm*ain)
                    else:
                        #print('tension',ain,'oooooo')
                        
                        sigma.append(-fctm*(1-0.1*((0.00015-ain)/(0.00015-0.9*(fctm*Ecm)))))
                        #print(ain,-fctm*(1-0.1*((0.00015-ain)/(0.00015-0.9*(fctm*Ecm)))))
                else:
                    if ain<0.002:
                        sigma.append(fctm/(0.002-0.00015)*ain-fctm-fctm/(0.002-0.00015)*0.00015)
                    else:
                        sigma.append(0)
                        '''      
                if (fctm / (Ecm/(1+k)) < ain):
                    if ain<0.001:
                        #sigma.append(fctm/(0.001-fctm / (Ecm/(1+k)))*ain-fctm-fctm/(0.001-fctm / (Ecm/(1+k)))*fctm / (Ecm/(1+k)))
                        sigma.append(0);
                    else:
                        sigma.append(0);
                else :
                    #fiberstress = -Ecm*1/(1+k) * ain;
                    #sigma.append(fiberstress);
                    sigma.append(0);

                
    return sigma  

#########################################################################################
def UserParabolicConcreteStress(x,y):
    sigma=[]
    fck=28.4

    ec1=-0.0022
    
    k=1
    
    ec2 = ec1*k
    eclim= -0.0035
    fctm=2.73
    Ecm=23400
    n = 2
    for i in range(0,len(x)):
        ain = -x[i]
        if (ain <= 0):
            if (ain < eclim):
                sigma.append(0)
            elif (ain > ec2):
                fiberstress = (fck / 1) *(1 - ((1 - (ain / (ec2)))**n));
                sigma.append(fiberstress)
    
            else:
                fiberstress = (fck / 1)
                sigma.append(fiberstress)
        else:     
                
                if ain<=fctm / Ecm :
                    if Ecm*ain<0.9*fctm:
                        sigma.append(-Ecm*ain)
                    else:
                        #print('tension',ain,'oooooo')
                        
                        sigma.append(-fctm*(1-0.1*((fctm / Ecm-ain)/(fctm / Ecm-0.9*(fctm*Ecm)))))
                        #print(ain,-fctm*(1-0.1*((0.00015-ain)/(0.00015-0.9*(fctm*Ecm)))))
                else:
                    if ain<0.0012:
                        sigma.append(fctm/(0.0012-fctm / Ecm)*ain-fctm-fctm/(0.0012-fctm / Ecm)*fctm / Ecm)
                        #sigma.append(0)
                    else:
                        sigma.append(0)
                '''            
                if (fctm / Ecm < ain):
                    sigma.append(0);
                else :
                    fiberstress = -Ecm*1/(k) * ain;
                    sigma.append(fiberstress);
                '''
                
    return sigma     

    
#########################################################################################

def concreteparameters(y):
    eclim=-0.0035
    fck=y
    df=8
    Ec0=21500
    
    # Mean value of tensile strength
    fctm=0.3*(fck)**(2/3)
    # Mean value of compressive strength
    fcm=fck+df
    
    
    # The modulus of elasticity in MPa at concrete age of 28 days
    Ecm=Ec0*(fcm/10)**(1/3)
    # A reduced modulus of elasticity
    alphai=min(0.8+0.2*(fcm/88),1)
    Ec=alphai*Ecm
    
    C=[12,16,20,25,30,35,40,45,50]
    ec1list=[-0.0019,-0.0020,-0.0021,-0.0022,-0.0023,-0.0023,-0.0024,-0.0025,-0.0026]
    indexfck=C.index(fck)
    ec1=ec1list[indexfck]
    Ec1=-fcm/ec1
    k=Ecm/Ec1
    
    return Ecm


#########################################################################################
def concretestress(x,y):
    sigma=[]
    if y==0:
        userdata = True
    else:
        userdata=False
        
    if userdata==True:
        fck=38.1  
        Ec0=24100
        fctm=2.04
        fcm=38.1
        Eci=24100
        Ecm=Ec0*(fcm/10)**(1/3)
        eclim=-0.0035
        ec1=-0.0019
    else:    
    #ec1=-0.0023
        eclim=-0.0035
        fck=y
        df=8
        Ec0=21500
        fctm=0.3*(fck)**(2/3)
        fcm=fck+df
        
        alphai=min(0.8+0.2*(fcm/88),1)
        Ecm=Ec0*(fcm/10)**(1/3)
        Ec=alphai*Ecm
        C=[12,16,20,25,30,35,40,45,50]
        ec1list=[-0.0019,-0.0020,-0.0021,-0.0022,-0.0023,-0.0023,-0.0024,-0.0025,-0.0026]
        indexfck=C.index(fck)
        ec1=ec1list[indexfck]
    
    eta=[]
    factor=1
    Ec1=-fcm/(ec1*factor)
    k=1.05*Ecm/Ec1
    
    
    #print(indexfck)
    #print(Ec1)
    #print(k)
    #ec1=ec1list(indexfck)
    c1=3
    c2=6.93
    Gf=0.073*fcm**(0.18)
    wc=5.14*Gf/fctm
    #print(Gf,Ecm,fcm)
    #print(eta)
    #print(x)
    for i in range(0,len(x)):
        ain=-x[i]
        #print(ain)
        if ain<=0:
            eta.append(ain/ec1)
            #print('compression', i)
            if ain<eclim:
                sigma.append(0)
            else: 
                if ain>ec1*factor:
                    fiberstress=((k*eta[i]-(eta[i]**2))/(1+(k-2)*eta[i]))*fcm
                    #fiberstress=-Ec1*ain
                    sigma.append(fiberstress)
                else:
                    fiberstress=fcm-fcm*(ec1+ain)
                    #fiberstress=((k*eta[i]-(eta[i]**2))/(1+(k-2)*eta[i]))*fcm
                    sigma.append(fiberstress)
        else:
            eta.append(-ain/ec1)
            #print('tension',ain)
            if ain<=0.00015:
                if Ecm*ain<0.9*fctm:
                    sigma.append(-Ecm*ain)
                else:
                    #print('tension',ain,'oooooo')
                    
                    sigma.append(-fctm*(1-0.1*((0.00015-ain)/(0.00015-0.9*(fctm*Ecm)))))
                    #print(ain,-fctm*(1-0.1*((0.00015-ain)/(0.00015-0.9*(fctm*Ecm)))))
            else:
                if ain<0.002:
                    sigma.append(fctm/(0.002-0.00015)*ain-fctm-fctm/(0.002-0.00015)*0.00015)
                else:
                    sigma.append(0)
                #print('c',ain,fctm/(0.008-0.00015)*ain-fctm-fctm/(0.008-0.00015)*0.0001)
                #sigma.append(0)
                #sigma.append=((1+(c1*w/wc))*np.exp(-c2*(w/wc))-(w/wc)*(1+c1**(3))*np.exp(-c2))
            
    #print(sigma)         
    return sigma



def steelstress(x):
    sigma=[]
    #fy=270.96774194
    fy=575
    es=0.0029
    Ess=(fy*10**6/(es))/10**6
    Et=1500
    #print(x)
    for i in range(0,len(x)):
        epsilon=round(float(x[i]),8)
        
        if epsilon>=0.0:
            if (epsilon)>=es:
                #print(epsilon,0)
                sigma.append(((fy)+Et*(epsilon-es)))
            else:
                #print(epsilon,1)
                sigma.append((Ess*epsilon))
        else:
            if (epsilon )<-es:
                #print(epsilon,2)
                sigma.append((-(fy)+Et*(epsilon+es)))
            else: 
                #print(epsilon,3)
                sigma.append((Ess*epsilon))
                
    return sigma

def UserSteelStress(x):
    sigma=[]
    fy=421
    es=0.00211
    Ess=(fy*10**6/(es))/10**6
    Et=1500
    #print(x)
    for i in range(0,len(x)):
        epsilon=round(float(x[i]),8)
        
        if epsilon>=0.0:
            if (epsilon)>=es:
                #print(epsilon,0)
                sigma.append(((fy)+Et*(epsilon-es)))
            else:
                #print(epsilon,1)
                sigma.append((Ess*epsilon))
        else:
            if (epsilon )<-es:
                #print(epsilon,2)
                sigma.append((-(fy)+Et*(epsilon+es)))
            else: 
                #print(epsilon,3)
                sigma.append((Ess*epsilon))
                
    return sigma


def CFRPstress(x,y):
    sigma=[]
    Ecfrp=y
    
    for i in range(0,len(x)):
        if x[i]>=0:
            if abs(x[i])>0.0055:
                sigma.append(0)
            else: 
                sigma.append(Ecfrp*x[i])
        else:
                sigma.append(Ecfrp*x[i])

            #sigma.append(float(Emod[i,i])*float(x[i])*0.000001)
    return sigma