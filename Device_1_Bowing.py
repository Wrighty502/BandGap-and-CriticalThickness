import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option("display.max.columns",10)

pi=math.pi

bA = 4          #Slip distnace (Magnititude of Burger's Vector) b in angstrograms
b=bA*(10**-10)

##Cladding
ABCD_E1 = "AlGaAsSb-Cladding"
E1_x=0.9
E1_y=0.07

AC_E1 = "AlAs"
AD_E1 = "AlSb"
BC_E1 = "GaAs"
BD_E1 = "GaSb"
Bowing_E1 = 0.48

##Barrier
ABCD_E2 = "AlGaAsSb-Barrier"
E2_x=0.35
E2_y=0.03

AC_E2 = "AlAs"
AD_E2 = "AlSb"
BC_E2 = "GaAs"
BD_E2 = "GaSb"
Bowing_E2 = 0.48

##QW
ABCD_E3 = "InGaAsSb-Well"
E3_x=0.55
E3_y=0.22

AC_E3 = "InAs"
AD_E3 = "InSb"
BC_E3 = "GaAs"
BD_E3 = "GaSb"
Bowing_E3 = 0.75

#Cap
AB_E0 = "GaSb"

NumberOFWells=10




def bandgap(a0,a,C_11,C_12,av,ac,D0,b_,EV_av,Eg): ## Calculates band gap using Krijn model

    

    aPara=a0
    D001=2*(C_12/C_11)
    aPerp=a*(1-D001*((aPara/a)-1))
    ePerp=(aPerp/a)-1
    ePara=(aPara/a)-1

    EHy_v=av*(2*ePara+ePerp)
    EHy_c=ac*(2*ePara+ePerp)
    
    dE=2*b_*(ePerp-ePara)
    ESH_hh=-0.5*dE
    ESH_lh=(-0.5*D0)+(0.25*dE)+(0.5*(D0**2+D0*dE+(9/4)*dE**2)**0.5)


    EV=EV_av+D0/3+EHy_v+max(ESH_hh,ESH_lh)
    EC=EV_av+D0/3+Eg+EHy_c
    return EV, EC



def QuarternaryBowing(x,y,AC,AD,BC,BD,C):
    QRT_ans = x*y*AC+x*(1-y)*AD +(1-x)*y*BC+(1-x)*(1-y)*BD-x*y*(1-y)*C-(1-x)*y*(1-y)*C-x*(1-x)*y*C-x*(1-x)*(1-y)*C
    return QRT_ans

def Quarternary(x,y,AC,AD,BC,BD):
    QRT_ans = x*y*AC+x*(1-y)*AD +(1-x)*y*BC+(1-x)*(1-y)*BD
    return QRT_ans



def CriticalLayer(E_x,E_y,AC_E,AD_E,BC_E,BD_E):


    ##Epilayer
    x=E_x
    y=E_y            

    C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E]),(Table.loc["C11",AD_E]),(Table.loc["C11",BC_E]),(Table.loc["C11",BD_E]))
    C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E]),(Table.loc["C12",AD_E]),(Table.loc["C12",BC_E]),(Table.loc["C12",BD_E]))
    LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E]),(Table.loc["Lattice Constant",AD_E]),(Table.loc["Lattice Constant",BC_E]),(Table.loc["Lattice Constant",BD_E]))


    ##Substrate
    C11_S_alloy = (Table.loc["C11",AB_E0])
    C12_S_alloy = (Table.loc["C12",AB_E0])
    LatCon_S_alloy = (Table.loc["Lattice Constant",AB_E0])


    v=C12_alloy/(C12_alloy+C11_alloy)
    #print("Poission Ratio: \n" + str(v) + "\n")


    A=LatCon_alloy
    a=A*(10**-10)
    #print("Lattice Constant Epilayer: \n"+ str(A) + "\n")
    #print("Lattice Constant Substrate: \n"+ str(LatCon_S_alloy) + "\n")


    f=abs((A-LatCon_S_alloy)/(LatCon_S_alloy))
    #print("Missfit Percentage: \n"+ str(f*100) + " \n")

    CritThickness = 50


    for i in range (11):  ##More then 12 means Crit-Thickness is negitive for any initial
        CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
        #print(CritThickness)
    CritThickness = CritThickness*10**9

    

    return CritThickness

def bandgapallocation(AC_E,AD_E,BC_E,BD_E,x,y,C):
    
    C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E]),(Table.loc["C11",AD_E]),(Table.loc["C11",BC_E]),(Table.loc["C11",BD_E]))
    C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E]),(Table.loc["C12",AD_E]),(Table.loc["C12",BC_E]),(Table.loc["C12",BD_E]))
    LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E]),(Table.loc["Lattice Constant",AD_E]),(Table.loc["Lattice Constant",BC_E]),(Table.loc["Lattice Constant",BD_E]))
    LatCon_S_alloy = (Table.loc["Lattice Constant",AB_E0])
    
    av=Quarternary(x,y,(Table.loc["av",AC_E]),(Table.loc["av",AD_E]),(Table.loc["av",BC_E]),(Table.loc["av",BD_E]))
    ac=Quarternary(x,y,(Table.loc["ac",AC_E]),(Table.loc["ac",AD_E]),(Table.loc["ac",BC_E]),(Table.loc["ac",BD_E]))
    D0=QuarternaryBowing(x,y,(Table.loc["D0",AC_E]),(Table.loc["D0",AD_E]),(Table.loc["D0",BC_E]),(Table.loc["D0",BD_E]),C)
    b_=Quarternary(x,y,(Table.loc["b_",AC_E]),(Table.loc["b_",AD_E]),(Table.loc["b_",BC_E]),(Table.loc["b_",BD_E]))
    Ev_av=Quarternary(x,y,(Table.loc["Ev_av",AC_E]),(Table.loc["Ev_av",AD_E]),(Table.loc["Ev_av",BC_E]),(Table.loc["Ev_av",BD_E]))
    Eg=QuarternaryBowing(x,y,(Table.loc["Eg",AC_E]),(Table.loc["Eg",AD_E]),(Table.loc["Eg",BC_E]),(Table.loc["Eg",BD_E]),C)


    EV_Epi,EC_Epi= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)
    return EV_Epi,EC_Epi





                ##([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",  "b_",   "Ev-av",   "Eg"])
InAsTable = pd.Series([     6.0583,         832.9,  452.6,   1.00,  -5.08,   0.39,  -1.8,     -0.59,  0.417])
GaAsTable = pd.Series([    5.65325,          1221,  566.0,   1.16,  -7.17,   0.34,  -2.0,     -0.80,  1.420])
InSbTable = pd.Series([     6.4794,         684.7,  373.5,   0.36,  -6.94,   0.81,  -2.0,     -0.00,  0.235])
AlASTable = pd.Series([     5.6611,          1250,  534.0,   2.47,  -5.64,   0.28,  -2.3,     -1.33,  3.099])
GaSbTable = pd.Series([     6.0959,         884.2,  402.6,    0.8,  -7.50,   0.76,  -2.0,     -0.03,  0.812])
AlSbTable = pd.Series([     6.1355,         876.9,  434.1,    1.4,  -4.50,  0.676, -1.35,     -0.41,  2.386])

Table=pd.DataFrame({"InAs":InAsTable,"GaAs":GaAsTable,"InSb":InSbTable,"AlAs":AlASTable,"GaSb":GaSbTable,"AlSb":AlSbTable})
Table.index =([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",   "b_",   "Ev_av",   "Eg"])
print("\n")
#print(Table)
print("\n")







##Epilayer
x=E1_x
y=E1_y
C=Bowing_E1

C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E1]),(Table.loc["C11",AD_E1]),(Table.loc["C11",BC_E1]),(Table.loc["C11",BD_E1]))
C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E1]),(Table.loc["C12",AD_E1]),(Table.loc["C12",BC_E1]),(Table.loc["C12",BD_E1]))
LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E1]),(Table.loc["Lattice Constant",AD_E1]),(Table.loc["Lattice Constant",BC_E1]),(Table.loc["Lattice Constant",BD_E1]))


##Substrate-GaSb


C11_S_alloy = (Table.loc["C11",AB_E0])
C12_S_alloy = (Table.loc["C12",AB_E0])
LatCon_S_alloy = (Table.loc["Lattice Constant",AB_E0])



v=C12_alloy/(C12_alloy+C11_alloy)
#print("Poission Ratio: \n" + str(v) + "\n")


A=LatCon_alloy
a=A*(10**-10)
#print("Lattice Constant Epilayer: \n"+ str(A) + "\n")
#print("Lattice Constant Substrate: \n"+ str(LatCon_S_alloy) + "\n")


f=abs((A-LatCon_S_alloy)/(LatCon_S_alloy))
#print("Missfit Percentage: \n"+ str(f*100) + " \n")

CritThickness = 50


for i in range (11):  ##More then 12 means Crit-Thickness is negitive for any initial
    CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
    #print(CritThickness)
CritThickness_E0_E1 = CritThickness*10**9

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E1 + ": " + str(CritThickness_E0_E1) +"nm")






##Epilayer
E_x=E2_x
E_y=E2_y
AC_E=AC_E2
AD_E=AD_E2
BC_E=BC_E2
BD_E=BD_E2



Critlayer_1=CriticalLayer(E_x,E_y,AC_E,AD_E,BC_E,BD_E)

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E2 + ": " + str(Critlayer_1) +"nm")



##Epilayer
E_x=E3_x
E_y=E3_y
AC_E=AC_E3
AD_E=AD_E3
BC_E=BC_E3
BD_E=BD_E3


Critlayer_2=CriticalLayer(E_x,E_y,AC_E,AD_E,BC_E,BD_E)

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E3 + ": " + str(Critlayer_2) +"nm")


##Epilayer
E_x=E2_x
E_y=E2_y
AC_E=AC_E2
AD_E=AD_E2
BC_E=BC_E2
BD_E=BD_E2


Critlayer_3=CriticalLayer(E_x,E_y,AC_E,AD_E,BC_E,BD_E)

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E2 + ": " + str(Critlayer_3) +"nm")




##Epilayer
E_x=E1_x
E_y=E1_y
AC_E=AC_E1
AD_E=AD_E1
BC_E=BC_E1
BD_E=BD_E1


Critlayer_4=CriticalLayer(E_x,E_y,AC_E,AD_E,BC_E,BD_E)

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E1 + ": " + str(Critlayer_4) +"nm")











n=E1_x
m=E1_y
av=(Table.loc["av",AB_E0])
ac=(Table.loc["ac",AB_E0])
D0=(Table.loc["D0",AB_E0])
b_=(Table.loc["b_",AB_E0])
Ev_av=(Table.loc["Ev_av",AB_E0])
Eg=(Table.loc["Eg",AB_E0])
LatCon_S_alloy = Quarternary(n,m,(Table.loc["Lattice Constant",AC_E1]),(Table.loc["Lattice Constant",AD_E1]),(Table.loc["Lattice Constant",BC_E1]),(Table.loc["Lattice Constant",BD_E1]))

EV_Epi_E0,EC_Epi_E0= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)
print("Ev for GaSb Cap :  "+ str(EV_Epi_E0))
print("Ec for GaSb Cap:  "+ str(EC_Epi_E0))

x=E1_x
y=E1_y
LatCon_S_alloy = (Table.loc["Lattice Constant","GaSb"])

av=Quarternary(x,y,(Table.loc["av",AC_E1]),(Table.loc["av",AD_E1]),(Table.loc["av",BC_E1]),(Table.loc["av",BD_E1]))
ac=Quarternary(x,y,(Table.loc["ac",AC_E1]),(Table.loc["ac",AD_E1]),(Table.loc["ac",BC_E1]),(Table.loc["ac",BD_E1]))
D0=QuarternaryBowing(x,y,(Table.loc["D0",AC_E1]),(Table.loc["D0",AD_E1]),(Table.loc["D0",BC_E1]),(Table.loc["D0",BD_E1]),C)
b_=Quarternary(x,y,(Table.loc["b_",AC_E1]),(Table.loc["b_",AD_E1]),(Table.loc["b_",BC_E1]),(Table.loc["b_",BD_E1]))
Ev_av=Quarternary(x,y,(Table.loc["Ev_av",AC_E1]),(Table.loc["Ev_av",AD_E1]),(Table.loc["Ev_av",BC_E1]),(Table.loc["Ev_av",BD_E1]))
Eg=QuarternaryBowing(x,y,(Table.loc["Eg",AC_E1]),(Table.loc["Eg",AD_E1]),(Table.loc["Eg",BC_E1]),(Table.loc["Eg",BD_E1]),C)

EV_Epi_E1,EC_Epi_E1= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)
print("Ev for " + ABCD_E1 + ":  "+ str(EV_Epi_E1))
print("Ec for " + ABCD_E1 + ":  "+ str(EC_Epi_E1))


#Barrier
C=Bowing_E2
AC_E=AC_E2
AD_E=AD_E2
BC_E=BC_E2
BD_E=BD_E2
x=E2_x
y=E2_y


EV_Epi_E2,EC_Epi_E2=bandgapallocation(AC_E,AD_E,BC_E,BD_E,x,y,C)

print("Ev for " + ABCD_E2 + ":  "+ str(EV_Epi_E2))
print("Ec for " + ABCD_E2 + ":  "+ str(EC_Epi_E2))


#Well
C=Bowing_E3
AC_E=AC_E3
AD_E=AD_E3
BC_E=BC_E3
BD_E=BD_E3
x=E3_x
y=E3_y

EV_Epi_E3,EC_Epi_E3=bandgapallocation(AC_E,AD_E,BC_E,BD_E,x,y,C)

print("Ev for " + ABCD_E3 + ":  "+ str(EV_Epi_E3))
print("Ec for " + ABCD_E3 + ":  "+ str(EC_Epi_E3))


#Barrier
C=Bowing_E2
AC_E=AC_E2
AD_E=AD_E2
BC_E=BC_E2
BD_E=BD_E2
x=E2_x
y=E2_y

EV_Epi_E22,EC_Epi_E22=bandgapallocation(AC_E,AD_E,BC_E,BD_E,x,y,C)

print("Ev for " + ABCD_E2 + ":  "+ str(EV_Epi_E22))
print("Ec for " + ABCD_E2 + ":  "+ str(EC_Epi_E22))


#Cladding
C=Bowing_E1
AC_E=AC_E1
AD_E=AD_E1
BC_E=BC_E1
BD_E=BD_E1
x=E1_x
y=E1_y


EV_Epi_E12,EC_Epi_E12=bandgapallocation(AC_E,AD_E,BC_E,BD_E,x,y,C)

print("Ev for " + ABCD_E1 + ":  "+ str(EV_Epi_E12))
print("Ec for " + ABCD_E1 + ":  "+ str(EC_Epi_E12))

#Cap
n=E1_x
m=E1_y
av=(Table.loc["av",AB_E0])
ac=(Table.loc["ac",AB_E0])
D0=(Table.loc["D0",AB_E0])
b_=(Table.loc["b_",AB_E0])
Ev_av=(Table.loc["Ev_av",AB_E0])
Eg=(Table.loc["Eg",AB_E0])
LatCon_S_alloy = Quarternary(n,m,(Table.loc["Lattice Constant",AC_E1]),(Table.loc["Lattice Constant",AD_E1]),(Table.loc["Lattice Constant",BC_E1]),(Table.loc["Lattice Constant",BD_E1]))

EV_Epi_E02,EC_Epi_E02= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)
print("Ev for GaSb Cap :  "+ str(EV_Epi_E0))
print("Ec for GaSb Cap:  "+ str(EC_Epi_E0))





EV_Epi_E3=abs(EV_Epi_E3)
EC_Epi_E3=abs(EC_Epi_E3)

BandGap_E3 = EC_Epi_E3 - EV_Epi_E3
print(BandGap_E3)

wavelength = ((6.62607015*10**-34)*(299792458))/(BandGap_E3*(1.6*10**-19))
wavelength = wavelength * 10**6
print("\n" + "Emits a wavelength of : " + str(wavelength) + " micrometers")



##Plotting



##Cap
z=pd.Series(np.arange(0,1000,1))
EVplt_Sub = np.where((z>=0)&(z<=100),EV_Epi_E0,np.nan)
ECplt_Sub = np.where((z>=0)&(z<=100),EC_Epi_E0,np.nan)
plt.vlines(100,EV_Epi_E0,EV_Epi_E1,color="k")
plt.vlines(100,EC_Epi_E0,EC_Epi_E1,color="k")
plt.plot(EVplt_Sub,color="k")
plt.plot(ECplt_Sub,color="k")

#Cladding
EVplt_Epi = np.where((z>=100)&(z<=200),EV_Epi_E1,np.nan)
ECplt_Epi = np.where((z>=100)&(z<=200),EC_Epi_E1,np.nan)
plt.vlines(200,EV_Epi_E1,EV_Epi_E2,color="k")
plt.vlines(200,EC_Epi_E1,EC_Epi_E2,color="k")
plt.plot(EVplt_Epi,color="k")
plt.plot(ECplt_Epi,color="k")


InitialPos=200
WellLength=100/NumberOFWells

for i in range (NumberOFWells-1):
    Position=InitialPos+i*WellLength*2


    #Barrier
    EVplt_Sub = np.where((z>=Position)&(z<=(Position+WellLength)),EV_Epi_E2,np.nan)
    ECplt_Sub = np.where((z>=Position)&(z<=(Position+WellLength)),EC_Epi_E2,np.nan)
    plt.vlines(Position+WellLength,EV_Epi_E2,EV_Epi_E3,color="k")
    plt.vlines(Position+WellLength,EC_Epi_E2,EC_Epi_E3,color="k")
    plt.plot(EVplt_Sub,color="k")
    plt.plot(ECplt_Sub,color="k")

    #Well
    EVplt_Epi = np.where((z>=(Position+WellLength))&(z<=(Position+2*WellLength)),EV_Epi_E3,np.nan)
    ECplt_Epi = np.where((z>=(Position+WellLength))&(z<=(Position+2*WellLength)),EC_Epi_E3,np.nan)
    plt.vlines(Position+2*WellLength,EV_Epi_E3,EV_Epi_E2,color="k")
    plt.vlines(Position+2*WellLength,EC_Epi_E3,EC_Epi_E2,color="k")
    plt.plot(EVplt_Epi,color="k")
    plt.plot(ECplt_Epi,color="k")
    
Position=InitialPos+(i+1)*WellLength*2

#Barrier
EVplt_Sub = np.where((z>=Position)&(z<=(Position+WellLength)),EV_Epi_E2,np.nan)
ECplt_Sub = np.where((z>=Position)&(z<=(Position+WellLength)),EC_Epi_E2,np.nan)
plt.vlines(Position+WellLength,EV_Epi_E2,EV_Epi_E3,color="k")
plt.vlines(Position+WellLength,EC_Epi_E2,EC_Epi_E3,color="k")
plt.plot(EVplt_Sub,color="k")
plt.plot(ECplt_Sub,color="k")

#Well
EVplt_Epi = np.where((z>=(Position+WellLength))&(z<=(Position+2*WellLength)),EV_Epi_E3,np.nan)
ECplt_Epi = np.where((z>=(Position+WellLength))&(z<=(Position+2*WellLength)),EC_Epi_E3,np.nan)
plt.vlines(Position+2*WellLength,EV_Epi_E3,EV_Epi_E22,color="k")
plt.vlines(Position+2*WellLength,EC_Epi_E3,EC_Epi_E22,color="k")
plt.plot(EVplt_Epi,color="k")
plt.plot(ECplt_Epi,color="k")

Position=InitialPos+(i+2)*WellLength*2
#Barrier
EVplt_Sub = np.where((z>=Position)&(z<=(Position+WellLength)),EV_Epi_E22,np.nan)
ECplt_Sub = np.where((z>=Position)&(z<=(Position+WellLength)),EC_Epi_E22,np.nan)
plt.vlines(Position+WellLength,EV_Epi_E22,EV_Epi_E12,color="k")
plt.vlines(Position+WellLength,EC_Epi_E22,EC_Epi_E12,color="k")
plt.plot(EVplt_Sub,color="k")
plt.plot(ECplt_Sub,color="k")



#Cladding
EVplt_Sub = np.where((z>=400+WellLength)&(z<=500+WellLength),EV_Epi_E12,np.nan)
ECplt_Sub = np.where((z>=400+WellLength)&(z<=500+WellLength),EC_Epi_E12,np.nan)
plt.plot(EVplt_Sub,color="k")
plt.plot(ECplt_Sub,color="k")

#Cap
EVplt_Sub = np.where((z>=500+WellLength)&(z<=650),EV_Epi_E02,np.nan)
ECplt_Sub = np.where((z>=500+WellLength)&(z<=650),EC_Epi_E02,np.nan)
plt.vlines(500+WellLength,EV_Epi_E0,EV_Epi_E12,color="k")
plt.vlines(500+WellLength,EC_Epi_E0,EC_Epi_E12,color="k")
plt.plot(EVplt_Sub,color="k")
plt.plot(ECplt_Sub,color="k")

plt.ylabel("Energy")
plt.show()





