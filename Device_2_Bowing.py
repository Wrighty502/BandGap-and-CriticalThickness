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
ABCD_E2 = "AlInGaAsSb-Barrier"
E2_x=0.2
E2_y=0.25
E2_z=0.24


AP_E2 = "AlAs"
BP_E2 = "InAs"
CP_E2 = "GaAs"
AQ_E2 = "AlSb"
BQ_E2 = "InSb"
CQ_E2 = "GaSb"
Bowing_E2 = 0.75

##QW
ABCD_E3 = "InGaAsSb-Well"
E3_x=0.54
E3_y=0.24


AC_E3 = "InAs"
AD_E3 = "InSb"
BC_E3 = "GaAs"
BD_E3 = "GaSb"
Bowing_E3 = 0.75

AB_E0 = "GaSb"

NumberOFWells=8




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



def Quarternary(x,y,AC,AD,BC,BD):
    QRT_ans = x*y*AC+x*(1-y)*AD +(1-x)*y*BC+(1-x)*(1-y)*BD
    return QRT_ans

def QuarternaryBowing(x,y,AC,AD,BC,BD,C):
    QRT_ans = x*y*AC+x*(1-y)*AD +(1-x)*y*BC+(1-x)*(1-y)*BD-x*y*(1-y)*C-(1-x)*y*(1-y)*C-x*(1-x)*y*C-x*(1-x)*(1-y)*C
    return QRT_ans

def Quinary(x,y,z,AP,BP,CP,AQ,BQ,CQ):
    QUIN_ans = x*z*AP+y*z*BP+(1-x-y)*z*CP+x*(1-z)*AQ+y*(1-z)*BQ+(1-x-y)*(1-z)*CQ
    return QUIN_ans

def QuinaryBowing(x,y,z,AP,BP,CP,AQ,BQ,CQ):
    QUIN_ans = x*z*AP+y*z*BP+(1-x-y)*z*CP+x*(1-z)*AQ+y*(1-z)*BQ+(1-x-y)*(1-z)*CQ
    return QUIN_ans







                ##([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",  "b_",   "Ev-av",   "Eg"]) ## Just to make adding data easier with labels
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


















##Epilayer-Barrier
x=E2_x
y=E2_y
z=E2_z


C11_alloy = Quinary(x,y,z,(Table.loc["C11",AP_E2]),(Table.loc["C11",BP_E2]),(Table.loc["C11",CP_E2]),(Table.loc["C11",AQ_E2]),(Table.loc["C11",BQ_E2]),(Table.loc["C11",CQ_E2]))
C12_alloy = Quinary(x,y,z,(Table.loc["C12",AP_E2]),(Table.loc["C12",BP_E2]),(Table.loc["C12",CP_E2]),(Table.loc["C12",AQ_E2]),(Table.loc["C12",BQ_E2]),(Table.loc["C12",CQ_E2]))
LatCon_alloy = Quinary(x,y,z,(Table.loc["Lattice Constant",AP_E2]),(Table.loc["Lattice Constant",BP_E2]),(Table.loc["Lattice Constant",CP_E2]),(Table.loc["Lattice Constant",AQ_E2]),(Table.loc["Lattice Constant",BQ_E2]),(Table.loc["Lattice Constant",CQ_E2]))


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
CritThickness_E1_E2 = CritThickness*10**9

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E2 + ": " + str(CritThickness_E1_E2) +"nm")








##Epilayer-Well
x=E3_x
y=E3_y           

C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E3]),(Table.loc["C11",AD_E3]),(Table.loc["C11",BC_E3]),(Table.loc["C11",BD_E3]))
C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E3]),(Table.loc["C12",AD_E3]),(Table.loc["C12",BC_E3]),(Table.loc["C12",BD_E3]))
LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E3]),(Table.loc["Lattice Constant",AD_E3]),(Table.loc["Lattice Constant",BC_E3]),(Table.loc["Lattice Constant",BD_E3]))


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


for i in range (6):  ##More then 10 or so means Crit-Thickness is negitive for any initial
    CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
    #print(CritThickness)
CritThickness_E2_E3 = CritThickness*10**9

print("Critical Thickness " + AB_E0 +" - "+ ABCD_E3 + ": " + str(CritThickness_E2_E3) +"nm")











##Epilayer-Barrier
x=E2_x
y=E2_y
z=E2_z

C11_alloy = Quinary(x,y,z,(Table.loc["C11",AP_E2]),(Table.loc["C11",BP_E2]),(Table.loc["C11",CP_E2]),(Table.loc["C11",AQ_E2]),(Table.loc["C11",BQ_E2]),(Table.loc["C11",CQ_E2]))
C12_alloy = Quinary(x,y,z,(Table.loc["C12",AP_E2]),(Table.loc["C12",BP_E2]),(Table.loc["C12",CP_E2]),(Table.loc["C12",AQ_E2]),(Table.loc["C12",BQ_E2]),(Table.loc["C12",CQ_E2]))
LatCon_alloy = Quinary(x,y,z,(Table.loc["Lattice Constant",AP_E2]),(Table.loc["Lattice Constant",BP_E2]),(Table.loc["Lattice Constant",CP_E2]),(Table.loc["Lattice Constant",AQ_E2]),(Table.loc["Lattice Constant",BQ_E2]),(Table.loc["Lattice Constant",CQ_E2]))


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
CritThickness_E2_E3 = CritThickness*10**9

print("2nd Critical Thickness " + AB_E0 +" - "+ ABCD_E1 + ": " + str(CritThickness_E2_E3) +"nm")












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
C=Bowing_E1
C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E1]),(Table.loc["C11",AD_E1]),(Table.loc["C11",BC_E1]),(Table.loc["C11",BD_E1]))
C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E1]),(Table.loc["C12",AD_E1]),(Table.loc["C12",BC_E1]),(Table.loc["C12",BD_E1]))
LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E1]),(Table.loc["Lattice Constant",AD_E1]),(Table.loc["Lattice Constant",BC_E1]),(Table.loc["Lattice Constant",BD_E1]))


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
x=E2_x
y=E2_y
z=E2_z

C11_alloy = Quinary(x,y,z,(Table.loc["C11",AP_E2]),(Table.loc["C11",BP_E2]),(Table.loc["C11",CP_E2]),(Table.loc["C11",AQ_E2]),(Table.loc["C11",BQ_E2]),(Table.loc["C11",CQ_E2]))
C12_alloy = Quinary(x,y,z,(Table.loc["C12",AP_E2]),(Table.loc["C12",BP_E2]),(Table.loc["C12",CP_E2]),(Table.loc["C12",AQ_E2]),(Table.loc["C12",BQ_E2]),(Table.loc["C12",CQ_E2]))
LatCon_alloy = Quinary(x,y,z,(Table.loc["Lattice Constant",AP_E2]),(Table.loc["Lattice Constant",BP_E2]),(Table.loc["Lattice Constant",CP_E2]),(Table.loc["Lattice Constant",AQ_E2]),(Table.loc["Lattice Constant",BQ_E2]),(Table.loc["Lattice Constant",CQ_E2]))


LatCon_S_alloy = (Table.loc["Lattice Constant","GaSb"])
av=Quinary(x,y,z,(Table.loc["av",AP_E2]),(Table.loc["av",BP_E2]),(Table.loc["av",CP_E2]),(Table.loc["av",AQ_E2]),(Table.loc["av",BQ_E2]),(Table.loc["av",CQ_E2]))
ac=Quinary(x,y,z,(Table.loc["ac",AP_E2]),(Table.loc["ac",BP_E2]),(Table.loc["ac",CP_E2]),(Table.loc["ac",AQ_E2]),(Table.loc["ac",BQ_E2]),(Table.loc["ac",CQ_E2]))
D0=Quinary(x,y,z,(Table.loc["D0",AP_E2]),(Table.loc["D0",BP_E2]),(Table.loc["D0",CP_E2]),(Table.loc["D0",AQ_E2]),(Table.loc["D0",BQ_E2]),(Table.loc["D0",CQ_E2]))
b_=Quinary(x,y,z,(Table.loc["b_",AP_E2]),(Table.loc["b_",BP_E2]),(Table.loc["b_",CP_E2]),(Table.loc["b_",AQ_E2]),(Table.loc["b_",BQ_E2]),(Table.loc["b_",CQ_E2]))
Ev_av=Quinary(x,y,z,(Table.loc["Ev_av",AP_E2]),(Table.loc["Ev_av",BP_E2]),(Table.loc["Ev_av",CP_E2]),(Table.loc["Ev_av",AQ_E2]),(Table.loc["Ev_av",BQ_E2]),(Table.loc["Ev_av",CQ_E2]))
Eg=Quinary(x,y,z,(Table.loc["Eg",AP_E2]),(Table.loc["Eg",BP_E2]),(Table.loc["Eg",CP_E2]),(Table.loc["Eg",AQ_E2]),(Table.loc["Eg",BQ_E2]),(Table.loc["Eg",CQ_E2]))



EV_Epi_E2,EC_Epi_E2= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)


print("Ev for " + ABCD_E2 + ":  "+ str(EV_Epi_E2))
print("Ec for " + ABCD_E2 + ":  "+ str(EC_Epi_E2))


#Well
C=Bowing_E3
x=E3_x
y=E3_y
C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E3]),(Table.loc["C11",AD_E3]),(Table.loc["C11",BC_E3]),(Table.loc["C11",BD_E3]))
C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E3]),(Table.loc["C12",AD_E3]),(Table.loc["C12",BC_E3]),(Table.loc["C12",BD_E3]))
LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E3]),(Table.loc["Lattice Constant",AD_E3]),(Table.loc["Lattice Constant",BC_E3]),(Table.loc["Lattice Constant",BD_E3]))


LatCon_S_alloy = (Table.loc["Lattice Constant","GaSb"])
av=Quarternary(x,y,(Table.loc["av",AC_E3]),(Table.loc["av",AD_E3]),(Table.loc["av",BC_E3]),(Table.loc["av",BD_E3]))
ac=Quarternary(x,y,(Table.loc["ac",AC_E3]),(Table.loc["ac",AD_E3]),(Table.loc["ac",BC_E3]),(Table.loc["ac",BD_E3]))
D0=QuarternaryBowing(x,y,(Table.loc["D0",AC_E3]),(Table.loc["D0",AD_E3]),(Table.loc["D0",BC_E3]),(Table.loc["D0",BD_E3]),C)
b_=Quarternary(x,y,(Table.loc["b_",AC_E3]),(Table.loc["b_",AD_E3]),(Table.loc["b_",BC_E3]),(Table.loc["b_",BD_E3]))
Ev_av=Quarternary(x,y,(Table.loc["Ev_av",AC_E3]),(Table.loc["Ev_av",AD_E3]),(Table.loc["Ev_av",BC_E3]),(Table.loc["Ev_av",BD_E3]))
Eg=QuarternaryBowing(x,y,(Table.loc["Eg",AC_E3]),(Table.loc["Eg",AD_E3]),(Table.loc["Eg",BC_E3]),(Table.loc["Eg",BD_E3]),C)
EV_Epi_E3,EC_Epi_E3= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)

print("Ev for " + ABCD_E3 + ":  "+ str(EV_Epi_E3))
print("Ec for " + ABCD_E3 + ":  "+ str(EC_Epi_E3))


#Barrier
x=E2_x
y=E2_y
z=E2_z
LatCon_S_alloy = (Table.loc["Lattice Constant","GaSb"])

C11_alloy = Quinary(x,y,z,(Table.loc["C11",AP_E2]),(Table.loc["C11",BP_E2]),(Table.loc["C11",CP_E2]),(Table.loc["C11",AQ_E2]),(Table.loc["C11",BQ_E2]),(Table.loc["C11",CQ_E2]))
C12_alloy = Quinary(x,y,z,(Table.loc["C12",AP_E2]),(Table.loc["C12",BP_E2]),(Table.loc["C12",CP_E2]),(Table.loc["C12",AQ_E2]),(Table.loc["C12",BQ_E2]),(Table.loc["C12",CQ_E2]))
LatCon_alloy = Quinary(x,y,z,(Table.loc["Lattice Constant",AP_E2]),(Table.loc["Lattice Constant",BP_E2]),(Table.loc["Lattice Constant",CP_E2]),(Table.loc["Lattice Constant",AQ_E2]),(Table.loc["Lattice Constant",BQ_E2]),(Table.loc["Lattice Constant",CQ_E2]))

av=Quinary(x,y,z,(Table.loc["av",AP_E2]),(Table.loc["av",BP_E2]),(Table.loc["av",CP_E2]),(Table.loc["av",AQ_E2]),(Table.loc["av",BQ_E2]),(Table.loc["av",CQ_E2]))
ac=Quinary(x,y,z,(Table.loc["ac",AP_E2]),(Table.loc["ac",BP_E2]),(Table.loc["ac",CP_E2]),(Table.loc["ac",AQ_E2]),(Table.loc["ac",BQ_E2]),(Table.loc["ac",CQ_E2]))
D0=Quinary(x,y,z,(Table.loc["D0",AP_E2]),(Table.loc["D0",BP_E2]),(Table.loc["D0",CP_E2]),(Table.loc["D0",AQ_E2]),(Table.loc["D0",BQ_E2]),(Table.loc["D0",CQ_E2]))
b_=Quinary(x,y,z,(Table.loc["b_",AP_E2]),(Table.loc["b_",BP_E2]),(Table.loc["b_",CP_E2]),(Table.loc["b_",AQ_E2]),(Table.loc["b_",BQ_E2]),(Table.loc["b_",CQ_E2]))
Ev_av=Quinary(x,y,z,(Table.loc["Ev_av",AP_E2]),(Table.loc["Ev_av",BP_E2]),(Table.loc["Ev_av",CP_E2]),(Table.loc["Ev_av",AQ_E2]),(Table.loc["Ev_av",BQ_E2]),(Table.loc["Ev_av",CQ_E2]))
Eg=Quinary(x,y,z,(Table.loc["Eg",AP_E2]),(Table.loc["Eg",BP_E2]),(Table.loc["Eg",CP_E2]),(Table.loc["Eg",AQ_E2]),(Table.loc["Eg",BQ_E2]),(Table.loc["Eg",CQ_E2]))

EV_Epi_E22,EC_Epi_E22= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)

print("Ev for " + ABCD_E2 + ":  "+ str(EV_Epi_E22))
print("Ec for " + ABCD_E2 + ":  "+ str(EC_Epi_E22))


#Cladding
C=Bowing_E1
x=E1_x
y=E1_y
LatCon_S_alloy = (Table.loc["Lattice Constant","GaSb"])
C11_alloy = Quarternary(x,y,(Table.loc["C11",AC_E1]),(Table.loc["C11",AD_E1]),(Table.loc["C11",BC_E1]),(Table.loc["C11",BD_E1]))
C12_alloy = Quarternary(x,y,(Table.loc["C12",AC_E1]),(Table.loc["C12",AD_E1]),(Table.loc["C12",BC_E1]),(Table.loc["C12",BD_E1]))
LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",AC_E1]),(Table.loc["Lattice Constant",AD_E1]),(Table.loc["Lattice Constant",BC_E1]),(Table.loc["Lattice Constant",BD_E1]))

av=Quarternary(x,y,(Table.loc["av",AC_E1]),(Table.loc["av",AD_E1]),(Table.loc["av",BC_E1]),(Table.loc["av",BD_E1]))
ac=Quarternary(x,y,(Table.loc["ac",AC_E1]),(Table.loc["ac",AD_E1]),(Table.loc["ac",BC_E1]),(Table.loc["ac",BD_E1]))
D0=QuarternaryBowing(x,y,(Table.loc["D0",AC_E1]),(Table.loc["D0",AD_E1]),(Table.loc["D0",BC_E1]),(Table.loc["D0",BD_E1]),C)
b_=Quarternary(x,y,(Table.loc["b_",AC_E1]),(Table.loc["b_",AD_E1]),(Table.loc["b_",BC_E1]),(Table.loc["b_",BD_E1]))
Ev_av=Quarternary(x,y,(Table.loc["Ev_av",AC_E1]),(Table.loc["Ev_av",AD_E1]),(Table.loc["Ev_av",BC_E1]),(Table.loc["Ev_av",BD_E1]))
Eg=QuarternaryBowing(x,y,(Table.loc["Eg",AC_E1]),(Table.loc["Eg",AD_E1]),(Table.loc["Eg",BC_E1]),(Table.loc["Eg",BD_E1]),C)
EV_Epi_E12,EC_Epi_E12= bandgap(LatCon_S_alloy,LatCon_alloy,C11_alloy,C12_alloy,av,ac,D0,b_,Ev_av,Eg)


print("Ev for " + ABCD_E1 + ":  "+ str(EV_Epi_E12))
print("Ec for " + ABCD_E1 + ":  "+ str(EC_Epi_E12))

#Cap
n=E1_x
m=E1_y
C12_alloy = (Table.loc["C12","GaSb"])
C11_alloy = (Table.loc["C11","GaSb"])
LatCon_alloy = (Table.loc["Lattice Constant","GaSb"])
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





