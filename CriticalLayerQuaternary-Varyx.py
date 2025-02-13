import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pi=math.pi
bA = 4  
b=bA*(10**-10)





##For InGaAsSb
##InAs-AC
LatCon_AC =    6.0583
C12_AC    =    452.6
C11_AC    =    832.9
##InSb-AD
LatCon_AD =    6.4794
C12_AD    =    373.5
C11_AD    =    684.7
##GaAs-BC
LatCon_BC =    5.65325
C12_BC    =    566.0
C11_BC    =    1221
##GaSb-BD
LatCon_BD =    6.0959
C12_BD    =    402.6
C11_BD    =    884.2


##Substrate-GaSb
LatCon_S =     6.0959
C11_S2    =    884.2
C12_S2    =    402.6


    

y = 0
Level = 20
for j in range (1,Level,1):       ## Calculates a range of x and y values, past 30% values tend and coallese
    y=j/Level
    x = pd.Series(np.arange(0.00,(1-y),0.1))


    #print(x)
    #print(y)

    C11_alloy = x*y*C11_AC+x*(1-y)*C11_AD +(1-x)*y*C11_BC+(1-x)*(1-y)*C11_BD  ##Calculates new average C11 and C12 values
    C12_alloy = x*y*C12_AC+x*(1-y)*C12_AD +(1-x)*y*C12_BC+(1-x)*(1-y)*C12_BD
    LatCon_alloy = x*y*LatCon_AC+x*(1-y)*LatCon_AD +(1-x)*y*LatCon_BC+(1-x)*(1-y)*LatCon_BD


    v=C12_alloy/(C12_alloy+C11_alloy)
    #print("Poission Ratio: \n" + str(v) + "\n")


    A=LatCon_alloy
    a=A*(10**-10)
   # print("Lattice Constant: \n"+ str(A) + "\n")


    f=abs(((A-LatCon_S)/LatCon_S))
   # print("Missfit Percentage: \n"+ str(f*100) + " \n")

    CritThickness = 20
    
    #for j in range (99):
        #CritThicknesstemp=1+CritThicknesstemp
        #ritThickness=CritThicknesstemp
        
    for i in range (99):
            CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))

    CritThickness= CritThickness*10**9
          

    print(str(CritThickness) +"nm")
    plt.plot(x,CritThickness, label = str(y))

#plt.title(str(y))
plt.yscale("log")
plt.ylim(1,200)
plt.ylabel("Critical thickness (nm)")
plt.xlabel("Percentage quantity of x")
plt.legend()
plt.show()










