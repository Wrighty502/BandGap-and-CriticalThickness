import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option("display.max.columns",10)

pi=math.pi

bA = 4          #Slip distnace (Magnititude of Burger's Vector) b in angstrograms
b=bA*(10**-10)




def bandgap(C_11,C_12,av,ac,D0,b,EV_av,Eg): ## Calculates band gap using Krijn model

    
    a0=LatCon_S1
    a=LatCon_C1
    aPara=a0
    D001=2*(C_12/C_11)
    aPerp=a*(1-D001*((aPara/a)-1))
    ePerp=(aPerp/a)-1
    ePara=(aPara/a)-1

    EHy_v=av*(2*ePara+ePerp)
    EHy_c=ac*(2*ePara+ePerp)
    
    dE=2*b*(ePerp-ePara)
    ESH_hh=-0.5*dE
    ESH_lh=(-0.5*D0)+(0.25*dE)+(0.5*(D0**2+D0*dE+(9/4)*dE**2)**0.5)


    EV=EV_av+D0/3+EHy_v+max(ESH_hh,ESH_lh)
    EC=EV_av+D0/3+Eg+EHy_c
    return EV, EC

def Ternary(x,AB,AC,C):
    Trny_ans = x*AB+(1-x)*AC+x*(x-1)*C
    return Trny_ans

def Quarternary(x,y,AB,AC,AD):
    QRT_ans = x*AB+y*AC+(1-x-y)*AD
    return QRT_ans





                ##([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",   "b",   "Ev-av",   "Eg"]) ## Just to make adding data easier with labels
InAsTable = pd.Series([     6.0583,         832.9,  452.6,   1.00,  -5.08,   0.39,  -1.8,     -0.59,  0.417])
GaAsTable = pd.Series([     5.6530,          1221,  566.0,   1.16,  -7.17,   0.34,  -2.0,     -0.80,  1.420])
InSbTable = pd.Series([     6.4794,         684.7,  373.5,   0.36,  -6.94,   0.81,  -2.0,     -0.00,  0.235])
InPTable  = pd.Series([     5.8697,          1011,  561.0,   0.60,  -6.00,   0.108, -2.0,     -0.94,  1.4236])
AlASTable = pd.Series([     5.6611,          1250,  534.0,   2.47,  -5.64,   0.28,  -2.3,     -1.33,  3.099])
GaPTable  = pd.Series([     5.4505,          1405,  620.3,   1.70,  -8.20,   0.08,  -1.6,     -1.27,  2.886])
GaSbTable = pd.Series([     6.0959,         884.2,  402.6,    0.8,  -7.50,   0.76,  -2.0,     -0.03,  0.812])
AlPTable  = pd.Series([     5.4672,          1330,    630,    3.0,  -5.70,   0.07,  -1.5,     -1.74,  3.630])
AlSbTable = pd.Series([     6.1355,         876.9,  434.1,    1.4,  -4.50,  0.676, -1.35,     -0.41,  2.386])

Table=pd.DataFrame({"InAs":InAsTable,"GaAs":GaAsTable,"InSb":InSbTable, "InP":InPTable,"AlAs":AlASTable,"GaP":GaPTable,"GaSb":GaSbTable,"AlP":AlPTable,"AlSb":AlSbTable})
Table.index =([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",   "b",   "Ev_av",   "Eg"])
print("\n")
print(Table)
print("\n")

#print(Table.loc["C11","InAs"])


Type = input("How many elements are in the epilayer? (i.e. Binary(2), Ternary(3), Quarternary(4)")
if Type == "2":
    EpiSelection = input("What compound is epilayer?")
    print("\n")
    #print(Table[str(EpiSelection)])
    #print("\n")
    print("\n")
    SubSelection = input("What compound is substrate?")
    print("\n")
    #print(Table[str(SubSelection)])
    #print("\n")
    print("\n")


    
    v=(Table.loc["C12",EpiSelection])/((Table.loc["C12",EpiSelection])+(Table.loc["C11",EpiSelection]))
    print("Poission Ratio: " + str(v) + "\n")
    
    A=(Table.loc["Lattice Constant",EpiSelection])
    a=A*(10**-10)
    print("Lattice Constant: "+ str(A) + "\n")
    
    f=abs((A-(Table.loc["Lattice Constant",SubSelection]))/(Table.loc["Lattice Constant",SubSelection]))
    print("Missfit Percentage: "+ str(f*100) + "% \n")

    CritThickness=2


    for i in range (99):
        if CritThickness > 0:
            CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
            #CritThickness = abs(CritThickness)
        else:
            print("Critical thickness became negitive.")
            break
    CritThickness= CritThickness*10**9
    print(str(CritThickness) +"nm")
    
    Is_BG_wanted = input("Do you want a band gap diagram? (Y/N)")
    if Is_BG_wanted == "Y":

        LatCon_C1=(Table.loc["Lattice Constant",EpiSelection])
        LatCon_S1=(Table.loc["Lattice Constant",SubSelection])
        C11_C1=(Table.loc["C11",EpiSelection])
        C12_C1=(Table.loc["C12",EpiSelection])
        av_C1=(Table.loc["av",EpiSelection])
        ac_C1=(Table.loc["ac",EpiSelection])
        D0_C1=(Table.loc["D0",EpiSelection])
        b_C1=(Table.loc["b",EpiSelection])
        EV_av_C1=(Table.loc["Ev_av",EpiSelection])
        Eg_C1=(Table.loc["Eg",EpiSelection])


        EV_Epi,EC_Epi= bandgap(C11_C1,C12_C1,av_C1,ac_C1,D0_C1,b_C1,EV_av_C1,Eg_C1)
        print("Ev for "+ EpiSelection +":"+ str(EV_Epi))
        print("Ec for "+ EpiSelection +":"+ str(EC_Epi))

        C11_S1=(Table.loc["C11",SubSelection])
        C12_S1=(Table.loc["C12",SubSelection])
        av_S1=(Table.loc["av",SubSelection])
        ac_S1=(Table.loc["ac",SubSelection])
        D0_S1=(Table.loc["D0",SubSelection])
        b_S1=(Table.loc["b",SubSelection])
        EV_av_S1=(Table.loc["Ev_av",SubSelection])
        Eg_S1=(Table.loc["Eg",SubSelection])


        EV_Sub,EC_Sub= bandgap(C11_S1,C12_S1,av_S1,ac_S1,D0_S1,b_S1,EV_av_S1,Eg_S1)
        print("Ev for "+ SubSelection +":"+ str(EV_Sub))
        print("Ec for "+ SubSelection +":"+ str(EC_Sub))

        if EV_Sub<EV_Epi and EC_Sub>EC_Epi:

            Is_QW_wanted = input("Do you want a Quantum Well? (Y/N)")

        
            if Is_QW_wanted == "Y":
                z=pd.Series(np.arange(-100,100,1))
                EVplt_Sub = np.where((z>=-50)&(z<=0),EV_Sub,np.nan)
                ECplt_Sub = np.where((z>=-50)&(z<=0),EC_Sub,np.nan)
                plt.vlines(100,EV_Epi,EV_Sub,color="k")
                plt.vlines(100,EC_Epi,EC_Sub,color="k")
                plt.plot(EVplt_Sub,color="k")
                plt.plot(ECplt_Sub,color="k")
                EVplt_Epi = np.where((z<=50)&(z>=0),EV_Epi,np.nan)
                ECplt_Epi = np.where((z<=50)&(z>=0),EC_Epi,np.nan)
                plt.plot(EVplt_Epi,color="k")
                plt.plot(ECplt_Epi,color="k")
                EVplt_Sub = np.where(z>=50,EV_Sub,np.nan)
                ECplt_Sub = np.where(z>=50,EC_Sub,np.nan)
                plt.vlines(150,EV_Epi,EV_Sub,color="k")
                plt.vlines(150,EC_Epi,EC_Sub,color="k")
                plt.plot(EVplt_Sub,color="k")
                plt.plot(ECplt_Sub,color="k")
                EVplt_Sub = np.where((z>=-50)&(z<=0),EV_Sub,np.nan)
                ECplt_Sub = np.where((z>=-50)&(z<=0),EC_Sub,np.nan)
                plt.plot(EVplt_Sub,color="k")
                plt.plot(ECplt_Sub,color="k")

                
                #plt.title()
                #plt.ylim()
                plt.ylabel("Energy")
                plt.xlabel(SubSelection + "   -----   " +EpiSelection + "   -----   " + SubSelection)
                #plt.legend()
                plt.show()


                
        else:
            z=pd.Series(np.arange(0,100,1))
            EVplt_Epi = np.where(z<=50,EV_Epi,np.nan)
            ECplt_Epi = np.where(z<=50,EC_Epi,np.nan)
            plt.vlines(50,EV_Epi,EV_Sub,color="k")
            plt.vlines(50,EC_Epi,EC_Sub,color="k")
            plt.plot(EVplt_Epi,color="k")
            plt.plot(ECplt_Epi,color="k")
            EVplt_Sub = np.where(z>=50,EV_Sub,np.nan)
            ECplt_Sub = np.where(z>=50,EC_Sub,np.nan)
            plt.plot(EVplt_Sub,color="k")
            plt.plot(ECplt_Sub,color="k")
            
            #plt.title()
            #plt.ylim()
            plt.ylabel("Energy")
            plt.xlabel(EpiSelection+ "   -----   " + SubSelection)
            #plt.legend()
            plt.show()


            
if Type == "3":
    Varyx = input(" Would you like to vary x to show how critical thickness changes? (Y/N)")
    if  Varyx == "N":
        EpiSelection = input("What compound for the epilayer? ABC")
        EpiSelectionAB = input("What compound for the epilayer AB, where the alloy is ABC? (e.g wanting InPAs, AB is InP and AC is InAs")
        EpiSelectionAC = input("What compound for the epilayer AC?")
        BowingParamE = input("What is the bowing parameter for C(Eg)? (input 0 if unsure)")
        BowingParamD = input("What is the bowing parameter for C(D0)? (input 0 if unsure)")
        x = input("What is x as a decimal")
        x=float(x)
        print("\n")
        SubSelection = input("What compound is substrate?")
        print("\n")
        print("\n")
        
        LatCon_alloy = Ternary(x,(Table.loc["Lattice Constant",EpiSelectionAB]),(Table.loc["Lattice Constant",EpiSelectionAC]),0.0)
        C11_alloy = Ternary(x,Table.loc["C11",EpiSelectionAB],Table.loc["C11",EpiSelectionAC],0)
        C12_alloy = Ternary(x,Table.loc["C12",EpiSelectionAB],Table.loc["C12",EpiSelectionAC],0)
        av_alloy = Ternary(x,Table.loc["av",EpiSelectionAB],Table.loc["av",EpiSelectionAC],0)
        ac_alloy = Ternary(x,Table.loc["ac",EpiSelectionAB],Table.loc["ac",EpiSelectionAC],0)
        D0_alloy = Ternary(x,Table.loc["D0",EpiSelectionAB],Table.loc["D0",EpiSelectionAC],float(BowingParamD))
        b_alloy = Ternary(x,Table.loc["b",EpiSelectionAB],Table.loc["b",EpiSelectionAC],0)
        Cv_av = 3*(Table.loc["av",EpiSelectionAB]-Table.loc["av",EpiSelectionAC])*((Table.loc["Lattice Constant",EpiSelectionAB]-Table.loc["Lattice Constant",EpiSelectionAC])/Table.loc["Lattice Constant",SubSelection])
        Ev_av_alloy = Ternary(x,Table.loc["Ev_av",EpiSelectionAB],Table.loc["Ev_av",EpiSelectionAC],Cv_av)
        Eg_alloy = Ternary(x,Table.loc["Eg",EpiSelectionAB],Table.loc["Eg",EpiSelectionAC],float(BowingParamE))

                           
        TrnyTabledata = pd.Series([LatCon_alloy,C11_alloy,C12_alloy,av_alloy,ac_alloy,D0_alloy,b_alloy,Ev_av_alloy,Eg_alloy])
        TrnyTable=pd.DataFrame({str(EpiSelection):TrnyTabledata})
        TrnyTable.index =([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",   "b",   "Ev_av",   "Eg"])
        print(TrnyTable)
        
        v=(TrnyTable.loc["C12",EpiSelection])/((TrnyTable.loc["C12",EpiSelection])+(TrnyTable.loc["C11",EpiSelection]))
        print("Poission Ratio: " + str(v) + "\n")
        
        A=(TrnyTable.loc["Lattice Constant",EpiSelection])
        a=A*(10**-10)
        print("Lattice Constant: "+ str(A) + "\n")
        
        f=abs((A-(Table.loc["Lattice Constant",SubSelection]))/(Table.loc["Lattice Constant",SubSelection]))
        print("Missfit Percentage: "+ str(f*100) + "% \n")

        CritThickness=2
        
        for i in range (99):
            if CritThickness > 0:
                CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
                #CritThickness = abs(CritThickness)
            else:
                print("Critical thickness became negitive.")
                break
        CritThickness= CritThickness*10**9
        print(str(CritThickness) +"nm")

        
        Is_BG_wanted = input("Do you want a band gap diagram? (Y/N)")
        if Is_BG_wanted == "Y":

            LatCon_C1=(TrnyTable.loc["Lattice Constant",EpiSelection])
            LatCon_S1=(Table.loc["Lattice Constant",SubSelection])
            C11_C1=(TrnyTable.loc["C11",EpiSelection])
            C12_C1=(TrnyTable.loc["C12",EpiSelection])
            av_C1=(TrnyTable.loc["av",EpiSelection])
            ac_C1=(TrnyTable.loc["ac",EpiSelection])
            D0_C1=(TrnyTable.loc["D0",EpiSelection])
            b_C1=(TrnyTable.loc["b",EpiSelection])
            EV_av_C1=(TrnyTable.loc["Ev_av",EpiSelection])
            Eg_C1=(TrnyTable.loc["Eg",EpiSelection])


            EV_Epi,EC_Epi= bandgap(C11_C1,C12_C1,av_C1,ac_C1,D0_C1,b_C1,EV_av_C1,Eg_C1)
            print("Ev for "+ EpiSelection +":"+ str(EV_Epi))
            print("Ec for "+ EpiSelection +":"+ str(EC_Epi))

            C11_S1=(Table.loc["C11",SubSelection])
            C12_S1=(Table.loc["C12",SubSelection])
            av_S1=(Table.loc["av",SubSelection])
            ac_S1=(Table.loc["ac",SubSelection])
            D0_S1=(Table.loc["D0",SubSelection])
            b_S1=(Table.loc["b",SubSelection])
            EV_av_S1=(Table.loc["Ev_av",SubSelection])
            Eg_S1=(Table.loc["Eg",SubSelection])


            EV_Sub,EC_Sub= bandgap(C11_S1,C12_S1,av_S1,ac_S1,D0_S1,b_S1,EV_av_S1,Eg_S1)
            print("Ev for "+ SubSelection +":"+ str(EV_Sub))
            print("Ec for "+ SubSelection +":"+ str(EC_Sub))

            if EV_Sub<EV_Epi and EC_Sub>EC_Epi:

                Is_QW_wanted = input("Do you want a Quantum Well? (Y/N)")

            
                if Is_QW_wanted == "Y":
                    z=pd.Series(np.arange(-100,100,1))
                    EVplt_Sub = np.where((z>=-50)&(z<=0),EV_Sub,np.nan)
                    ECplt_Sub = np.where((z>=-50)&(z<=0),EC_Sub,np.nan)
                    plt.vlines(100,EV_Epi,EV_Sub,color="k")
                    plt.vlines(100,EC_Epi,EC_Sub,color="k")
                    plt.plot(EVplt_Sub,color="k")
                    plt.plot(ECplt_Sub,color="k")
                    EVplt_Epi = np.where((z<=50)&(z>=0),EV_Epi,np.nan)
                    ECplt_Epi = np.where((z<=50)&(z>=0),EC_Epi,np.nan)
                    plt.plot(EVplt_Epi,color="k")
                    plt.plot(ECplt_Epi,color="k")
                    EVplt_Sub = np.where(z>=50,EV_Sub,np.nan)
                    ECplt_Sub = np.where(z>=50,EC_Sub,np.nan)
                    plt.vlines(150,EV_Epi,EV_Sub,color="k")
                    plt.vlines(150,EC_Epi,EC_Sub,color="k")
                    plt.plot(EVplt_Sub,color="k")
                    plt.plot(ECplt_Sub,color="k")
                    EVplt_Sub = np.where((z>=-50)&(z<=0),EV_Sub,np.nan)
                    ECplt_Sub = np.where((z>=-50)&(z<=0),EC_Sub,np.nan)
                    plt.plot(EVplt_Sub,color="k")
                    plt.plot(ECplt_Sub,color="k")

                    
                    #plt.title()
                    #plt.ylim()
                    plt.ylabel("Energy")
                    plt.xlabel(SubSelection + "   -----   " +EpiSelection + "   -----   " + SubSelection)
                    #plt.legend()
                    plt.show()


                
            else:
                z=pd.Series(np.arange(0,100,1))
                EVplt_Epi = np.where(z<=50,EV_Epi,np.nan)
                ECplt_Epi = np.where(z<=50,EC_Epi,np.nan)
                plt.vlines(50,EV_Epi,EV_Sub,color="k")
                plt.vlines(50,EC_Epi,EC_Sub,color="k")
                plt.plot(EVplt_Epi,color="k")
                plt.plot(ECplt_Epi,color="k")
                EVplt_Sub = np.where(z>=50,EV_Sub,np.nan)
                ECplt_Sub = np.where(z>=50,EC_Sub,np.nan)
                plt.plot(EVplt_Sub,color="k")
                plt.plot(ECplt_Sub,color="k")
                
                #plt.title()
                #plt.ylim()
                plt.ylabel("Energy")
                plt.xlabel(EpiSelection+ "   -----   " + SubSelection)
                #plt.legend()
                plt.show()
                
    if Varyx=="Y":
        EpiSelection = input("What compound for the epilayer? ABC")
        EpiSelectionAB = input("What compound for the epilayer AB, where the alloy is ABC? (e.g wanting InPAs, AB is InP and AC is InAs")
        EpiSelectionAC = input("What compound for the epilayer AC?")
        BowingParamE = input("What is the bowing parameter for C(Eg)? (input 0 if unsure)")
        BowingParamD = input("What is the bowing parameter for C(D0)? (input 0 if unsure)")
        x=pd.Series(np.arange(0,1,0.01))
        #print(x)
        SubSelection = input("What compound is substrate?")

        LatCon_alloy = Ternary(x,(Table.loc["Lattice Constant",EpiSelectionAB]),(Table.loc["Lattice Constant",EpiSelectionAC]),0.0)
        C11_alloy = Ternary(x,Table.loc["C11",EpiSelectionAB],Table.loc["C11",EpiSelectionAC],0)
        C12_alloy = Ternary(x,Table.loc["C12",EpiSelectionAB],Table.loc["C12",EpiSelectionAC],0)

        v=(C12_alloy)/((C12_alloy)+(C11_alloy))
        print("Poission Ratio: " + str(v) + "\n")
        
        A=(LatCon_alloy)
        a=A*(10**-10)
        print("Lattice Constant: "+ str(A) + "\n")
        
        f=abs((A-(Table.loc["Lattice Constant",SubSelection]))/(Table.loc["Lattice Constant",SubSelection]))
        print("Missfit Percentage: "+ str(f*100) + "% \n")

        CritThickness=2
        
        for i in range (99):
                CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
                #CritThickness = abs(CritThickness)
        CritThickness= CritThickness*10**9
        print(str(CritThickness) +"nm")
        plt.plot(x,CritThickness)
        plt.yscale("log")
        plt.ylabel("Critical thickness")
        plt.xlabel("Percentage quantity of x")
        x=x[CritThickness.idxmax()]

        print ("Max Crit, x at:" + str(x))
        plt.show()
        

        LatCon_alloy = Ternary(x,(Table.loc["Lattice Constant",EpiSelectionAB]),(Table.loc["Lattice Constant",EpiSelectionAC]),0.0)
        C11_alloy = Ternary(x,Table.loc["C11",EpiSelectionAB],Table.loc["C11",EpiSelectionAC],0)
        C12_alloy = Ternary(x,Table.loc["C12",EpiSelectionAB],Table.loc["C12",EpiSelectionAC],0)
        av_alloy = Ternary(x,Table.loc["av",EpiSelectionAB],Table.loc["av",EpiSelectionAC],0)
        ac_alloy = Ternary(x,Table.loc["ac",EpiSelectionAB],Table.loc["ac",EpiSelectionAC],0)
        D0_alloy = Ternary(x,Table.loc["D0",EpiSelectionAB],Table.loc["D0",EpiSelectionAC],float(BowingParamD))
        b_alloy = Ternary(x,Table.loc["b",EpiSelectionAB],Table.loc["b",EpiSelectionAC],0)
        Cv_av = 3*(Table.loc["av",EpiSelectionAB]-Table.loc["av",EpiSelectionAC])*((Table.loc["Lattice Constant",EpiSelectionAB]-Table.loc["Lattice Constant",EpiSelectionAC])/Table.loc["Lattice Constant",SubSelection])
        Ev_av_alloy = Ternary(x,Table.loc["Ev_av",EpiSelectionAB],Table.loc["Ev_av",EpiSelectionAC],Cv_av)
        Eg_alloy = Ternary(x,Table.loc["Eg",EpiSelectionAB],Table.loc["Eg",EpiSelectionAC],float(BowingParamE))

                           
        TrnyTabledata = pd.Series([LatCon_alloy,C11_alloy,C12_alloy,av_alloy,ac_alloy,D0_alloy,b_alloy,Ev_av_alloy,Eg_alloy])
        TrnyTable=pd.DataFrame({str(EpiSelection):TrnyTabledata})
        TrnyTable.index =([   "Lattice Constant",  "C11",  "C12",   "av",   "ac",   "D0",   "b",   "Ev_av",   "Eg"])
        print(TrnyTable)
        


        
        Is_BG_wanted = input("Do you want a band gap diagram? (Y/N)")
        if Is_BG_wanted == "Y":

            LatCon_C1=(TrnyTable.loc["Lattice Constant",EpiSelection])
            LatCon_S1=(Table.loc["Lattice Constant",SubSelection])
            C11_C1=(TrnyTable.loc["C11",EpiSelection])
            C12_C1=(TrnyTable.loc["C12",EpiSelection])
            av_C1=(TrnyTable.loc["av",EpiSelection])
            ac_C1=(TrnyTable.loc["ac",EpiSelection])
            D0_C1=(TrnyTable.loc["D0",EpiSelection])
            b_C1=(TrnyTable.loc["b",EpiSelection])
            EV_av_C1=(TrnyTable.loc["Ev_av",EpiSelection])
            Eg_C1=(TrnyTable.loc["Eg",EpiSelection])


            EV_Epi,EC_Epi= bandgap(C11_C1,C12_C1,av_C1,ac_C1,D0_C1,b_C1,EV_av_C1,Eg_C1)
            print("Ev for "+ EpiSelection +":"+ str(EV_Epi))
            print("Ec for "+ EpiSelection +":"+ str(EC_Epi))

            C11_S1=(Table.loc["C11",SubSelection])
            C12_S1=(Table.loc["C12",SubSelection])
            av_S1=(Table.loc["av",SubSelection])
            ac_S1=(Table.loc["ac",SubSelection])
            D0_S1=(Table.loc["D0",SubSelection])
            b_S1=(Table.loc["b",SubSelection])
            EV_av_S1=(Table.loc["Ev_av",SubSelection])
            Eg_S1=(Table.loc["Eg",SubSelection])


            EV_Sub,EC_Sub= bandgap(C11_S1,C12_S1,av_S1,ac_S1,D0_S1,b_S1,EV_av_S1,Eg_S1)
            print("Ev for "+ SubSelection +":"+ str(EV_Sub))
            print("Ec for "+ SubSelection +":"+ str(EC_Sub))

            if EV_Sub<EV_Epi and EC_Sub>EC_Epi:

                Is_QW_wanted = input("Do you want a Quantum Well? (Y/N)")

            
                if Is_QW_wanted == "Y":
                    z=pd.Series(np.arange(-100,100,1))
                    EVplt_Sub = np.where((z>=-50)&(z<=0),EV_Sub,np.nan)
                    ECplt_Sub = np.where((z>=-50)&(z<=0),EC_Sub,np.nan)
                    plt.vlines(100,EV_Epi,EV_Sub,color="k")
                    plt.vlines(100,EC_Epi,EC_Sub,color="k")
                    plt.plot(EVplt_Sub,color="k")
                    plt.plot(ECplt_Sub,color="k")
                    EVplt_Epi = np.where((z<=50)&(z>=0),EV_Epi,np.nan)
                    ECplt_Epi = np.where((z<=50)&(z>=0),EC_Epi,np.nan)
                    plt.plot(EVplt_Epi,color="k")
                    plt.plot(ECplt_Epi,color="k")
                    EVplt_Sub = np.where(z>=50,EV_Sub,np.nan)
                    ECplt_Sub = np.where(z>=50,EC_Sub,np.nan)
                    plt.vlines(150,EV_Epi,EV_Sub,color="k")
                    plt.vlines(150,EC_Epi,EC_Sub,color="k")
                    plt.plot(EVplt_Sub,color="k")
                    plt.plot(ECplt_Sub,color="k")
                    EVplt_Sub = np.where((z>=-50)&(z<=0),EV_Sub,np.nan)
                    ECplt_Sub = np.where((z>=-50)&(z<=0),EC_Sub,np.nan)
                    plt.plot(EVplt_Sub,color="k")
                    plt.plot(ECplt_Sub,color="k")

                    
                    #plt.title()
                    #plt.ylim()
                    plt.ylabel("Energy")
                    plt.xlabel(SubSelection + "   -----   " +EpiSelection + "   -----   " + SubSelection)
                    #plt.legend()
                    plt.show()


                
            else:
                z=pd.Series(np.arange(0,100,1))
                EVplt_Epi = np.where(z<=50,EV_Epi,np.nan)
                ECplt_Epi = np.where(z<=50,EC_Epi,np.nan)
                plt.vlines(50,EV_Epi,EV_Sub,color="k")
                plt.vlines(50,EC_Epi,EC_Sub,color="k")
                plt.plot(EVplt_Epi,color="k")
                plt.plot(ECplt_Epi,color="k")
                EVplt_Sub = np.where(z>=50,EV_Sub,np.nan)
                ECplt_Sub = np.where(z>=50,EC_Sub,np.nan)
                plt.plot(EVplt_Sub,color="k")
                plt.plot(ECplt_Sub,color="k")
                
                #plt.title()
                #plt.ylim()
                plt.ylabel("Energy")
                plt.xlabel(EpiSelection+ "   -----   " + SubSelection)
                #plt.legend()
                plt.show()
        



        
            
if Type == "4":
    EpiSelection = input("What compound for the epilayer? ABCD")
    EpiSelectionAB = input("What compound for the epilayer AB, where the alloy is ABD? (e.g wanting wanting InPAs, AB is InP and AC is InAs")
    EpiSelectionAC = input("What compound for the epilayer AC?")
    EpiSelectionAD = input("What compound for the epilayer AD?")
    x = input("What is x as a decimal")
    y = input("What is y as a decimal")
    x=float(x)
    y=float(y)
    print("\n")
    SubSelection = input("What compound is substrate?")
    print("\n")
    print("\n")

    
    LatCon_alloy = Quarternary(x,y,(Table.loc["Lattice Constant",EpiSelectionAB]),(Table.loc["Lattice Constant",EpiSelectionAC]),(Table.loc["Lattice Constant",EpiSelectionAC]))
    C11_alloy = Quarternary(x,y,Table.loc["C11",EpiSelectionAB],Table.loc["C11",EpiSelectionAC],(Table.loc["C11",EpiSelectionAD]))
    C12_alloy = Quarternary(x,y,Table.loc["C12",EpiSelectionAB],Table.loc["C12",EpiSelectionAC],(Table.loc["C12",EpiSelectionAD]))

    QRTTabledata = pd.Series([LatCon_alloy,C11_alloy,C12_alloy])
    QRTTable=pd.DataFrame({str(EpiSelection):QRTTabledata})
    QRTTable.index =([   "Lattice Constant",  "C11",  "C12"])
    print(QRTTable)
    
    v=(QRTTable.loc["C12",EpiSelection])/((QRTTable.loc["C12",EpiSelection])+(QRTTable.loc["C11",EpiSelection]))
    print("Poission Ratio: " + str(v) + "\n")
    
    A=(QRTTable.loc["Lattice Constant",EpiSelection])
    a=A*(10**-10)
    print("Lattice Constant: "+ str(A) + "\n")
    
    f=abs((A-(Table.loc["Lattice Constant",SubSelection]))/(Table.loc["Lattice Constant",SubSelection]))
    print("Missfit Percentage: "+ str(f*100) + "% \n")


    CritThickness=2
    
    for i in range (99):
        if CritThickness > 0:
            CritThickness = ((1 - v)/(1 + v))*(1/(16 * pi* np.sqrt(2)))*((b**2)/a)*((1/((f)**2)))*(np.log(CritThickness/b))
            #CritThickness = abs(CritThickness)
        else:
            print("Critical thickness became negitive.")
            break
    CritThickness= CritThickness*10**9
    print(str(CritThickness) +"nm")

    





