# Also check SHOCKWAVE NOTES PDF File in the current directory
#FLOW TABLES ARE OBTAINED FROM URL :
#hTcp://www.aerodynamics4students.com/table4.php
#https://www.codeproject.com/Tips/864704/Python-Line-Intersection-for-Pygame
#https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
#http://www.aerodynamics4students.com/moc/moc_flow.php?problem=1&alt=0&gamma=1.4&Pco=500&tto=3150&mach=2.35&vplot=1&oPcion=Display
#https://tbc-python.fossee.in/convert-notebook/Modern_Compressible_Flow_with_historical_perspective_by_John_D_Anderson/4._Oblique_Shock_and_Expansion_Waves.ipynb

from mpl_toolkits import mplot3d
from scipy.integrate import odeint
import numpy as np
import math
from math import atan
from math import asin
from math import sqrt
from math import pow
import numpy as np
from pprint import pprint 
import matplotlib.pyplot as plt
import scipy.interpolate
import collections
import json
import csv

def clean_sin(sin_angle):
    return min(1,max(sin_angle,-1))
                     
def u_PM(M):
#M = Mach number of the gas
    M =clean_sin(1/M)
    return CDegrees(math.asin(M))
def CDegrees(Angle):
    return (180/math.pi)*Angle

def CRadians(Angle):
    return (math.pi/180)*Angle

def v_PM(Gama,M):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M= float(M)    
    v1 = math.sqrt((Gama+1.0)/(Gama-1.0))
    v2 = v1* CDegrees(math.atan(math.sqrt(abs(((Gama-1.000)/(Gama + 1.000))*((M ** 2)-1)))))
    v3 = v2-CDegrees(math.atan(math.sqrt(abs((M ** 2)-1.0))))
    return float(v3)

def M_No(Gama,M1,Beta,Theta):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M1 = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1= float(M1)
    Theta=float(Theta)
    Beta= float(Beta)
    if M1 and Gama and Theta and Beta and Beta-Theta:
        v1 = 1 + (((Gama-1)/2) * ( M1 ** 2 ) * ( math.sin(math.radians(Beta)) ** 2 ))
        v2 = (Gama * ( M1 ** 2 ) * ( math.sin(math.radians(Beta)) ** 2 ))-((Gama-1)/2)
        v3 = (1/math.sin(math.radians(Beta-Theta))) * math.sqrt(abs(v1/v2))
        return v3

def slope(p1, p2) :
   if (p2[0] - p1[0])  == 0:
        p2[0] = 1
   return (p2[1] - p1[1]) * 1. / (p2[0] - p1[0])
   
def y_intercept(slope, p1) :
   return p1[1] - 1. * slope * p1[0]

def line_intersect(m1, b1, m2, b2):
    if m1 == m2:
        print ("These lines are parallel!!!")
        return None
    # y = mx + b
    # Set both lines equal to find the intersection point in the x direction
    # m1 * x + b1 = m2 * x + b2
    # m1 * x - m2 * x = b2 - b1
    # x * (m1 - m2) = b2 - b1
    # x = (b2 - b1) / (m1 - m2)
    x = (b2 - b1) / (m1 - m2)
    # Now solve for y -- use either line, because they are equal here
    # y = mx + b
    y = m1 * x + b1
    return x,y

def intersect(line1, line2) :
    
   min_allowed = 1e-5   # guard against overflow
   big_value = 1e10     # use instead (if overflow would have occurred)
   m1 = slope(line1[0], line1[1])
   #print( 'm1: %d' % m1 )
   b1 = y_intercept(m1, line1[0])
   #print( 'b1: %d' % b1 )
   m2 = slope(line2[0], line2[1])
   #print( 'm2: %d' % m2 )
   b2 = y_intercept(m2, line2[0])
   #print( 'b2: %d' % b2 )
    
   if abs(m1 - m2) < min_allowed :
      x = big_value
   else :
      x = (b2 - b1) / (m1 - m2)
   y = m1 * x + b1
   y2 = m2 * x + b2
   #print( '(x,y,y2) = %d,%d,%d' % (x, y, y2))
   return float(x),float(y)

def segment_intersect(line1, line2) :
   intersection_pt = intersect(line1, line2)

   #print( line1[0][0], line1[1][0], line2[0][0], line2[1][0], intersection_pt[0] )
   #print( line1[0][1], line1[1][1], line2[0][1], line2[1][1], intersection_pt[1] )
   
   if (line1[0][0] < line1[1][0]) :
      if intersection_pt[0] < line1[0][0] or intersection_pt[0] > line1[1][0] :
         #print( 'exit 1' )
         return None
   else :
      if intersection_pt[0] > line1[0][0] or intersection_pt[0] < line1[1][0] :
         #print( 'exit 2' )
         return None
         
   if (line2[0][0] < line2[1][0]) :
      if intersection_pt[0] < line2[0][0] or intersection_pt[0] > line2[1][0] :
         #print( 'exit 3' )
         return None
   else :
      if intersection_pt[0] > line2[0][0] or intersection_pt[0] < line2[1][0] :
         #print( 'exit 4' )
         return None

   return intersection_pt

def combine_remove_and_sort(list1, list2):
    return sorted(list(set(list1+list2)))

def Beta_X(u,v):
    return (CDegrees(atan(u/v)))

def Mass_Rate(Gama,Pr_Throat,Temp_Throat,R_Gas,Area_Throat):
    Gama=float(Gama)
    Pr_Throat= float(Pr_Throat)
    Temp_Throat=float(Temp_Throat)
    R_Gas= float(R_Gas)
    Area_Throat= float(Area_Throat)
    v1 = ((Area_Throat * Pr_Throat)/math.sqrt(abs(Temp_Throat)))
    v2 = math.sqrt(Gama/R_Gas)
    v3 = ( (Gama + 1 ) / 2 ) ** (-1 * ((Gama + 1 )/ (2 * ( Gama - 1 ))))
    return v1 * v2 * v3

def Area_Ratio(Gama , Mach_Exit):
    #returns Area Ratio Ae / At
    Gama= float(Gama)
    Mach_Exit=float(Mach_Exit)
    v1 = ((Gama + 1 ) / 2) ** (-1 * ((Gama + 1 )/ (2 * ( Gama - 1 ))))
    v2 = ((1 + ((( Gama - 1 )/2) * (Mach_Exit ** 2) ))** ((Gama + 1 )/ (2 * ( Gama - 1 ))))/Mach_Exit
    return v1 * v2 
    
def Temp_Ratio(Gama,Mach_Exit):
    #returns Area Ratio Te / Tc
    Gama= float(Gama)
    Mach_Exit=float(Mach_Exit)
    v1 = ( 1 + ((Gama - 1) / 2) * (Mach_Exit ** 2 ) ) ** (-1)
    return v1

def Pr_Ratio(Gama,Mach_Exit):
    #returns Area Ratio Pe / Pt
    Gama= float(Gama)
    Mach_Exit=float(Mach_Exit)
    v1 = ( 1 + ((Gama - 1) / 2) * (Mach_Exit ** 2 ) ) ** (-Gama/(Gama-1))
    return v1

def PoP_Ratio(Gama,M1):
# Po is the stagnation pressure of the stream
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1=float(M1)
    v1 =  (1+ (((Gama-1)/2)*(M1**2)))
    v2 =  math.pow(v1,(Gama/(Gama-1)))
    return v2

def Exit_Velocity(Gama,Mach_Exit,Temp_Exit,R_Gas):
    #returns Area Ratio Pe / Pt
    Gama= float(Gama)
    R_Gas= float(R_Gas)
    Mach_Exit=float(Mach_Exit)
    Temp_Exit=float(Temp_Exit)
    v1 = Mach_Exit * math.sqrt(abs(Gama * R_Gas * Temp_Exit))
    return v1

def Exit_Thrust(Mass_Rate,Exit_Velocity,Exit_Pr,P_Atm,Exit_Area):
    Mass_Rate= float(Mass_Rate)
    Exit_Velocity= float(Exit_Velocity)
    Exit_Pr=float(Exit_Pr)
    P_Atm=float(P_Atm)
    Exit_Area=float(Exit_Area)

    v1 = (Mass_Rate * Exit_Velocity) + ((Exit_Pr - P_Atm) * Exit_Area)
    return v1

def Pr_Ratio(Gama,M1,Beta):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1=float(M1)
    Beta= float(Beta)
    v1 =  1 + ((( 2 * Gama)/(Gama + 1 )) * ((( M1 ** 2 ) * ( math.sin(math.radians(Beta)) ** 2 )) -1 ))
    #P2/P1 = v1
    return v1


def D_Ratio(Gama,M1,Beta):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1=float(M1)
    Beta= float(Beta)
    v1 = (Gama + 1)*(M1 ** 2 ) * ( math.sin(math.radians(Beta)) ** 2 )
    v2 =  ((Gama - 1)*(M1 ** 2 ) * ( math.sin(math.radians(Beta)) ** 2 )) + 2
    # D1/D1 = v1/v2
    return v1/v2

def Do_Ratio(Gama,M1):
# Do = Stagnation Density of stream
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1=float(M1)
    v1 =  (1+ (((Gama-1)/2)*(M1**2)))
    v2 =  math.pow(v1,(1/(Gama-1)))
    return v2

def T_Ratio(Gama,M1,Beta):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1=float(M1)
    Beta= float(Beta)
    v1 =  Pr_Ratio(Gama,M1,Beta)
    v2 =  D_Ratio(Gama,M1,Beta)
    v3 = v1/v2
    # T2/T1 = P1/P1 x D1/D2 = v3
    return v3

def To_Ratio(Gama,M1):
# To = Stagnation Temperature of stream
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
    Gama=float(Gama)
    M1=float(M1)
    v1 =  (1+ (((Gama-1)/2)*(M1**2)))
    return v1

def Cf(Gama,Pc,Pe,Patm,Rt,Re):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# Pc = Combustion chamber Pressue
# Pe = Exit Pressure of expander nozzle 
# Patm = Atmospheric Pressure at that altitude
# Rt = Throat Area
# Re = Exit Area
    Gama=float(Gama)
    Pc=float(Pc)
    Pe=float(Pe)
    Rt=float(Rt)
    Re=float(Re)

    exp1 =  abs(( 2 * (Gama ** 2))/(Gama - 1))
    exp2 = abs( 2 /(Gama + 1)) ** (( Gama + 1 ) / (Gama - 1 ))
    exp3 =abs (1 - ( ( Pe / Pc ) ** (( Gama - 1) / Gama )))


    Area_Throat = (math.pi * ( Rt ** 2 ))
    Area_Exit = (math.pi * ( Re ** 2 ))

    exp4 = ( (Pe - Patm)* Area_Exit )
    exp5 = ( Pc * Area_Throat )
    
    value = math.sqrt( exp1 * exp2 * exp3 ) + ( exp4 / exp5 )

    return value

def Co_Ratio(Gama,M1):
# Co / C = To / T = Co is Stagnation Speed Of Sound
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required

    Gama=float(Gama)
    M1=float(M1)
    v1 =  (1+ (((Gama-1)/2)*(M1**2)))
    v2 =  math.sqrt(v1)
    return v2

def P01P02_Ratio(Gama,M1):
# G = Gama Adiabatic work doing index of the gas = cp/cv
# M = Mach Number of the Gas for which prandtl meyer expansion function v1 is required
   if M1 and Gama:
        Gama=float(Gama)
        M1=float(M1)
        v1 =  (((Gama+1)/2)*(M1 ** 2)/(1 + ((Gama-1)/2)*(M1 ** 2)))
        v1X=math.pow(v1,(Gama/(Gama-1)))
        v2 =  (((2*Gama)/(Gama+1))*(M1 ** 2))-((Gama-1)/(Gama+1))
        v2X = math.pow(1/v2,1/(Gama-1))
        return v1X*v2X
    
def Calc_Ideal_MachNo(Gama,M1,pt):
    p = P01P02_Ratio(Gama,M1)* pt
    q = (Gama/2) * p * (M1**2)
    s = (2 / (Gama -1 ))
    pwr=(Gama- 1)/Gama
    t =( (( q / p) + 1 )**pwr) - 1
    Machno= math.sqrt(s * t)
    print(Machno)
    
def Calc_Throat_Sound_VelocityRatio(Gama):
    # C_Throat / Co
    return (2/(Gama+1))

def Calc_Throat_Pr_Ratio(Gama):
    # P_Throat / Po
    Power = (Gama/(Gama-1))
    EXR = ((2/(Gama+1)) ** Power)
    return EXR

def Calc_Throat_T_Ratio(Gama):
    # T_Throat / To
    EXR = (2/(Gama+1))
    return EXR

def Calc_Throat_D_Ratio(Gama):
    # D_Throat / Do
    Power = (1/(Gama-1))
    EXR = ((2/(Gama+1)) ** Power)
    return EXR

    
    
    
#================================================================================
#   NOZZLE DESIGN DATA
#================================================================================
#COMBUSTION PROPERTIES
Pc = 500 * (101325/14.7)    #Bar 500/10
Tc = 2000                 # K
R_Gas = 287                #J / Kgm - K
Dc = 1.225                  #Kgm/m3
P_Atm = 101325              # Pascal N /m2 
#scale factor for the nozzle calculations
Scale= 1.00
#Gama , Adiabatic bulk mod of elasticity
G = 1.4


#Number of Points to Draw Contour on is given be N_X
N_X= int(22/ Scale) #33 140

#THE INTERVAL FOR THRUST EVALUATION FOR NOZZLE [ METHOD OF CHARCTERISTICS ]
N_K = N_X - (int(N_X/2)- 2)

#instrument movement inside nozzle to find shock wave conditions
KX = 0 # Any Number between 1 to 100[or higher]

#Mach Number M_Sonic of Nozzle [ Enterance ]
#ASSUMED AVERAGE THROAT MACH NUMBER , M = 1.1 [ CHOCKED NOZZLE ]
# P = G . h => h = P / G
#Density Air = P = D . R . T
D_Air = ( Pc / (R_Gas * Tc))
#Pressure head in terms of manometer height , H_max
H_max= Pc/ ( G * D_Air )
#Free Stream Velocity before nozzle Gets Choked
# v2  2 . g = h => V_Max = sqrt( 2 . g . H_Max )
V_Max = sqrt(2 * 9.81 * H_max)

#Sound Velocity , C 
# C = Sqrt ( G .  R . T )
C_Max = sqrt(G * R_Gas * Tc)

# Mach Number  M_Sonic = V_Max / C
M_Sonic = float(V_Max / C_Max )
if ( M_Sonic > 1 ):
    M_Sonic = float(1.0001)

print(M_Sonic)
    
M2_No = 0
Calc_Ideal_MachNo(G,0.8,Pc)
#Exit Mach Number M Of Nozzle 
M_Exit = float( 2.0 / Scale) #2.355

# Deflection Angle of the given contour
Theta_Max = 0  # Degrees

#μ=Mach angle, (deg), (sin-1(1/M1)) of minimal oblique wave for M1
ux = 0

#Prandtl Meyer Expansion Angle After Shock Wave or Expansion Wave , v1
#Prandtl Meyer Angle at Sonic Line or Enterance = CONSTANT 
v_M=0         #Degrees
v_PMA  = float(v_PM(G,M_Sonic))
Beta =v_PMA

#Prandtl Meyer Expansion Angle After Shock Wave or Expansion Wave , v2
# theta = v_M2 - v_M1
v_M2=0         #Degrees

v_M = v_PM(G,M_Exit)
ux = u_PM(M_Exit)



Theta_Max =  float( v_M/2 )

#===================================================================================
#NOZZLE GRID DEVELOPMENT
#L/R = Length of expander from throat to Throat Radius ratio : 308 /23.8 = 12.94  
#=================================================================================
L_Rt_Ratio = 6.07
D_Throat = 2 #foot
R_Throat = float(((D_Throat/2) * 12) * 25.4)/Scale # foot to mm
A_Throat = (math.pi * (R_Throat/1000) * (R_Throat/1000))
#ideal nozzle radius , Pc = Ft / At , At = Ft / Pc  = 

X_Max = L_Rt_Ratio
X_Max = ((X_Max *12)*25.4)/Scale #foot to mm
Exit_Angle= 0
RcRt= 2
R_C = RcRt * R_Throat 
X_MX = [-1,
0,
0.01,
0.02,
0.03,
0.04,
0.05,
0.06,
0.07,
0.08,
0.09,
0.1,
0.11,
0.12,
0.13,
0.14,
1.3353,
1.7003,
1.9942,
2.2583,
2.5072,
2.7502,
2.9913,
3.2323,
3.4754,
3.722,
3.9732,
4.2294,
4.4902,
4.7424,
5.074,
]
Y_MX1 = [1,
1,
1.0001,
1.0003,
1.0006,
1.0011,
1.0016,
1.0022,
1.003,
1.0039,
1.0048,
1.0059,
1.0071,
1.0084,
1.0098,
1.0113,
1.2047,
1.2597,
1.3008,
1.3349,
1.3643,
1.3904,
1.4137,
1.4343,
1.4526,
1.4684,
1.4818,
1.4928,
1.5012,
1.5066,
1.5101,
]

Xa_M = []
Ya_M = []
X_M  = []
Y_M1 = []
Y_M2 = []

i = 0 

#NOTE :
# if we use 2 x N_X THEN THE FLOW PATcEREN CHANGES THROUGH OUT THE NOZZLE
# if we use N_X then the flow patteren is different Through the nozzle
# N_X is the minimum value for the interval of x given below

#Theta_Max = abs( math.degrees(math.atan(Y_MX1[0]/X_MX[0])))+0.5
Angle = float(Theta_Max )
x = np.linspace(0,X_Max,int(N_X + KX *N_X))
b=math.tan(math.radians(Angle))-math.tan(math.radians(Exit_Angle)) /2 
a= -(b / (2 * X_Max))
c= -math.pow(b,2)/4*a
pprint("Theta_Max : " + str(Angle))
pprint("a : " + str(a))
pprint("b : " + str(b))
pprint("c : " + str(c))

Angle_Circle= np.linspace(0.0,Angle,int(N_X/2 + KX * N_X/2))
for ar in Angle_Circle:
    xm= (R_C * math.sin(math.radians(ar)))
    y= (R_C + (R_C - R_C * math.cos(math.radians(ar))))/2.07
    print("x(mm) : " + str(xm) + " , " + "y(mm) : " + str(y))
    X_M.append(xm)
    Y_M1.append(y)
    Y_M2.append(-y)
    Xa_M.append(xm)
    Ya_M.append(y)


#Xa = -1*25.4*12
#Ya = 1*25.4*12
Xa = xm
Ya = y
i = 0
for xi in x:
    xm1= Xa + xi
    ym1= Ya + (((a* xi**2) + (b * xi) + c))/2.07
    if(xm1 <= X_Max +1  ):
        X_M.append(xm1)
        Y_M1.append(ym1)
        Y_M2.append(-ym1)
        Xa_M.append(xm)
        Ya_M.append(y)
        print("x(mm) : " + str(xm1) + " , " + "y(mm) : " + str(ym1))
            

    

#fig = plt.figure()
#ax = plt.axes(projection='3d')


#ax.plot3D(y,Uxy,Uxy)
plt.plot(X_M,Y_M1, 'ro')
plt.plot(X_M,Y_M2, 'ro')

i = 0
#Points for Analysis
Y_r1 = []
Y_r2 = {}
Y_r3 = {}
X_r1 = []

#two lines have points
AX1 = []
BY1 = []

AX2 = []
BY2 = []
C = []

#Angle of Grid 
Grid_Angles = []
ThetaX = 0
Theta0 = 0
i = 0
j = 0
m = 0
l = 0
xx1,xx2 =0,0
yy1,yy2 =0,0

x1,x2 =0,0
y1,y2 =0,0

XXT1 =  0
YYT1 =  0
XXP1,YYP1=0,0
XTR  = []
YTR  = []
C    = []

xxxx1 = []
yyyy1 = []

xxxx2 = []
yyyy2 = []

xxxx1X = []
yyyy1X = []

xxxx2X = []
yyyy2X = []

def connectpoints(x,y,p1,p2):
    
    k = N_X      #Points  
    x1, x2 = x[p1], x[p2]
    y1, y2 = y[p1],y[p2]
    #plt.plot([x1,x2],[y1,y2],'k-')

    #JOINING POINTS OF EXPANSION FAN TO SONIC LINE
    #plt.plot([0,x2],[y2,0],'-b')
    #plt.plot([0,x2],[-y2,0],'-b')
    
    #DRAWING SONIC LINE OF THE NOZZLE 
    plt.plot([0,0],[0,y2],'--g')
    plt.plot([0,0],[0,-y2],'--g')
    plt.plot([0,x[-1]],[0,0],'-b')

    #change value of j and s for finding intersection points down stream of nozzle
    #and 0 for throat region calculations 
    j  = int(KX/2)
    i  = 0
    r = 1
    jx = 0
    for j in range(0,len(X_M),KX+1):
       if ( k <= 0 ):
            break
       if (X_M[j] <= Xa ):
              
       #LEFT RUNNING CHARACTERISTICS C- : K _
            xx1= X_M[j]
            yy1 = Y_M1[j]
            xx2=   X_M[j + abs(int(len(X_M)/5+ KX)) ]
            yy2=  Y_M2[j + abs(int(len(X_M)/5+ KX)) ]
            #print(j ,",",j +  int(len(X_M)/2) )
            AX1=[xx1 , xx2 ]
            BY1=[yy1 , yy2 ]
            plt.plot([xx1,xx2],[yy1,yy2])
              
    #RIGHT RUNNING CHARACTERISTICS C+ : K +    
            x1=X_M[j]
            y1=Y_M2[j]
            x2= X_M[j + abs(int(len(X_M)/5+ KX)) ]
            y2= Y_M1[j + abs(int(len(X_M)/5+ KX)) ]
            AX2 = [x1 , x2]
            BY2 = [y1 , y2]
            plt.plot([x1,x2],[y1,y2])

            
                    
                  
    #POINT OF INTERSECTION OF TWO CHARACTERISTICS
            x_interp = scipy.interpolate.interp1d(BY1, AX1)
            # X intercePc where Line crosses the axis , Y = 0
            XXP1 = x_interp(0) 
            YYP1=0

    #change value of j and s for finding intersection points down stream of nozzle
    #and 0 for throat region calculations 
            s = int(KX/2)
             
    #LEFT RUNNING CHARACTERISTICS C- : K _
            for s in range(0,int(len(X_M)/2),KX+1 ):
                if (X_M[s] <= Xa ):
                    xxx1= X_M[s]
                    yyy1 = Y_M1[s]
                    xxx2=   X_M[s + abs(int(len(X_M)/5+ KX )) ]
                    yyy2=  Y_M2[s + abs(int(len(X_M)/5+ KX )) ]

                    xxxx1.append(xxx1)
                    yyyy1.append(yyy1)
                    
                    xxxx2.append(xxx2)
                    yyyy2.append(yyy2)
                    #plt.plot([xxx1,xxx2],[yyy1,yyy2],'-g')
                    
                    xR1=X_M[j]
                    yR1=Y_M2[j]
                    xR2= X_M[j + abs(int(len(X_M)/5+ KX)) ]
                    yR2= Y_M1[j + abs(int(len(X_M)/5+ KX)) ]
            

                    xxxx1X.append(xR1)
                    yyyy1X.append(yR1)
                        
                    xxxx2X.append(xR2)
                    yyyy2X.append(yR2)
                           
                    A1 = [xxx1,yyy1]
                    A2 = [xxx2,yyy2]
           
                    B1 = [x1,y1]
                    B2 = [x2,y2]
                  

                    #print(A1,B1,A2,B2)
                    
                    #slope_A = slope(A1, A2)
                    #slope_B = slope(B1, B2)
                    
                    #y_int_A = y_intercePc(A1, slope_A)
                    #y_int_B = y_intercePc(B1, slope_B)
                    if segment_intersect((A1,A2),(B1,B2)) is not None:
                        intersect_point=segment_intersect((A1,A2),(B1,B2) )                   
                        C.append(intersect_point)
                        #plt.plot(intersect_point[0],intersect_point[1],'bo')

                    C1 = [xx1,yy1]
                    C2 = [xx2,yy2]
           
                    D1 = [x1,y1]
                    D2 = [x2,y2]
                  
                    if segment_intersect((C1,C2),(D1,D2)) is not None:
                        intersect_point=segment_intersect((C1,C2),(D1,D2) )
                        C.append(intersect_point)
                        #plt.plot(intersect_point[0],intersect_point[1],'ro')
                    
                    dY = np.linspace(0,R_Throat,len(X_M)*10)
                    E1 = [X_M[s],Y_M1[s]]
                    E2 = [X_M[s+1],dY[s]]
               
                    F1 = [xxxx1[j],yyyy1[j]]
                    F2 = [xxxx2[j],yyyy2[j]]
                    #plt.plot([x1,xx1],[y1,yy1],'r')
                    #plt.plot([x1,x2],[y1,y2],'r')
                   
                    if segment_intersect((E1,E2),(F1,F2)) is not None:
                        intersect_point=segment_intersect((E1,E2),(F1,F2) )
                        C.append(intersect_point)
                        #plt.plot(intersect_point[0],intersect_point[1],'ro')

                    dY = np.linspace(0,R_Throat,len(X_M)*10)
                    G1 = [X_M[s],Y_M2[s]]
                    G2 = [X_M[s+1],dY[s]]
               
                    H1 = [xxxx1X[j],yyyy1X[j]]
                    H2 = [xxxx2X[j],yyyy2X[j]]
                                         
                        
                    if segment_intersect((G1,G2),(H1,H2)) is not None:
                        intersect_point=segment_intersect((G1,G2),(H1,H2) )
                        C.append(intersect_point)
                        plt.plot(intersect_point[0],intersect_point[1],'ro')
                                    
                                          
                    #pprint(line_intersect(slope_A, y_int_A, slope_B, y_int_B))
                         
                    #s = s + 1
                      
            #jx = jx + 1             
            AX=np.linspace(0,x2,k)
            AY=np.linspace(0,y2,k)

            #plt.plot(XXP1,YYP1,'bo')
            
            
            Theta0=Angle - math.floor(Angle)
            ThetaM= math.ceil(math.floor(Angle) / (len(X_M)-KX))
            #ThetaX=math.floor(Angle)/k
            
            
            i = 0
            r = 1              
            for axt in AX:
                if( Theta0 + (j * ThetaM) > Angle):
                    break
                elif ( i == 0 ):
                    #Y_r1.append([ chr(65+j) + "-" + str(i), axt, AY[i],0 ])
                    pass
                elif ( i == 1 and j == 0):
                    Theta0=Angle - math.floor(Angle)
                    Y_r1.append([chr(65+j) + "-" + str(i),axt,AY[i],Theta0])
                    X_r1=[axt,AY[i]]
                elif ( i == 1 and j > 0):
                    Theta0=Angle - math.floor(Angle)
                    ThetaX=np.linspace(0,math.floor(Angle),k)
                    X_r1=[axt,AY[i]]
                    Y_r1.append([chr(65+j) + "-" + str(i),axt,AY[i],Theta0 + (j * ThetaM)])
                    
                
                    
                else:
                    #Y_r1.append([chr(65+j) + "-" + str(i),axt,AY[i],0])
                    pass
                
                i = i + 1
            
             
            #k = k - 1
            #j = j  + 1
           
def connectpointsX():
    k = N_X      #Points  
    j  = 0
    i  = 0
    s = 0 

     #LEFT RUNNING CHARACTERISTICS C- : K _
    xx1= X_M[j]
    yy1 = Y_M1[j]
    xx2=   X_M[j +  N_K ]
    yy2=  Y_M2[j +  N_K ]
    AX1=[xx1 , xx2 ]
    BY1=[yy1 , yy2 ]
    plt.plot([xx1,xx2],[yy1,yy2])
              
    #RIGHT RUNNING CHARACTERISTICS C+ : K +    
    x1=X_M[j]
    y1=Y_M2[j]
    x2= X_M[j + N_K]
    y2= Y_M1[j +  N_K]
    AX2 = [x1 , x2]
    BY2 = [y1 , y2]
    plt.plot([x1,x2],[y1,y2])
            
    #POINT OF INTERSECTION OF TWO CHARACTERISTICS
    x_interp = scipy.interpolate.interp1d(BY1, AX1)
    # X intercePc where Line crosses the axis , Y = 0
    XXP1 = x_interp(0) 
    YYP1=0

    dx1=  XXP1 
    dy1 = (yy1-y1)
    dz1 = math.sqrt((dx1 ** 2) + ( dy1 ** 2 ))
    PcX = np.linspace(0,dx1,N_X)
    PTY = np.linspace(0,dy1,N_X)
    PTZ = np.linspace(0,dz1,N_X)
    
    pprint(dy1)
    #plt.plot([xx1,XXP1],[yy1,YYP1],'--b')
    #plt.plot([xx2,XXP1],[yy2,YYP1],'--b')
    #plt.plot(XXP1,YYP1,'ro')

   
print("System is Collecting Intersection Points .... plz wait ")            
connectpoints(X_M,Y_M1,0,1)
print("Intersection Points Collected .... Finding Thermodynamic Properties .... plz wait ")
#pprint(Y_r1)
ch  =  0
v   =  0
div =  []
for i in range(0,len(C),1):
    if(C[i][1] >= 0 ):
        if ( C[i][1] == 0 ):
            ch = ch + 1
            div.append(v)
            v = 0
        else:
            div.append(1)
            
    
    v = v + 1

v      = 1
p      = 0 
ch     = 0 
ThetaM = 0
count  = 0

Mach_Avg = 0
M_Rate_Avg = 0
V_Exit_Avg = 0
P_Exit_Avg = 0
T_Exit_Avg = 0
ThrustX_Avg = 0

Mach_No_Array = 0
M_Rate_Array = 0
V_Exit_Array = 0
P_Exit_Array = 0
T_Exit_Array = 0
ThrustX_Array = 0
Correction_Factor_Avg = 0

Q = 0
R=  0


# YOU CAN USE REQUIRED MACH NUMBER [DESIGN MACH NO ] FOR MAX EXIT PRESSURE OF NOZZLE
Pr_X=abs(Pr_Ratio(G,M_Exit,Beta)* Pc)
R_Exit = Y_M1[ len(X_M)-1 ]
#MAX AREA RATIO FOR THESE THEMODYNAMIC PROPERTIES CAN BE CALCULATED USING EPSILON = A_Exit /A_throat
#area ratio formula is given above in the function A_Ratio
# A_Ratio = A_Exit / A_Throat
#A_Throat = A_Exit / A_Ratio
#PI . R ^ 2 = A_Exit / A_Ratio
#R_Throat = sqrt ( A_Exit / A_Ratio x PI )
#Area Ratio [ MAX ] epsilon = A_Exit/A_throat
A_Ratio = abs(Area_Ratio(G , M_Exit))
A_Ratio_MAX = A_Ratio
A_Exit = math.pi * (( R_Exit/1000)**2 )
R_t_MAX = math.sqrt(A_Exit/(A_Ratio * math.pi))

#MINIMUM LENGTH OF THE NOZZLE CAN BE FOUND BELOW :
L_Min = X_M[ len(X_M)-1 ]

#You can select a Mach Number and Area Correction from 1 to 2.5
#THRUST CORRECTION FACTOR Cf can be calculated as below
#Correction_Factor = Cf(G,Pt,Pr_X , P_Atm, R_Throat/1000 , R_Exit/1000 )
#if Correction_Factor >= 2.5 :
#    Correction_Factor = float(2.5)

#--------------------------------------------------------------------------------------
#
#
#
# THROAT PROPERTIES CALCULATION
#
#
#
#
#-------------------------------------------------------------------------------
C_Throat = Calc_Throat_Sound_VelocityRatio(G)*C_Max
Pr_Throat = Calc_Throat_Pr_Ratio(G)*Pc
T_Throat = Calc_Throat_T_Ratio(G) * Tc
D_Throat = Calc_Throat_D_Ratio(G) * Dc

V_Throat = Exit_Velocity(G,M_Sonic,T_Throat,R_Gas)
M_Sonic = V_Throat/C_Throat
A_Ratio_Throat = Area_Ratio(G , M_Sonic)
A_Throat_X = A_Ratio_Throat * A_Throat

M_Rate = Mass_Rate(G,Pr_Throat,T_Throat,R_Gas,A_Throat_X)
Correction_Factor = Cf(G,Pc,Pr_Throat , P_Atm, R_Throat/1000 , R_Exit/1000 )
Thrust_Throat=Exit_Thrust(M_Rate,V_Throat,Pr_Throat,P_Atm,A_Throat_X)

print("Velocity Sound :" , C_Throat)
print("Pr Throat :" ,Pr_Throat)
print("Temperature Throat :" ,T_Throat)
print("Density Throat :" ,D_Throat)

print("Velocity Throat :" ,V_Throat)
print("Mach No Throat :" ,M_Sonic)
print("A Ratio Throat :" ,A_Ratio_Throat)
print("Area Throat :" ,A_Throat_X)

print("Mass Rate :" ,M_Rate)
print("Correction Factor Cf :" ,Correction_Factor)
print("Thrust at Throat  :" ,Thrust_Throat)

Correction_Factor = float(Correction_Factor)
Correction_FactorX2 = float(Correction_Factor)
print("\n\n")
#--------------------------------------------------------------------------------------
#
#
#
# THROAT STAGNATION PROPERTIES CALCULATION
#
#
#
#
#-------------------------------------------------------------------------------
C_stag_Throat = Co_Ratio(M_Sonic,G)*C_Max
V_stag_Throat = Exit_Velocity(G,M_Sonic,Tc,R_Gas)
M_stag_Throat = V_stag_Throat/C_stag_Throat

T_stag_Throat = To_Ratio(M_stag_Throat,G) * Tc
A_Ratio_Throat = Area_Ratio(G , M_stag_Throat)
A_Throat_X = A_Ratio_Throat * A_Throat

Pr_stag_Throat = (1/PoP_Ratio(M_stag_Throat,G))*Pc

D_stag_Throat = Do_Ratio(M_stag_Throat,G) * Dc


M_stag_Rate = Mass_Rate(G,Pr_stag_Throat,T_stag_Throat,R_Gas,A_Throat_X)
Correction_stag_Factor = Cf(G,Pc,Pr_stag_Throat , P_Atm, R_Throat/1000 , R_Exit/1000 )
Thrust_stag_Throat=Exit_Thrust(M_stag_Rate,V_stag_Throat,Pr_stag_Throat,P_Atm,A_Throat_X)

print("stag Velocity Sound :" , C_stag_Throat)
print("stag Pr Throat :" ,Pr_stag_Throat)
print("stag Temperature Throat :" ,T_stag_Throat)
print("stag Density Throat :" ,D_stag_Throat)

print("stag Velocity Throat :" ,V_stag_Throat)
print("stag Mach No Throat :" ,M_stag_Throat)
print("stag A Ratio Throat :" ,A_Ratio_Throat)
print("stag Area Throat :" ,A_Throat_X)

print("stag Mass Rate :" ,M_stag_Rate)
print("stag Correction Factor Cf :" ,Correction_Factor)
print("stag Thrust at Throat  :" ,Thrust_stag_Throat)

C_Throat = C_stag_Throat
Pr_Throat = Pr_stag_Throat
T_Throat = T_stag_Throat
D_Throat = D_stag_Throat

V_Throat = V_stag_Throat
M_Sonic = V_Throat/C_Throat
A_Ratio_Throat = Area_Ratio(G , M_Sonic)
A_Throat_X = A_Ratio_Throat * A_Throat

M_Rate = M_stag_Rate
Thrust_Throat=Thrust_stag_Throat

#REMOVING DUPLICATE INTESECTION POINTS FROM LIST 'C'
C = list(dict.fromkeys(C))


for i in range(0,len(C),KX+1):
    
    #if(C[i][1] >= 0 ):
       
        if ( i == 0):
            
            
            if ( v_PMA <= 0 or Beta <= 0  ):
                v_PMA = Beta
                Beta = v_PMA
            if M2_No is None :    
                M2_No = float(M_Sonic)
                
            v_PMA = v_PM(G,M2_No)    
            Theta0=Angle - math.floor(Angle)
            Q = v_PMA + Theta0
            R= v_PMA - Theta0
            Theta =( Q - R ) / 2
            Beta = ( Q + R ) / 2
            if( Theta <= 0 ):
                Theta = Theta0

            M2_No = M_Sonic    
            M2_No = abs(M_No(G,M2_No,Beta,Theta))
            if(M2_No >=  M_Exit):
                M2_No =   2 * M_Exit
                
            
            v_PMA = v_PM(G,M2_No)
            Pr_X=Pr_Ratio(G,M2_No,Beta) * Pr_Throat
            D_X = D_Ratio(G,M2_No,Beta) * D_Throat
            Temp_X=T_Ratio(G,M2_No,Beta) * T_Throat
            V_Exit = Exit_Velocity(G,M2_No,Temp_X,R_Gas)
            A_Ratio = Area_Ratio(G , M2_No)
            if(A_Ratio > A_Ratio_MAX):
                A_Ratio = A_Ratio_MAX
                         
            # Parametric area for flow through conic section of shock wave ,
            # Area = 2 x Pi x r x dr
            #2 and Pi are constants of integration so integration of r.dr
            # gives , Flow Area = 2 . Pi . [r^2 /2 ] = Pi . r ^ 2
            
            A_Exit = A_Ratio * A_Throat
            
            M_Rate = Mass_Rate(G,Pr_X,Temp_X,R_Gas,A_Throat)
            Correction_Factor = Cf(G,Pc,Pr_X , P_Atm, R_Throat/1000 , R_Exit/1000 )
            ThrustX=(Exit_Thrust(M_Rate,V_Exit,Pr_X,P_Atm,A_Exit))

            Mach_Avg=(Mach_Avg + M2_No)
            V_Exit_Avg = V_Exit_Avg + V_Exit   
            M_Rate_Avg = M_Rate_Avg + M_Rate
            P_Exit_Avg = P_Exit_Avg + (Pr_X)
            T_Exit_Avg = T_Exit_Avg + (Temp_X)
            ThrustX_Avg = ThrustX_Avg + ThrustX
            Correction_Factor_Avg = Correction_Factor_Avg + Correction_Factor
            
            Y_r2[p] = { chr(65+ch)+ "-"+ str(count): { "x" : C[i][0],
                                   "y" : C[i][1],
                                   "Theta": Theta0,
                                                                                           
                                   }
                     }
            Y_r3[p] = { chr(65+ch)+ "-"+ str(count): { "x" : abs(C[i][0]),
                                   "y" : abs(C[i][1]),
                                   "Q" : Q,
                                   "R" : R,
                                   "Theta": ( Q - R )/2,
                                   "Beta": ( Q + R )/2,
                                   "v_PM": v_PMA,
                                   "Mach":M2_No,
                                   "Pr_Exit ": Pr_X,
                                   "Temp Exit":Temp_X,
                                   "V_Exit ": V_Exit ,
                                   "Epsilon=A_e/A_t" : A_Ratio,
                                   "Exit_Area":A_Exit,
                                   "M_Rate ": M_Rate,
                                   "Thrust " :ThrustX,
                                   }
                     }
            
            #pprint(str(C[i][0]) + "," + str(C[i][1]))
            plt.plot(C[i][0],C[i][1],'yo')
            p = p + 1
            #count = count + 1 

        
            
        elif ( i > 0 and C[i][1] == 0 ):
            ch = ch + 1
            count = 0
            if ( v_PMA <= 0 or Beta <= 0  ):
                v_PMA = Beta
                Beta = v_PMA
            if M2_No is None :    
                M2_No = float(M_Sonic)
            
            v_PMA = v_PM(G,M2_No)    
            Theta0=Angle - math.floor(Angle)
            Q = v_PMA + Theta0
            R= v_PMA - Theta0
            Theta =( Q - R ) / 2
            Beta = ( Q + R ) / 2
            if( Theta <= 0 ):
                Theta = Theta0
            M2_No = M_Sonic
            M2_No = abs(M_No(G,M2_No,Beta,Theta))
            if(M2_No >=  M_Exit):
                M2_No =   2 * M_Exit
                
                
            v_PMA = v_PM(G,M2_No)
            Pr_X=Pr_Ratio(G,M2_No,Beta) * Pr_Throat
            D_X = D_Ratio(G,M2_No,Beta) * D_Throat
            Temp_X=T_Ratio(G,M2_No,Beta) * T_Throat
            V_Exit = Exit_Velocity(G,M2_No,Temp_X,R_Gas)
            A_Ratio = Area_Ratio(G , M2_No)
            if(A_Ratio > A_Ratio_MAX):
                A_Ratio = A_Ratio_MAX
                         
            # Parametric area for flow through conic section of shock wave ,
            # Area = 2 x Pi x r x dr
            #2 and Pi are constants of integration so integration of r.dr
            # gives , Flow Area = 2 . Pi . [r^2 /2 ] = Pi . r ^ 2
            
            A_Exit = A_Ratio * A_Throat
            M_Rate = Mass_Rate(G,Pr_X,Temp_X,R_Gas,A_Throat)
            Correction_Factor = Cf(G,Pc,Pr_X , P_Atm, R_Throat/1000 , R_Exit/1000 )
            ThrustX=(Exit_Thrust(M_Rate,V_Exit,Pr_X,P_Atm,A_Exit))
                               
            Mach_Avg=(Mach_Avg + M2_No)
            M_Rate_Avg = M_Rate_Avg + M_Rate
            V_Exit_Avg = V_Exit_Avg + V_Exit 
            P_Exit_Avg = P_Exit_Avg + Pr_X
            T_Exit_Avg = T_Exit_Avg + Temp_X
            ThrustX_Avg = ThrustX_Avg + ThrustX
            Correction_Factor_Avg = Correction_Factor_Avg + Correction_Factor
            
            Y_r2[p] = { chr(65+ch)+ "-" +str(count) : { "x" : C[i][0],
                                   "y" : C[i][1],
                                   "Theta": Theta0 ,
                                   
                                   }
                     }
            Y_r3[p] = { chr(65+ch)+ "-"+ str(count): { "x" : abs(C[i][0]),
                                   "y" : abs(C[i][1]),
                                   "Q" : Q,
                                   "R" : R,
                                   "Theta": ( Q - R )/2,
                                   "Beta": ( Q + R )/2,
                                   "v_PM": v_PMA,
                                   "Mach":M2_No,
                                   "Pr_Exit ": Pr_X,
                                   "Temp Exit":Temp_X,
                                   "V_Exit ": V_Exit ,
                                   "Epsilon=A_e/A_t" : A_Ratio,
                                   "Exit_Area":A_Exit,
                                   "M_Rate ": M_Rate,
                                   "Thrust " :ThrustX,
                                   }
                     }
            
            #pprint(str(C[i][0]) + "," + str(C[i][1]))
            plt.plot(C[i][0],C[i][1],'yo')
            p = p + 1 
            v =1
                    
        elif (i > 0 and len(div) > 0 and (Theta0 + (v * ThetaM)) <= Angle ):
            count = count + 1
            if ( v_PMA <= 0 or Beta <= 0  ):
                v_PMA = Beta
                Beta = v_PMA
                           
            if M2_No is None :    
                M2_No = float(M_Sonic)
                
            v_PMA = v_PM(G,M2_No)    
            Theta0=Angle - math.floor(Angle)
            ThetaM= math.ceil(math.floor(Angle) / (div[v]+1))
            Q = v_PMA + (Theta0 + (v * ThetaM))
            R= v_PMA - (Theta0 + (v * ThetaM))
            Theta =( Q - R ) / 2
            Beta = ( Q + R ) / 2
            if( Theta <= 0 ):
                Theta = (Theta0 + (v * ThetaM))
                
            M2_No = abs(M_No(G,M2_No,Beta,Theta))
            if(M2_No >=  M_Exit):
                M2_No =   2 * M_Exit
                
                
            
            Pr_X=Pr_Ratio(G,M2_No,Beta) * Pr_Throat
            D_X = D_Ratio(G,M2_No,Beta) * D_Throat
            Temp_X=T_Ratio(G,M2_No,Beta) * T_Throat
            V_Exit = Exit_Velocity(G,M2_No,Temp_X,R_Gas)
            A_Ratio = Area_Ratio(G , M2_No)
            if(A_Ratio > A_Ratio_MAX):
                A_Ratio = A_Ratio_MAX
                         
            # Parametric area for flow through conic section of shock wave ,
            # Area = 2 x Pi x r x dr
            #2 and Pi are constants of integration so integration of r.dr
            # gives , Flow Area = 2 . Pi . [r^2 /2 ] = Pi . r ^ 2
            
            A_Exit = A_Ratio * A_Throat
            M_Rate = Mass_Rate(G,Pr_X,Temp_X,R_Gas,A_Throat)
            Correction_Factor = Cf(G,Pc,Pr_X , P_Atm, R_Throat/1000 , R_Exit/1000 )
            ThrustX=( Exit_Thrust(M_Rate,V_Exit,Pr_X,P_Atm,A_Exit))

            Mach_Avg=(Mach_Avg + M2_No)
            M_Rate_Avg = M_Rate_Avg + M_Rate
            V_Exit_Avg = V_Exit_Avg + V_Exit 
            P_Exit_Avg = P_Exit_Avg + Pr_X
            T_Exit_Avg = T_Exit_Avg + Temp_X
            ThrustX_Avg = ThrustX_Avg + ThrustX
            Correction_Factor_Avg = Correction_Factor_Avg + Correction_Factor
            
            
            #M2_No =abs(M2_No)
            
            Y_r2[p] = { chr(65+ch)+ "-" +str(count) : { "x" : C[i][0],
                                   "y" : C[i][1],
                                   "Theta": Theta0 + (v * ThetaM),
                                                        
                                   }
                     }
            Y_r3[p] = { chr(65+ch)+ "-"+ str(count): { "x" : abs(C[i][0]),
                                   "y" : abs(C[i][1]),
                                   "Q" : Q,
                                   "R" : R,
                                   "Theta": ( Q - R )/2,
                                   "Beta": ( Q + R )/2,
                                   "v_PM": v_PMA,
                                   "Mach":M2_No,
                                   "Pr_Exit ": Pr_X,
                                   "Temp Exit":Temp_X,
                                   "V_Exit ": V_Exit ,
                                   "Epsilon=A_e/A_t" : A_Ratio,
                                   "Exit_Area":A_Exit,
                                   "M_Rate ": M_Rate,
                                   "Thrust " :ThrustX,
                                   }
                     }
            
            #pprint(str(C[i][0]) + "," + str(C[i][1]))
            plt.plot(C[i][0],C[i][1],'yo')             
            v = v + 1                 
            p = p + 1
        else :
            count = count + 1
            if ( v_PMA <= 0 or Beta <= 0  ):
                v_PMA = Beta
                Beta = v_PMA
                           
            if M2_No is None :    
                M2_No = float(M_Sonic)
                
            v_PMA = v_PM(G,M2_No)    
            Theta0=Angle - math.floor(Angle)
            ThetaM= math.ceil(math.floor(Angle) / (div[v]+1))
            Q = v_PMA + (Theta0 + (v * ThetaM))
            R= v_PMA - (Theta0 + (v * ThetaM))
            Theta =( Q - R ) / 2
            Beta = ( Q + R ) / 2
            if( Theta <= 0 ):
                Theta = (Theta0 + (v * ThetaM))
                
            M2_No = abs(M_No(G,M2_No,Beta,Theta))
            if(M2_No >=  M_Exit):
                M2_No =   2 * M_Exit
                
                
            
            Pr_X=Pr_Ratio(G,M2_No,Beta) * Pr_Throat
            D_X = D_Ratio(G,M2_No,Beta) * D_Throat
            Temp_X=T_Ratio(G,M2_No,Beta) * T_Throat
            V_Exit = Exit_Velocity(G,M2_No,Temp_X,R_Gas)
            A_Ratio = Area_Ratio(G , M2_No)
            if(A_Ratio > A_Ratio_MAX):
                A_Ratio = A_Ratio_MAX
                
            # Parametric area for flow through conic section of shock wave ,
            # Area = 2 x Pi x r x dr
            #2 and Pi are constants of integration so integration of r.dr
            # gives , Flow Area = 2 . Pi . [r^2 /2 ] = Pi . r ^ 2
            
            A_Exit = A_Ratio * A_Throat
            M_Rate = Mass_Rate(G,Pr_X,Temp_X,R_Gas,A_Throat)
            Correction_Factor = Cf(G,Pc,Pr_X , P_Atm, R_Throat/1000 , R_Exit/1000 )
            ThrustX=( Exit_Thrust(M_Rate,V_Exit,Pr_X,P_Atm,A_Exit))

            Mach_Avg=(Mach_Avg + M2_No)
            M_Rate_Avg = M_Rate_Avg + M_Rate
            V_Exit_Avg = V_Exit_Avg + V_Exit 
            P_Exit_Avg = P_Exit_Avg + Pr_X
            T_Exit_Avg = T_Exit_Avg + Temp_X
            ThrustX_Avg = ThrustX_Avg + ThrustX
            Correction_Factor_Avg = Correction_Factor_Avg + Correction_Factor
            
                        
            #M2_No =abs(M2_No)
            
            Y_r2[p] = { chr(65+ch)+ "-" +str(count) : { "x" : C[i][0],
                                   "y" : C[i][1],
                                   "Theta": Theta0 + (v * ThetaM),
                                                        
                                   }
                     }
            Y_r3[p] = { chr(65+ch)+ "-"+ str(count): { "x" : abs(C[i][0]),
                                   "y" : abs(C[i][1]),
                                   "Q" : Q,
                                   "R" : R,
                                   "Theta": ( Q - R )/2,
                                   "Beta": ( Q + R )/2,
                                   "v_PM": v_PMA,
                                   "Mach":M2_No,
                                   "Pr_Exit ": Pr_X,
                                   "Temp Exit":Temp_X,
                                   "V_Exit ": V_Exit ,
                                   "Epsilon=A_e/A_t" : A_Ratio,
                                   "Exit_Area":A_Exit,
                                   "M_Rate ": M_Rate,
                                   "Thrust " :ThrustX,
                                   }
                     }
            #pprint(str(C[i][0]) + "," + str(C[i][1]))
            plt.plot(C[i][0],C[i][1],'yo')             
            v = v + 1                 
            p = p + 1
            
        
        

print(p)
print(len(C))
Mach_No_Array=(Mach_Avg/(1 * (len(C)-1)))
#Mach_Avg = 0
M_Rate_Array = (M_Rate_Avg/(1 * (len(C)-1)))
#M_rate_Avg = 0
V_Exit_Array = (V_Exit_Avg /(1 * (len(C)-1)))
P_Exit_Array = (P_Exit_Avg/(1 * (len(C)-1)))
#P_Exit_Avg = 0
T_Exit_Array =(T_Exit_Avg/(1 * (len(C)-1)))
#T_Exit_Avg = 0
ThrustX_Array =abs(ThrustX_Avg/(1 * (len(C)-1)))
#ThrustX_Avg = 0
Correction_Factor_Array = (Correction_Factor_Avg/(1 * (len(C)-1)))

pprint("Correction Factor M_No and A_Ratio : " + str(Correction_Factor))
pprint("Average Correction Factor M_No and A_Ratio : " + str(Correction_Factor_Array))
pprint("Average Mach Number : " + str(Mach_No_Array))
pprint("Area Ratio [MAX] : " + str(A_Ratio_MAX))
pprint("Exit Radius [MIN] : " + str(R_Exit))
pprint("Exit Length [MIN] : " + str(L_Min))
pprint("Average Mass FlowRate (Kgm/sec): " + str(M_Rate_Array))
pprint("Average Exit Velocity (m/sec): " + str(V_Exit_Array))
pprint("Average Exit Pressure (N/m2) : " + str(P_Exit_Array))
pprint("Average Exit Temperature (K)  : " + str(T_Exit_Array))
pprint("Average Exit Thrust (N) : " + str(ThrustX_Array))


#pprint(Y_r2)
#pprint(Y_r3)

ch=0
count = 0
print(len(Y_r3),"no of records")
jsonElements = []               
array = json.dumps(Y_r3)
x = json.loads(array)
with open('./test_website_valuesX3.csv', 'w+', newline='') as outfile:
    f = csv.writer(outfile)
    # Write CSV Header, If you dont need that, remove this line
    f.writerow(["x" ,"y","Q","R","Theta","Beta","v_PM", "Mach", "Pr_exit" , "T_exit","V_exit", "A_Ratio","Exit_Area","M_rate", "Thrust"])

#pprint("x" + "  " + "y" + " " + "Q" + " " + "R" + "  " + "Theta" + " " + "v_PM" + " " + "Mach" + "  " + "Pr_exit" + " "+ "T_exit" + " " + "V_exit" + "  " + "A_Ratio" + " "+ "M_rate" + " " + "Thrust")
    for i in range(0,len(Y_r3),1):
        #pprint(Y_r3[i])
        for key,value in Y_r3[i].items():
            for x,value1 in value.items():
                #f.writerow(x["x"]),x["y"],x["Q"],x["R"])
                jsonElements.append(value1)
                #print(x,value1,"\n")
            f.writerow(jsonElements)
            jsonElements = []                    
plt.axis('equal')
plt.show()

