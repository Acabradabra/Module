#!/usr/bin/env python3
#-*- coding: utf-8 -*-
from numpy import *
from collections import OrderedDict
# from scipy.integrate import simps
from scipy.interpolate import Rbf
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.optimize    import minimize
from scipy.optimize    import minimize_scalar
#from PIL import Image
# from math import *
# import Utilities as util
# from scipy.interpolate import griddata
# import Porosity as po
import Utilities as util
# import matplotlib.transforms as tr
import os
import sys
# import time
# import h5py as h5
# try : import h5py as h5
# except : 
	# os.system('echo Load h5py ; source ~/moduleH5PY.sh')
	# import h5py as h5

(plt,mtp)=util.Plot0()

#======> DATA_MANIP
DATA_MANIP='/Users/fmuller/Soft/DATA-MANIP/'

#======> Water Flow meter
Qfmet=[1.8,32]
Vfmet=[2,10]

#======> Pyrometer
Tpyro=[250,1200]
Vpyro=[0.8,4]
# Vpyro=[0.5,4]

#======> HW0
# Phw=[  -8.08864063,   50.13406312, -111.688949,    108.0523443,   -38.70044368] #===> All0 , U : 0.02,0.5 , E : 1.42,1.57
# Chw=[   1.31373307,    1.62208518                                             ] #===> All0 , U : 0.02,0.5 , E : 1.42,1.57
# Phw=[  37.42890761, -210.34514692,  499.43946881, -555.37747293,  233.89735974] #===> All1 , U : 16  ,60  , E : 1.68,2.1§
# Chw=[   1.35793283,    0.43463635                                             ] #===> All1 , U : 16  ,60  , E : 1.68,2.1§

#======> HW1
# Phw=[ 0.92821247, -2.33414951,  0.87927599,  1.79213319, -1.31992371] #===> Noze0
# Chw=[ 1.41652803,  2.01172194                                       ] #===> Noze0
# Phw=[  37.13967391, -218.51417678,  539.85684914, -623.40018703,  272.75803063] #===> Noze 1    , U : 0.5,58.5 , E : ,2.2
# Chw=[   1.47220749,    0.45002202                                             ] #===> Noze 1    , U : 0.5,58.5 , E : ,2.2
# Phw=[  22.31391098, -137.27233784,  353.02101289, -422.2841403,   190.96603416] #===> Noze 1 UP , U : 0.5,60   , E : 1.37,2.44
# Chw=[   1.51889546,    0.58906092                                             ] #===> Noze 1 UP , U : 0.5,60   , E : 1.37,2.44
# Phw=[ -34.78849233,  228.57351992, -532.6906731 ,  536.12538639, -199.33624338] #===> Noze 1 UP 2 , U : 0.27,2 , E : 1.32,1.53
# Chw=[   1.37157144,    0.6901588                                              ] #===> Noze 1 UP 2 , U : 0.27,2 , E : 1.32,1.53
# Phw=[  9.65743756,  -42.34049485,  87.64569601, -97.01439692,  43.78526429] #===> Noze 1 UP 3 , U : 0.5,57.6 , E : 1.37,2,45
# Chw=[  1.5221189 ,    0.60624963                                          ] #===> Noze 1 UP 3 , U : 0.5,57.6 , E : 1.37,2,45
# Phw=[  -0.14341575,    3.03305622, -9.42878159, 10.71677025, -4.23407232] #===> Noze 0 UP , U : 0.02,0.5 , E : 1.32,1.73
# Chw=[   1.4074094 ,    2.25114063                                       ] #===> Noze 0 UP , U : 0.02,0.5 , E : 1.32,1.73

#======> HW2
# Phw=[ -4613.7064021  ,  25988.92770468 , -54875.07226792 , 51479.21605271 , -18105.10210283 ] #===> Noze 1
# Chw=[     1.58733931 ,      0.72729577                                                      ] #===> Noze 1

#======> HW3
Phw=[ -641.66057336 , 3997.38563227 , -9321.61969109 , 9649.2191077 , -3742.59324519] #===> Noze 1 , U : 0.178,1.018 , E : 
Chw=[    1.76203906 ,    0.82538195                                                 ] #===> Noze 1 , U : 0.178,1.018 , E : 
# Phw=[ -18.2504436 ,  129.29429818, -323.49397333,  347.93777799, -137.69650722] #===> Noze 1 , U : 0.332,2.021
# Chw=[   1.75497399,    0.87368903                                             ] #===> Noze 1 , U : 0.332,2.021

#======> PCI
PCI_CH4=50e6  # J/Kg
PCI_H2=120e6  # J/kg
#======> Normal
P0=101325 # Pa
T0=273.15 # K
#======> Ambiant
Pamb=101325 # Pa
Tamb=298.15 # K
#======> Universal gas constant
R=8.314 # J/mol K
#======> Air
XO2_a=0.21
aN2=(1-XO2_a)/XO2_a
#======> Molar mass
WCH4=(12+4*1 )*1e-3 # Kg/mol
WH2O=(2*1+16 )*1e-3 # Kg/mol
WH2 =(2*1    )*1e-3 # Kg/mol
WO2 =(2*16   )*1e-3 # Kg/mol
WN2 =(2*14   )*1e-3 # Kg/mol
WNO =(14+16  )*1e-3 # Kg/mol
WNO2=(14+2*16)*1e-3 # Kg/mol
Wa  =0.21*WO2+0.79*WN2

Rho0_a  =(Wa  *P0  )/(R*T0  ) # Kg/m3
Rho0_H2 =(WH2 *P0  )/(R*T0  ) # Kg/m3
Rho0_CH4=(WCH4*P0  )/(R*T0  ) # Kg/m3

RhoA_a  =(Wa  *Pamb)/(R*Tamb)
RhoA_H2 =(WH2 *Pamb)/(R*Tamb)
RhoA_CH4=(WCH4*Pamb)/(R*Tamb)

#======> Planck
kb=1.38064900e-23 # Boltzman cst
hp=6.62607015e-34 # Planck cst
cl=299792458      # Celerity
si=5.670374e-8    # Steph Boltzman cst

#---------------------------------------------------------------------
######################################             Uncertainties     #
#---------------------------------------------------------------------

#=====> Temperature
Ea_T = 1

#=====> Flow rates
Ea_Qa  =0.003*130
Ea_QH2 =0.003*20
Ea_QCH4=0.001*19.05
Er_Qa  =0.005
Er_QH2 =0.005
Er_QCH4=0.005

Ea_h   =1

#=====> Divers
Ea_r=1e-3 # Radius

#---------------------------------------------------------------------
######################################             Rotameter         #
#---------------------------------------------------------------------

H_GT1000=[20  ,250   ]
Q_GT1000=[0.41,410/60]
(VQH_GT1000,EP,detP)=util.InterP(Q_GT1000,H_GT1000,2,len(H_GT1000))
(VHQ_GT1000,EP,detP)=util.InterP(H_GT1000,Q_GT1000,2,len(Q_GT1000))

Conv_CH4=1.34

Ea_QH2O=Ea_h*VQH_GT1000[1]

# mm , nL/min
Calib_Rot1355=array([ # 1355
[ 20.81, 6.02],
[ 22.25, 6.63],
[ 24.05, 7.23],
[ 25.86, 7.84],
[ 27.48, 8.44],
[ 29.46, 9.05],
[ 30.91, 9.65],
[ 32.71,10.26],
[ 34.69,10.86],
[ 35.95,11.47],
[ 37.76,12.07],
[ 39.38,12.68],
[ 41.00,13.28],
[ 42.80,13.89],
[ 44.25,14.49],
[ 46.05,15.10],
[ 47.49,15.70],
[ 49.11,16.31],
[ 50.02,16.67],
[ 50.56,16.91],
[ 51.82,17.27],
[ 52.54,17.52],
[ 53.26,17.88],
[ 53.98,18.12],
[ 55.06,18.48],
[ 55.42,18.73],
[ 56.51,19.09],
[ 57.23,19.33],
[ 57.95,19.69],
[ 58.67,19.94],
[ 59.75,20.30],
[ 60.47,20.54],
[ 61.19,20.90],
[ 61.91,21.15],
[ 62.82,21.51],
[ 63.54,21.75],
[ 64.44,22.11],
[ 64.98,22.36],
[ 65.88,22.72],
[ 66.60,22.96],
[ 67.50,23.32],
[ 68.04,23.57],
[ 68.94,23.93],
[ 69.49,24.17],
[ 70.75,24.53],
[ 71.11,24.78],
[ 72.01,25.14],
[ 72.73,25.38],
[ 73.63,25.74],
[ 74.17,25.99],
[ 74.89,26.35],
[ 75.62,26.59],
[ 76.70,26.95],
[ 77.42,27.20],
[ 78.14,27.56],
[ 78.68,27.80],
[ 79.76,28.16],
[ 80.30,28.41],
[ 81.20,28.77],
[ 81.74,29.01],
[ 82.65,29.37],
[ 83.19,29.62],
[ 84.27,29.98],
[ 84.63,30.22],
[ 85.71,30.58],
[ 86.07,30.83],
[ 86.97,31.19],
[ 87.69,31.43],
[ 88.60,31.79],
[ 89.14,32.04],
[ 90.22,32.40],
[ 90.76,32.64],
[ 91.48,33.00],
[ 92.20,33.25],
[ 92.92,33.61],
[ 93.64,33.85],
[ 94.36,34.21],
[ 95.09,34.46],
[ 95.99,34.82],
[ 96.53,35.06],
[ 97.61,35.42],
[ 97.97,35.67],
[ 99.05,36.03],
[ 99.41,36.27],
[100.31,36.63],
[100.85,36.88],
[101.94,37.24],
[102.30,37.48],
[103.20,37.84],
[103.74,38.09],
[104.64,38.45],
[105.18,38.69],
[106.26,39.05],
[106.62,39.30],
[107.53,39.66],
[108.07,39.90],
[108.97,40.26],
[109.51,40.51],
[110.59,40.87],
[110.95,41.11],
[112.03,41.47],
[112.39,41.72],
[113.29,42.08],
[113.84,42.32],
[114.74,42.68],
[115.28,42.93],
[116.18,43.29],
[116.72,43.53],
[117.44,43.89],
[118.16,44.14],
[119.06,44.50],
[119.60,44.74],
[120.33,45.10],
[121.05,45.35],
[121.77,45.71],
[122.49,45.95],
[123.21,46.31],
[123.75,46.56],
[124.65,46.92],
[125.01,47.16],
[125.91,47.52],
[126.46,47.77],
[127.36,48.13],
[127.72,48.37],
[128.44,48.73],
[128.98,48.98],
[129.70,49.34],
[130.42,49.58],
[131.32,49.94],
[131.68,50.19],
[132.40,50.55],
[132.95,50.79],
[133.67,51.15],
[134.03,51.40],
[134.93,51.76],
[135.47,52.00],
[136.19,52.36],
[136.55,52.61],
[137.45,52.97],
[137.99,53.21],
[138.71,53.57],
[138.89,53.82],
[139.80,54.18],
[140.52,54.42],
[141.24,54.78],
[141.78,55.03],
[142.50,55.39],
[142.86,55.63],
[143.76,55.99],
[144.12,56.24],
[144.84,56.60],
[145.38,56.84],
[146.29,57.20],
[146.65,57.45],
[147.37,57.81],
[147.91,58.05],
[148.63,58.41],
[149.17,58.66],
[149.71,59.02],
[150.00,60.00]
])

Calib_Rot1357=array([ # 1357
# [ 14.50389,0.5307575],
[ 14.50389,0.542212333],
[ 15.79927,0.609259667],
[ 16.46228,0.62423],
[ 17.30537,0.649339167],
[ 18.32408,0.695995333],
[ 19.07125,0.731570833],
[ 20.69268,0.806869667],
[ 22.38585,0.878755333],
[ 24.37279,0.9583025],
[ 26.40915,1.033007667],
[ 27.98662,1.111987833],
[ 29.66278,1.18726],
[ 31.62298,1.2632125],
[ 33.54046,1.332942833],
[ 35.01937,1.405197833],
[ 36.61702,1.4728595],
[ 38.17481,1.545107333],
[ 40.14896,1.614642],
[ 42.21812,1.701085167],
[ 44.04663,1.7764855],
[ 45.70515,1.852008333],
[ 47.86682,1.930382167],
[ 50.15892,2.017749333],
[ 51.94168,2.102126333],
[ 53.80318,2.179249],
[ 55.76186,2.253214833],
[ 57.40508,2.334285],
[ 59.79926,2.425679333],
[ 61.36479,2.505921167],
[ 63.14754,2.5684665],
[ 65.18497,2.661226833],
[ 67.2224 ,2.738353],
[ 69.25983,2.8326185],
[ 71.29726,2.909322167],
[ 73.33852,2.995318167],
[ 75.43942,3.081331833],
[ 77.18285,3.155858333],
[ 78.87337,3.222214333],
[ 80.53154,3.2953755],
[ 82.50313,3.3718265],
[ 84.79524,3.471970167],
[ 86.6857 ,3.535919667],
[ 88.61542,3.618157667],
[ 90.14912,3.689576],
[ 92.18092,3.766103167],
[ 94.21835,3.851813333],
[ 95.84415,3.919065833],
[ 97.62689,3.993846667],
[ 99.82129,4.083776333],
[101.64607,4.160902333],
[103.89015,4.252929],
[105.89737,4.331583667],
[107.67671,4.397913667],
[109.70993,4.479186],
[111.51929,4.560303167],
[113.52802,4.646644167],
[115.06641,4.717925],
[117.13943,4.791102167],
[119.18216,4.880651833],
[121.2143 ,4.964196333],
[123.15947,5.043569],
[124.69703,5.111998],
[126.76541,5.190915],
[128.8744 ,5.2758235],
[130.56816,5.3434235],
[132.28183,5.420629833],
[134.34547,5.504735833],
[136.34352,5.592021667],
[138.68013,5.677459667],
[140.48904,5.762918167],
[142.77695,5.8479695],
[144.62372,5.926235],
[146.86172,6.013979167],
[148.69643,6.087800667],
[150.23929,6.157613667],
[152.28511,6.234739667],
[154.06786,6.311134333],
[155.80145,6.3804225],
[157.83459,6.4557765],
[159.71498,6.534675167],
[161.56396,6.610013333],
[163.40714,6.688927167],
[165.39069,6.774622667],
[167.22266,6.83461],
[169.34857,6.916718],
[171.60294,7.014571167],
[173.51056,7.093102],
[175.12267,7.1616485],
[177.26888,7.241688],
[179.15608,7.323075833],
[181.29527,7.408536833],
[183.3942 ,7.484057167],
[185.41928,7.572216833],
[187.40939,7.652927667],
[189.40097,7.73004],
[191.42467,7.812608333],
[193.03371,7.8800985],
[195.32581,7.967410333],
[197.02968,8.0429205],
[198.85922,8.111476667],
[200.67407,8.184992],
[202.7115 ,8.269827833],
[204.79372,8.358624333],
[206.80179,8.440377833],
[208.82379,8.515244],
[210.86124,8.605381333],
[213.40802,8.702876833],
[215.44545,8.7865295],
[217.07695,8.852972333],
[218.92986,8.9258505],
[220.68141,8.990195667],
[222.26529,9.0590555],
[224.03291,9.138584333],
[225.63664,9.205330333],
[227.67003,9.2822825],
[229.79936,9.374270667],
[231.74489,9.443512],
[233.8423 ,9.525456667],
[235.81976,9.611152167],
[238.11185,9.688621],
[239.93117,9.770178167],
[241.77457,9.8418955],
[243.89752,9.927750667],
[246.0069 ,10.0013815],
[247.61794,10.07201033]
])

AirData=array([
# T 	ρ 	    μ 	      ν   	      Cp 	    λ 	      a 	    Pr 
# K 	kg m−3 	kg m−1 s−1 	m2 s−1 	J kg−1 K−1 	W m−1 K−1 m2 s−1 	- 
[  250 , 1.413 , 1.60e-5 ,  0.949e-5 , 1005 , 0.0223 ,  1.32e-5 , 0.722 ],
[  300 , 1.177 , 1.85e-5 ,  1.570e-5 , 1006 , 0.0262 ,  2.22e-5 , 0.708 ],
[  350 , 0.998 , 2.08e-5 ,  2.080e-5 , 1009 , 0.0300 ,  2.98e-5 , 0.697 ],
[  400 , 0.883 , 2.29e-5 ,  2.590e-5 , 1014 , 0.0337 ,  3.76e-5 , 0.689 ],
[  450 , 0.783 , 2.48e-5 ,  2.890e-5 , 1021 , 0.0371 ,  4.22e-5 , 0.683 ],
[  500 , 0.705 , 2.67e-5 ,  3.690e-5 , 1030 , 0.0404 ,  5.57e-5 , 0.680 ],
[  550 , 0.642 , 2.85e-5 ,  4.430e-5 , 1039 , 0.0436 ,  6.53e-5 , 0.680 ],
[  600 , 0.588 , 3.02e-5 ,  5.130e-5 , 1055 , 0.0466 ,  7.51e-5 , 0.680 ],
[  650 , 0.543 , 3.18e-5 ,  5.850e-5 , 1063 , 0.0495 ,  8.58e-5 , 0.682 ],
[  700 , 0.503 , 3.33e-5 ,  6.630e-5 , 1075 , 0.0523 ,  9.67e-5 , 0.684 ],
[  750 , 0.471 , 3.48e-5 ,  7.390e-5 , 1086 , 0.0551 , 10.80e-5 , 0.686 ],
[  800 , 0.441 , 3.63e-5 ,  8.230e-5 , 1098 , 0.0578 , 12.00e-5 , 0.689 ],
[  850 , 0.415 , 3.77e-5 ,  9.070e-5 , 1110 , 0.0603 , 13.10e-5 , 0.692 ],
[  900 , 0.392 , 3.90e-5 ,  9.930e-5 , 1121 , 0.0628 , 14.30e-5 , 0.696 ],
[  950 , 0.372 , 4.02e-5 , 10.800e-5 , 1132 , 0.0653 , 15.50e-5 , 0.699 ],
[ 1000 , 0.352 , 4.15e-5 , 11.800e-5 , 1142 , 0.0675 , 16.80e-5 , 0.702 ],
[ 1100 , 0.320 , 4.40e-5 , 13.700e-5 , 1161 , 0.0723 , 19.50e-5 , 0.706 ],
[ 1200 , 0.295 , 4.63e-5 , 15.700e-5 , 1179 , 0.0763 , 22.00e-5 , 0.714 ],
[ 1300 , 0.271 , 4.85e-5 , 17.900e-5 , 1197 , 0.0803 , 24.80e-5 , 0.722 ]
])

YounisData=array([
# dp     C       m
[ 1.52*1e-3 , 0.146 , 0.96 ],
[ 0.94*1e-3 , 0.139 , 0.92 ],
[ 0.76*1e-3 , 0.456 , 0.70 ],
[ 0.42*1e-3 , 0.485 , 0.55 ],
[ 0.29*1e-3 , 0.638 , 0.42 ]
])
fY_C=interp1d( YounisData[:,0] , YounisData[:,1] , kind='cubic' )
fY_m=interp1d( YounisData[:,0] , YounisData[:,2] , kind='cubic' )

#---------------------------------------------------------------------
######################################             Radiation         #
def Planck(l,T) : return( (2*hp*cl**2/l**5)/( exp(hp*cl/(l*kb*T))-1 ) )
#---------------------------------------------------------------------
def ConvPor( por_s ) : return( 25.4e-3/por_s )
#---------------------------------------------------------------------
######################################             Function          #
#---------------------------------------------------------------------
def circle(X0,x0,r) : return(sqrt(r**2-(X0-x0)**2)) # y=f(x)
def Circle(P0,r,N) :
    (x0,y0)=P0
    X0=linspace(x0-r,x0+r,N)
    Y0=circle(X0,x0,r)
    return( X0,y0+Y0,y0-Y0,0*X0 ) # x , ytop , ybot , z
#####################################################################
def Correlation(T0,T1) :
    MeanT0=mean(T0) ; T0p=T0-MeanT0 ; StdT0=sqrt( mean(T0p**2) )
    MeanT1=mean(T1) ; T1p=T1-MeanT1 ; StdT1=sqrt( mean(T1p**2) )
    # Corel=sum(T0p*T1p)/( sum(abs(T0p))*sum(abs(T1p)) )
    Corel=mean(T0p*T1p)/( StdT0*StdT1 )
    return( [MeanT0,T0p,StdT0],[MeanT1,T1p,StdT1],Corel )
#####################################################################
def CheckStability(MT) : 
    MaxMT=MT.max(axis=0)
    MinMT=MT.min(axis=0)
    DT   =        MaxMT[1:12]-MinMT[1:12] ; DTmax=max(DT) ; idmax=list(DT).index(DTmax)
    if   all( abs(MaxMT[1:12]-MinMT[1:12])<5    ) : return( 'std' , MinMT , MaxMT )
    elif all( abs(MT[-1,1:12]-MaxMT[1:12])<5e-2 ) : return( 'inc' , MinMT , MaxMT )
    elif all( abs(MT[-1,1:12]-MinMT[1:12])<5e-2 ) : return( 'dec' , MinMT , MaxMT )
    else : print('DTmax : {0:.2f} , Idmax : {1:.0f}'.format(DTmax,idmax)) ; return( 'not' , MinMT , MaxMT )
#####################################################################
def ColState(state) :
    col='r'*(state=='dec')+'b'*(state=='std')+'g'*(state=='inc')+'y'*(state=='not')
    # util.Section(state,0,1,col)
    return(col)
#####################################################################
def DebitFinder(name) :
    if   'KW' in name : return( float(name.split('-')[1]) )
    elif 'Qa' in name : return( float(name[3:])           )
    else              : return( float(name.split('-')[0]) )
#####################################################################
def PicNameConverter(name) :
    if   '(' in name : return(name.split(' ')[0  ])
    elif '_' in name : return(name.split('_')[0  ])
    else             : return(name           [:-7])
#---------------------------------------------------------------------
######################################             Error             #
#---------------------------------------------------------------------
def DaS(r) : return( 2*pi*r*Ea_r )
#---------------------------------------------------------------------
def DaQa  (Qa         ) : return( max( Ea_Qa   , Er_Qa  *Qa   )  )
def DaQH2 (   QH2     ) : return( max( Ea_QH2  , Er_QH2 *QH2  )  )
def DaQCH4(       QCH4) : return( max( Ea_QCH4 , Er_QCH4*QCH4 )  )
def DaQtot(Qa,QH2,QCH4) : return( DaQa(Qa)+DaQH2(QH2)+DaQCH4(QCH4) )
# def DaQtot(Qa,QH2,QCH4) : return( DaQa(Qa)+DaQH2(QH2)+DaQCH4(QCH4) )
#---------------------------------------------------------------------
def DaPhi_CH4(Qa    ,QCH4) : return( Phi_CH4(Qa    ,QCH4)*(DaQa(Qa)/Qa+DaQCH4(QCH4)/QCH4) )
def DaPhi_H2 (Qa,QH2     ) : return( Phi_H2 (Qa,QH2     )*(DaQa(Qa)/Qa+DaQH2 (QH2 )/QH2 ) )
def DaPhi_tot(Qa,QH2,QCH4) : 
	out=0
	if QH2 >0 : out+=DaPhi_H2 (Qa,QH2)
	if QCH4>0 : out+=DaPhi_CH4(Qa,QCH4)
	return( out )
#---------------------------------------------------------------------
def DrU( Qa,QH2,QCH4,r ) : return( DaS(r)/(2*pi*r**2) + DaQtot(Qa,QH2,QCH4)/(Qa+QH2+QCH4) )
#---------------------------------------------------------------------
######################################             Debit             #
#---------------------------------------------------------------------
def Vit(Q,D,Tu) : r=0.5*D ; S=pi*r**2 ; Q2=Q/6e4 ; Q3=Q2*(Tu/T0) ; return( Q3/S )
#---------------------------------------------------------------------
def Deb(V,D,Tu) : r=0.5*D ; S=pi*r**2 ; Q0=V*S   ; Q1=Q0*6e4 ; Q2=Q1*T0/Tu ; return(Q2)
#---------------------------------------------------------------------
def Alpha(hyb) :
    if   hyb==0 : return(0)
    elif hyb==1 : return(1)
    else        :
        rh=hyb/(1-hyb)
        rq=rh*(Rho0_CH4*PCI_CH4)/(Rho0_H2 *PCI_H2 )
        return( rq/(1+rq) )
#---------------------------------------------------------------------
def Q_H2 ( phi,Qa)       : return( phi*Qa*0.21/0.5 )
def Q_CH4( phi,Qa)       : return( phi*Qa*0.21/2   )
def Qa_tot(phi,QH2,QCH4) : return( (0.5*QH2+2*QCH4)/(0.21*phi) )
#---------------------------------------------------------------------
def Phi_H2 (Qa,Qf)       : 
	if all(Qa>0) : return( (0.5/0.21)*(Qf/Qa) )
	else    : return( 0 )
def Phi_CH4(Qa,Qf)       : 
	if all(Qa>0) : return( (2  /0.21)*(Qf/Qa) )
	else    : return( 0 )
def Phi_tot(Qa,QH2,QCH4) : return( Phi_H2(Qa,QH2)+Phi_CH4(Qa,QCH4) )
#---------------------------------------------------------------------
def MassFraction(Fuel,phi) :
	XO2air=0.21 ; a=3.76 #(1-XO2air)/XO2air
	sch4=2
	sh2 =0.5
	if   Fuel=='CH4' : s=sch4 ; Wf=WCH4
	elif Fuel=='H2'  : s=sh2  ; Wf=WH2
	else : sys.exit(2*'\n'+('\033[31m =====> Fuel  not build yet \033[0m'+2*'\n'))

	sy=s*WO2/Wf
	ay=a*WN2/WO2
	b=( 1+(1+ay)*sy/phi )

	# print('aN2 : ',a )
	# print('s   : ',sy)
	# print('a   : ',ay)
	# print('b   : ',b )

	Yf =1/b
	YO2=Yf*sy/phi
	YN2=YO2*a*WN2/WO2

	# print('Yf  : ',Yf )
	# print('YO2 : ',YO2)
	# print('YN2 : ',YN2)

	return( Yf,YO2,YN2 )
#---------------------------------------------------------------------
def Hyb_Power(p,a) :
	Pch4=(1-a)*p ; Qch4=(Pch4*1e3)/(Rho0_CH4*PCI_CH4) ; QCH4=Qch4*(60*1000) # nL/min 
	Ph2 =   a *p ; Qh2 =(Ph2 *1e3)/(Rho0_H2 *PCI_H2 ) ; QH2 =Qh2 *(60*1000) # nL/min 
	return(QCH4,QH2)
#---------------------------------------------------------------------
def Compo_CH4(a,phi) :
    #Reactifs
    NCH4_r=1
    NO2_r =2/phi
    NN2_r =2*a/phi
    Ntot_r=NCH4_r+NO2_r+NN2_r
    #Produits
    NCO2_p=1
    NH2O_p=2
    NO2_p =(1-phi)*2/phi
    NN2_p =2*a/phi
    # Ntot_p=NCO2_p+NH2O_p+NO2_p+NN2_p
    Ntot_p=NCO2_p+NO2_p+NN2_p
    # return( [NCH4_r/Ntot_r,NO2_r/Ntot_r,NN2_r/Ntot_r],[NCO2_p/Ntot_p,NH2O_p/Ntot_p,NO2_p/Ntot_p,NN2_p/Ntot_p] )
    return( [NCH4_r/Ntot_r,NO2_r/Ntot_r,NN2_r/Ntot_r],[NCO2_p/Ntot_p,NO2_p/Ntot_p,NN2_p/Ntot_p] )
#---------------------------------------------------------------------
######################################             HW                #
#---------------------------------------------------------------------
def VelocityProfile(d,Q0,D,ax1,ax2) :
    S=pi*(0.5*D)**2
    vit0=Vit(Q0,D)                   ; print('=> Wmoy : ',vit0)
    Ees=EtalonEU(vit0,Chw[0],Chw[1]) ; print('=> Emoy : ',Ees)    
    Ray=[]
    E_min,E_moy,E_max,E_dev=[],[],[],[]
    U_min,U_moy,U_max,U_dev=[],[],[],[]
    for r in os.popen('ls '+d+'/*.dat').read().split('\n')[:-1] :
    	if not 'Ref' in r :
	        Ray.append(float(r.split('/')[-1].split('.')[0][1:]))
	        
	        Ehw=ReadHw(r)
	        E_min.append(min (Ehw))
	        E_moy.append(mean(Ehw))
	        E_max.append(max (Ehw))
	        E_dev.append(std (Ehw))
            # E_dev.append(std (Ehw))

	        Uhw=polyval(Phw,Ehw)
	        # Uhw=EtalonUE(Ehw,Chw[0],Chw[1])
	        U_min.append(min (Uhw))
	        U_moy.append(mean(Uhw))
	        U_max.append(max (Uhw))
	        U_dev.append(std (Uhw))

    if Ray :
        Ray  =array(Ray)*1e-4 #0.5e-3
        E_min=array(E_min)
        E_moy=array(E_moy)
        E_max=array(E_max)
        U_min=array(U_min)
        U_moy=array(U_moy)
        U_max=array(U_max)

        print( '=> R  : {0:.4f}  ,  {1:.4f}  ,  {2:.4f}'.format(min(Ray  ),mean(Ray  ),max(Ray  )) )
        print( '=> E0 : {0:.4f}  ,  {1:.4f}  ,  {2:.4f}'.format(min(E_moy),mean(E_moy),max(E_moy)) )
        print( '=> U0 : {0:.4f}  ,  {1:.4f}  ,  {2:.4f}'.format(min(U_moy),mean(U_moy),max(U_moy)) )

        Sel=U_moy>0.7
        # Sel=Ray<10e-3
        # Ray2=abs(Ray-45e-3)
        RI,U0=Ray[Sel],U_moy[Sel] ; RI-=RI[0] ; DRI=RI[-1]-RI[0] 
        Ray2=abs(RI-0.5*RI[-1])
        # UR=U_moy*Ray2
        UI=U0*Ray2
        IntR=simps(UI,RI)
        Int0=simps(U0,RI)       ; vit1=Int0/DRI ; Ev1=abs(vit1-vit0)/vit0
        Q1=pi*IntR              ; vit2=Q1/S     ; Ev2=abs(vit2-vit0)/vit0
        Q2=Q1*6e4 ; Eq2=abs(Q2-Q0)/Q0
        print('=> DR : {0:.2e} m'.format( DRI ))
        print('=> Q0 : {0:.2f} L/min  ,  Q1 : {1:.2e} m3/s ,  Q2 : {2:.2f} L/min'.format(Q0,Q1,Q2) )
        print('=> V0 : {0:.2f} m/s    ,  V1 : {1:.2f} m/s  ,  V2 : {2:.2f} m/s  '.format(vit0,vit1,vit2) )
        print('=> Eq2 : {0:.2f} %     ,  Ev1 : {1:.2f} %   ,  Ev2 : {2:.2f} %   '.format(Eq2*100,Ev1*100,Ev2*100) )
        # print('=> Q0 : {0:.2e} m3/s  ,  Q2 : {1:.2f} L/min ,  V : {2:.2f} m/s  ,  Erreur : {3:.4f} %'.format(Q1,Q2,vit2,Eq2*100) )

        # ax1.plot( Ray*1e3 ,      E_min , ':'  +'r')
        ax1.plot( Ray*1e3 ,      E_moy , '-o' +'k')
        # ax1.plot( Ray*1e3 ,      E_max , ':'  +'r')
        ax1.plot( Ray*1e3 , mean(E_moy)+0*Ray, '--' +'g')
        ax11=ax1.twinx()
        ax11.plot( Ray*1e3 , E_dev/E_moy , ':'  +'b')
        # ax1.plot( Ray*1e3 , E_moy+E_dev , ':' +'b')
        # ax1.plot( Ray*1e3 , E_moy-E_dev , ':' +'b')
        # ax1.xlim(0,90)
        # ax1.ylim(0,0.5)
    
        # ax2.plot( Ray*1e3 , U_min , ':' +'r')
        ax2.plot( Ray*1e3 , U_moy , '-o'+'k') ; ax2.plot( Ray[Sel]*1e3 , U_moy[Sel] , 'x'+'r')
        # ax2.plot( Ray*1e3 , U_max , ':' +'r')
        ax2.plot( Ray*1e3 , mean(U_moy)+0*Ray, '--' +'g')
        ax22=ax2.twinx()
        ax22.plot( Ray*1e3 , U_dev/U_moy , ':' +'b')
        # ax2.plot( Ray*1e3 , U_moy+U_dev , ':' +'b')
        # ax2.plot( Ray*1e3 , U_moy-U_dev , ':' +'b')

        # ax2.plot( [Ray[0]*1e3,Ray[-1]*1e3] , 2*[vit0] , '--g' )
        # ax2.set_xlim(0,90)
        # ax2.set_ylim(0,0.4)
#---------------------------------------------------------------------
def Calib_HW( hw ) :
	if   hw=='HW0' :
		Phw=[  -8.08864063,   50.13406312, -111.688949,    108.0523443,   -38.70044368] #===> All0 , U : 0.02,0.5 , E : 1.42,1.57
		Chw=[   1.31373307,    1.62208518                                             ] #===> All0 , U : 0.02,0.5 , E : 1.42,1.57
		# Phw=[  37.42890761, -210.34514692,  499.43946881, -555.37747293,  233.89735974] #===> All1 , U : 16  ,60  , E : 1.68,2.1 
		# Chw=[   1.35793283,    0.43463635                                             ] #===> All1 , U : 16  ,60  , E : 1.68,2.1 
	elif hw=='HW1' :
		# Phw=[ 0.92821247, -2.33414951,  0.87927599,  1.79213319, -1.31992371] #===> Noze0
		# Chw=[ 1.41652803,  2.01172194                                       ] #===> Noze0
		# Phw=[  37.13967391, -218.51417678,  539.85684914, -623.40018703,  272.75803063] #===> Noze 1    , U : 0.5,58.5 , E : ,2.2
		# Chw=[   1.47220749,    0.45002202                                             ] #===> Noze 1    , U : 0.5,58.5 , E : ,2.2
		# Phw=[  22.31391098, -137.27233784,  353.02101289, -422.2841403,   190.96603416] #===> Noze 1 UP , U : 0.5,60   , E : 1.37,2.44
		# Chw=[   1.51889546,    0.58906092                                             ] #===> Noze 1 UP , U : 0.5,60   , E : 1.37,2.44
		Phw=[ -34.78849233,  228.57351992, -532.6906731 ,  536.12538639, -199.33624338] #===> Noze 1 UP 2 , U : 0.27,2 , E : 1.32,1.53
		Chw=[   1.37157144,    0.6901588                                              ] #===> Noze 1 UP 2 , U : 0.27,2 , E : 1.32,1.53
		# Phw=[  9.65743756,  -42.34049485,  87.64569601, -97.01439692,  43.78526429] #===> Noze 1 UP 3 , U : 0.5,57.6 , E : 1.37,2,45
		# Chw=[  1.5221189 ,    0.60624963                                          ] #===> Noze 1 UP 3 , U : 0.5,57.6 , E : 1.37,2,45
		# Phw=[  -0.14341575,    3.03305622, -9.42878159, 10.71677025, -4.23407232] #===> Noze 0 UP , U : 0.02,0.5 , E : 1.32,1.73
		# Chw=[   1.4074094 ,    2.25114063                                       ] #===> Noze 0 UP , U : 0.02,0.5 , E : 1.32,1.73
	elif hw=='HW2' :
		Phw=[ -4613.7064021  ,  25988.92770468 , -54875.07226792 , 51479.21605271 , -18105.10210283 ] #===> Noze 1 , U : 0.160,0.503 # 20/09/2021
		Chw=[     1.58733931 ,      0.72729577                                                      ] #===> Noze 1 , U : 0.160,0.503 # 20/09/2021
	elif hw=='HW3a' :
		Phw=[ -641.66057336 , 3997.38563227 , -9321.61969109 , 9649.2191077 , -3742.59324519] #===> Noze 1 , U : 0.178,1.018 # 22/09/2021
		Chw=[    1.76203906 ,    0.82538195                                                 ] #===> Noze 1 , U : 0.178,1.018 # 22/09/2021
	elif hw=='HW3b' :
		Phw=[ -18.2504436 ,  129.29429818, -323.49397333,  347.93777799, -137.69650722] #===> Noze 1 01 , U : 0.332,2.021 # 14/10/2021
		Chw=[   1.75497399,    0.87368903                                             ] #===> Noze 1 01 , U : 0.332,2.021 # 14/10/2021
	return( Phw,Chw )
#---------------------------------------------------------------------
def EtalonEU(U,A,B) : return(sqrt( A+B*sqrt(U) ))
def EtalonUE(E,A,B) : return( ((E**2-A)/B)**2   )
#---------------------------------------------------------------------
def ReadHw(f) : return( array([float(v) for v in os.popen("sed -e '1,12d' "+f+" | cut -f1 | tr , .").read().split('\n')[:-1]]) )
#---------------------------------------------------------------------
######################################             Conversion        #
#---------------------------------------------------------------------
def ConvPu(V) :    return(                           1000*  array(V)           /10                   )
def ConvQa(V) :    return(                            130*  array(V)           /5                    )
def ConvQh(V) :    return(                             20*  array(V)           /5                    )
def ConvHr(V) :    return(                                  array(V)                                 )
def ConvTs(V) :    return(  Tpyro[0]+(Tpyro[-1]-Tpyro[0])*( array(V)-Vpyro[0] )/(Vpyro[1]-Vpyro[0])  )
def ConvQe(V) :    return(  Qfmet[0]+(Qfmet[-1]-Qfmet[0])*( array(V)-Vfmet[0] )/(Vfmet[1]-Vfmet[0])  )
#---------------------------------------------------------------------
def CompM(Qa,Qh2,Qch4) : return( Rho0_a*Qa+Rho0_H2*Qh2+Rho0_CH4*Qch4 )
#---------------------------------------------------------------------
def ConvG1(h) : return( VQH_GT1000[0]+VQH_GT1000[1]*h )
#---------------------------------------------------------------------
######################################             Process           #
#---------------------------------------------------------------------
def FlamePos(Zav,Tav) :
	DTav=diff(Tav) ; ZDav=Zav[:-1]+0.5*diff(Zav)
	return( ZDav[DTav==max(DTav)][0] )
#---------------------------------------------------------------------
def Emit(Zav,Tav) :
	hv=1e4/2
	DZav=diff(Zav)*1e-3
	DTav=diff(Tav)
	DTDZ=DTav/DZav
	epl_u=(DTDZ[ 0]-hv*(Tav[0]-20))/(si*(Tav[ 0]+273.15)**4)
	epl_b=-DTDZ[-1]/(si*(Tav[-1]+273.15)**4)
	util.Section('Emissivity =>   epl_u : {0:.3f}  ,  epl_b : {1:.3f}'.format(epl_u,epl_b),1,0,'y' )
	return(epl_u,epl_b)
#---------------------------------------------------------------------
######################################             Reading           #
#---------------------------------------------------------------------
def LastFile0 (dir    ) : return( os.popen('ls    '+dir+' | tail -n1 ').read().split('\n')[:-1] )
def ListFiles0(dir    ) : return( os.popen('ls    '+dir               ).read().split('\n')[:-1] )
def ListFiles (dir,sup) : return( os.popen('ls    '+dir+sup           ).read().split('\n')[:-1] )
def ListDirs0 (dir    ) : return( os.popen('ls -d '+dir               ).read().split('\n')[:-1] )
def ListDirs  (dir,sup) : return( os.popen('ls -d '+dir+sup           ).read().split('\n')[:-1] )
#####################################################################
def PicSel(Dir,root,Ny,Nz) :
		if   os.path.exists(root+ '--PicSel.txt') : picsel=root+ '--PicSel.txt' #; print('local')
		elif os.path.exists(Dir +   'PicSel.txt') : picsel=Dir +   'PicSel.txt' #; print('direc')
		elif os.path.exists(Dir +'../PicSel.txt') : picsel=Dir +'../PicSel.txt' #; print('gener')
		else                                      : picsel=''                   #; print('nothi')
		if picsel :
			NSel1=array([ int(n) for n in os.popen("sed -n '1p' "+picsel).read().split(',') ])
			NSel2=array([ int(n) for n in os.popen("sed -n '2p' "+picsel).read().split(',') ])
			rot=os.popen('grep Rotate {0} | cut -d: -f2'.format(picsel)).read().split('\n')[0] #; print(rot)
			rvy=os.popen('grep Revery {0} | cut -d: -f2'.format(picsel)).read().split('\n')[0] #; print(rvy)
			rvz=os.popen('grep Reverz {0} | cut -d: -f2'.format(picsel)).read().split('\n')[0] #; print(rvz)
		else      :
			NSel1=array([0,Nz])
			NSel2=array([0,Ny])
			rot=rvy=rvz=''
		return(NSel1,NSel2,rot,rvy,rvz)
#####################################################################
def MarkPic(Dir,IM,root,NSEL,arg) :
	(Nz,Ny,Nc)=IM.shape #; print('IM : ',Nz,Ny,Nc) #; return(IM)
	[[x0,x1],[y0,y1]]=NSEL
	Dcut=int(round(0.9*(x1-x0)/5,0))
	Dpix=50 ; Dpiy=30 
	# Dpix=100 ; Dpiy=100 
	Mark_x0,Mark_x1= max(Dcut+x0-Dpix,0) , min(x1+Dpix,Nz-1)
	# Mark_x0,Mark_x1= 0 , min(x1+Dpix,Nz-1)
	Mark_y0,Mark_y1= max(     y0-Dpiy,0) , min(y1+Dpix,Ny-1)
	# IM3=IM[ Dcut+x0-Dpix : x1+Dpix , y0-Dpiy : y1+Dpix , : ]
	IM3=IM[ Mark_x0:Mark_x1 , Mark_y0:Mark_y1 , : ]
	plt.imsave( root+'--PicCut'+arg+'.JPG' , IM3[:,::-1,:] )
	IM2=copy(IM)
	#print(IM2.shape)
	#print(x0,x1,y0,y1)
	IM2[ x0    , y0:y1 , : ]=255
	IM2[   x1  , y0:y1 , : ]=255
	IM2[ x0:x1 , y0    , : ]=255
	IM2[ x0:x1 ,    y1 , : ]=255
	for x in linspace(x0,x1,6) : IM2[ int( round(x,0) ) , y0:y1 , : ]=255 
	plt.imsave( root+'--PicCol'+arg+'.JPG' , IM2 ) #, cmap=cm.binary_r )
#####################################################################
def CutPic(Dir,IM,root,Arg,Dy) :
		(Nz,Ny,Nc)=IM.shape #; print('IM : ',Nz,Ny,Nc) #; return(IM)
		(NSel1,NSel2,rot,rvy,rvz)=PicSel(Dir,root,Ny,Nz) #; print((rot=='0'),rvy,rvz)
		if (Arg[4] or rot=='1') and ( not rot=='0' ) : IM=flip(transpose(IM,axes=[1,0,2]),axis=0) #; print( 'Rotate  ' )
		if (Arg[6] or rvz=='1') and ( not rvz=='0' ) : IM=IM[::-1, :  ,:]                         #; print( 'Revers z' )
		if (          rvy=='1') and ( not rvy=='0' ) : IM=IM[ :  ,::-1,:]                         #; print( 'Revers y' )
		(Nz,Ny,Nc)=IM.shape #; print('IM : ',Nz,Ny,Nc) #; return(IM)
		if   True                 :  x0,x1=NSel1[0],NSel1[1] ; y0,y1=NSel2[0],NSel2[1] #return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		elif (Nz,Ny)==(3712,5568) :  x0,x1=NSel1[0],NSel1[1] ; y0,y1=NSel2[0],NSel2[1] #return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		elif (Nz,Ny)==(5568,3712) :  x0,x1=NSel1[0],NSel1[1] ; y0,y1=NSel2[0],NSel2[1] #return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		elif (Nz,Ny)==(4128,2752) :  x0,x1=NSel1[0],NSel1[1] ; y0,y1=NSel2[0],NSel2[1] #return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		elif (Nz,Ny)==(1080,1920) :  x0,x1=NSel1[0],NSel1[1] ; y0,y1=NSel2[0],NSel2[1] #return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		elif (Nz,Ny)==(1920,1080) :  x0,x1=NSel2[0],NSel2[1] ; y0,y1=NSel1[0],NSel1[1] #return( IM[ NSel2[0]:NSel2[1] , NSel1[0]:NSel1[1] , : ] )
		elif False : #(Nz,Ny)==(1080,1920) :
			top=IM[  0,:,:] ; (x0,x1)=array(1280*NSel1/3712,dtype=int) #; print(NSel1)   NSel1
			bot=IM[ -6,:,:] ; (y0,y1)=array(1920*NSel2/5568,dtype=int) #; print(NSel2)   NSel2
			IM=insert(IM,[ 1]*100,top,axis=0)
			IM=insert(IM,[-1]*100,bot,axis=0)
			# IM=insert(IM,[Nz+100-1]*100,flip(IM[Nz+100-1,:,:],axis=1),axis=0)
			# return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		elif False : #(Nz,Ny)==(3128,5568) : 
			Dpix=3712-Nz ; Dpix2=int((3712-Nz)*0.5)
			return( IM[ NSel1[0]-Dpix2:NSel1[1]-Dpix2 , NSel2[0]:NSel2[1] , : ] )
		else :
			# [x0,x1]=array([0,Nz],dtype=int)
			# [y0,y1]=array([0,Ny],dtype=int)
			[x0,x1]=array(Nz*NSel1/3712,dtype=int)
			[y0,y1]=array(Ny*NSel2/5568,dtype=int)
			# NSel1=array(Nz*NSel1/3712,dtype=int)
			# NSel2=array(Ny*NSel2/5568,dtype=int)
			# return( IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ] )
		#print(x0,x1,y0,y1)
		arg=''
		if Dy>0 : cy=0.5*(y0+y1) ; y0,y1=int(round(cy-0.5*Dy,0)),int(round(cy+0.5*Dy,0)) ; arg='-Band'
		MarkPic(Dir,IM,root,[[x0,x1],[y0,y1]],arg)
		os.system( 'cp {0} {1}/PICCOL/ '.format(root+'--PicCol.JPG',Dir) )
		return( IM[ x0:x1+1 , y0:y1+1 , : ] )
#####################################################################
def CutProf(SetCut,Zpic,Tpic) :
	[Ncut,Dcut]=SetCut
	Tpcut=zeros(Ncut)
	Zpcut=zeros(Ncut)
	for n in range(Ncut) :
		z=n*Dcut
		Diff=abs(Zpic-z) ; sel=Diff==min(Diff)
		Tpcut[n]=Tpic[ sel ,0][0]
		Zpcut[n]=z
	Xf=FlamePos(Zpcut,Tpcut) ; XT=Zpcut[Tpcut==max(Tpcut)][0]
	return( array(Zpcut) , array(Tpcut) , Xf , XT )
#####################################################################
def AlignPic(IM,ny,Lz) :
	(Nz,Ny,Nc)=IM.shape ; Nb=int(Ny/ny)+1
	IM=mean(IM,axis=2)
	VZ=linspace(0,Lz,Nz)
	VI=array(range(Nz))
	To,Zo,Io=[],[],[]
	fig,ax=plt.subplots()
	for n in range(Nb) : 
		n0,n1=n*ny,min((n+1)*ny,Ny-1)
		Tb=mean( IM[:,n0:n1],axis=1 ) ; P=ones(Nz) #; print(len(Tb),' IM : ',Nz,Ny,Nc)
		(Ts,P)=util.SmoothMeanLoop(Tb,P,50,2) ; Ts=array(Ts) ; MTs=max(Ts) ; SelT=(Ts==MTs)
		Zmax=VZ[SelT]
		ax.plot(VZ,Tb,'k',VZ,Ts,'g',Zmax,MTs,'.r')
		To.append(MTs)
		Zo.append(Zmax)
		Io.append(VI[SelT])
	plt.show()
	sys.exit('by by')
#####################################################################
def SelPic(Pictures) : pass
#####################################################################
def ReadPic(Dir,SetCut,root,Arg,Nsmooth) :
	util.Section('Pic  '+root,0,3,'g')
	# Pictures=[ p for p in ListFiles0(root+'*.JPG') if not '--PicSel.JPG' in p ]
	Pictures=[ p for p in ListFiles0(root+'*.JPG') if not '--Pic' in p ]
	[Ncut,Dcut]=SetCut ; Lz=(Ncut-1)*Dcut
	if Pictures :
		#=====> Reading
		pic=Pictures[-1] ; print(pic)
		IM=plt.imread(pic)                                     #; print('IM :',IM.shape)
		if Arg[9] : 
			IM3=array(1*IM,dtype=uint8)
			IM3[IM3>255]=255
			print( IM3.max() , IM.dtype , IM3.dtype )
			plt.imsave( root+'--PicNorm.JPG' , IM3 )
		####if Arg[4] : IM=flip(transpose(IM,axes=[1,0,2]),axis=0) #; print('IM :',IM.shape)
		####if Arg[6] : IM=IM[::-1,:,:]
		# if Arg[6] : IM=IM[::-1,::-1,:]
		#=====> Cuting
		IMsel=CutPic(Dir,IM,root,Arg,0) ; (nz,ny,nc)=IMsel.shape #IM[ NSel1[0]:NSel1[1] , NSel2[0]:NSel2[1] , : ]
		# AlignPic(IMsel,200,Lz)
		plt.imsave( root+'--PicSel.JPG' , IMsel ) #, cmap=cm.binary_r )
		#=====> Marking
		# MarkPic(Dir,IM,root)
		#=====> Fields
		IMmoy=IMsel # CutPic(Dir,IM,root,Arg,200) ; (nz,ny,nc)=IMmoy.shape
		Zpic =linspace(0,Lz,nz)
		Tpic=flip(mean( IMmoy , axis=1 ),axis=0) # Switch Color and Z => 3Color
		Tav =     mean( Tpic  , axis=1 )         # 1 Color
		#=====> Smoothing
		# Tsmooth=Tpic[:,0]
		Tsmooth=Tav
		for n in range(Nsmooth) : [Tsmooth,Ddum]=util.SmoothMean(Tsmooth,ones(nz),10)
		return( CutProf(SetCut,Zpic,Tpic) , [Zpic,Tpic,Tav,array(Tsmooth),Zpic[Tsmooth==max(Tsmooth)][0],FlamePos(Zpic,Tsmooth)] )
	else        : 
		print('=> No Pic')
		return([[],[],[],[]],[[],[],[],[],[],[]])
#####################################################################
def ListHybTu (M,Nd) : return( sort(list(set(M[:,0]))) , sort(list(set(M[:,1]))) , sort(list(set(M[:,2]))) )
def ListHybTuM(M,Nd) : return( sort(list(set(M[0,:]))) , sort(list(set(M[1,:]))) )
    # hyb,Tu,Phi=[],[],[]
    # for n in range(Nd) :
            # h=M[n,0]
            # p=M[n,1]
            # T=M[n,2]
            # if not h in hyb : hyb.append(h)
            # if not p in Phi : Phi.append(p)
            # if not T in Tu  : Tu. append(T)
    # return(hyb,Phi,Tu)
#####################################################################
def AdiaDataM(M,phi,Tu,Np,Nt) :
	# phi=round(phi,Np)
	Tu =round(Tu ,Nt)
	# sel=(M[0,:]==phi)*(M[1,:]==Tu)
	selT=(M[1,:]==Tu)
	Phi=M[0,selT]
	f_Ta=interp1d( Phi , M[2,selT] , kind='cubic' )
	f_Sl=interp1d( Phi , M[3,selT] , kind='cubic' )
	f_DL=interp1d( Phi , M[4,selT] , kind='cubic' )
	f_Ro=interp1d( Phi , M[5,selT] , kind='cubic' )
	f_IH=interp1d( Phi , M[6,selT] , kind='cubic' )
	f_MH=interp1d( Phi , M[7,selT] , kind='cubic' )
	# if sum(sel) : return( [ s[0] for s in M[:,sel] ] )
	if sum(selT) and min(Phi)<=phi<=max(Phi) : return( [ phi,mean(M[1,selT]),f_Ta(phi),f_Sl(phi),f_DL(phi),f_Ro(phi),f_IH(phi),f_MH(phi) ] )
	else        : 
		print('\033[31m Adia Not Found => phi : {0:.3f}  ,  Tu : {1:.0f}  ,  NselT : {2}  ,  Phi : {3:.3f}  {4:.3f} \033[0m'.format(phi,Tu,sum(selT),min(Phi),max(Phi)) )
		return( zeros(8)              )
#####################################################################
def Mat_Adia( data_adia ) :
	Ph,Tu,Ta,Sl,DL,Ro,IH,MH=[],[],[],[],[],[],[],[]
	# print('Adia Title :')
	# os.system("sed -n '1p' "+data_adia)
	for l in os.popen("sed -e '1d' "+data_adia).read().split('\n')[:-1] : #[:-1] :
		split=[ float(a) for a in l.split(',')[1:] ]
		Ph.append( split[0] )
		Tu.append( split[1] )
		Ta.append( split[2] )
		Sl.append( split[3] )
		DL.append( split[4] )
		Ro.append( split[5] )
		IH.append( split[6] )
		MH.append( split[7] )
	return( array([Ph,Tu,Ta,Sl,DL,Ro,IH,MH]) )
#####################################################################
def Mat_Diff( data_adia ) :
	Ph,Tu,Ta,Sl,DL,Ro,IH,MH,LA,CP=[],[],[],[],[],[],[],[],[],[]
	# print('Adia Title :')
	# os.system("sed -n '1p' "+data_adia)
	for l in os.popen("sed -e '1d' "+data_adia).read().split('\n')[:-1] : #[:-1] :
		split=[ float(a) for a in l.split(',')[1:] ]
		Ph.append( split[0] )
		Tu.append( split[1] )
		Ta.append( split[2] )
		Sl.append( split[3] )
		DL.append( split[4] )
		Ro.append( split[5] )
		IH.append( split[6] )
		MH.append( split[7] )
		LA.append( split[8] )
		CP.append( split[9] )
	return( array([Ph,Tu,Ta,Sl,DL,Ro,IH,MH,LA,CP]) )
#####################################################################
def Adia_Data_M( M_adia,phi,Tu ) :
	SelPh = (M_adia[0,:]==round(phi,2) ) #; print(phi, float(round(phi,12)) )
	SelTu = (M_adia[1,:]==round(Tu ,2) ) #; print(Tu , float(round(Tu ,12)) )
	# print(M_adia[:, SelPh*SelTu ][0])
	# print(phi,sum(SelPh))
	# print(Tu ,sum(SelTu))
	# return( M_adia[:, SelPh*SelTu ] ) # Phi Tu Tad Sl Dl Rho IH MH
	return( [ d[0] for d in M_adia[:,SelPh*SelTu] ] )
#####################################################################
def DSl( M_adia,hyb,phi0,Tu0 ) :
	LTu=list( OrderedDict.fromkeys(M_adia[1]) )
	( Tu,ITu )=util.Nearest( LTu , Tu0 ) ; Sel_Tu = M_adia[1,:]==Tu ; N=len(M_adia[3,:][Sel_Tu])
	( Ph,IPh )=util.Nearest(  M_adia[0,:][Sel_Tu] , phi0 )
	if   0<IPh<N-1 : return( (M_adia[3,:][Sel_Tu][IPh+1]-M_adia[3,:][Sel_Tu][IPh-1])/(M_adia[0,:][Sel_Tu][IPh+1]-M_adia[0,:][Sel_Tu][IPh-1]) )
	elif 0<IPh     : return( (M_adia[3,:][Sel_Tu][IPh  ]-M_adia[3,:][Sel_Tu][IPh-1])/(M_adia[0,:][Sel_Tu][IPh  ]-M_adia[0,:][Sel_Tu][IPh-1]) )
	elif   IPh<N-1 : return( (M_adia[3,:][Sel_Tu][IPh+1]-M_adia[3,:][Sel_Tu][IPh  ])/(M_adia[0,:][Sel_Tu][IPh+1]-M_adia[0,:][Sel_Tu][IPh  ]) )
#####################################################################
def Dad( M_adia,hyb,phi0,Tu0,n ) :
	LTu=list( OrderedDict.fromkeys(M_adia[1]) )
	( Tu,ITu )=util.Nearest( LTu , Tu0 ) ; Sel_Tu = M_adia[1,:]==Tu ; N=len(M_adia[3,:][Sel_Tu])
	( Ph,IPh )=util.Nearest( M_adia[0,:][Sel_Tu] , phi0 )
	if   0<IPh<N-1 : return( (M_adia[n,:][Sel_Tu][IPh+1]-M_adia[n,:][Sel_Tu][IPh-1])/(M_adia[0,:][Sel_Tu][IPh+1]-M_adia[0,:][Sel_Tu][IPh-1]) )
	elif 0<IPh     : return( (M_adia[n,:][Sel_Tu][IPh  ]-M_adia[n,:][Sel_Tu][IPh-1])/(M_adia[0,:][Sel_Tu][IPh  ]-M_adia[0,:][Sel_Tu][IPh-1]) )
	elif   IPh<N-1 : return( (M_adia[n,:][Sel_Tu][IPh+1]-M_adia[n,:][Sel_Tu][IPh  ])/(M_adia[0,:][Sel_Tu][IPh+1]-M_adia[0,:][Sel_Tu][IPh  ]) )
#####################################################################
def found_adia(data_adia,lgrep,Tu0,Phi0) :
	# print(lgrep)
	LSplit=os.popen(lgrep).read().split('\n')[:-1] #; print(LSplit)    # hyb,phi,Tu,Tad,Sl,Thick,Rhou,IHrr,MHrr
	N=len(LSplit)
	if   N==0         : return([  ]) # Tad,Sl,Thick,Rhou,IHrr,MHrr
	if False : #N==0 : #and False :
		Ph,Tu,Ta,Sl,DL,Ro,IH,MH=[],[],[],[],[],[],[],[]
		for l in os.popen("sed -e '1d' "+data_adia).read().split('\n')[:-1] : #[:-1] :
			# print('N=0 : ',l)
			split=[ float(a) for a in l.split(',')[1:] ]
			Ph.append( split[0] )
			Tu.append( split[1] )
			Ta.append( split[2] )
			Sl.append( split[3] )
			DL.append( split[4] )
			Ro.append( split[5] )
			IH.append( split[6] )
			MH.append( split[7] )
		Ph=array(Ph) ; print(Ph[0],Ph[-1],Phi0)
		Tu=array(Tu) ; print(Tu[0],Tu[-1],Tu0 )
		Ta=array(Ta) ; f_Ta=interp2d( Ph,Tu,Ta , kind='cubic' ) # 'linear'
		Sl=array(Sl) ; f_Sl=interp2d( Ph,Tu,Sl , kind='cubic' ) # 'linear'
		DL=array(DL) ; f_DL=interp2d( Ph,Tu,DL , kind='cubic' ) # 'linear'
		Ro=array(Ro) ; f_Ro=interp2d( Ph,Tu,Ro , kind='cubic' ) # 'linear'
		IH=array(IH) ; f_IH=interp2d( Ph,Tu,IH , kind='cubic' ) # 'linear'
		MH=array(MH) ; f_MH=interp2d( Ph,Tu,MH , kind='cubic' ) # 'linear'
		if not Tu[0]<=Tu0<=Tu[-1] : 
			print('\033[31m=> Tu0 not in Tu => Tu : {0:.3f} , {1:.3f}      Tu0 : {2:.3f} \033[0m'.format(Tu[0],Tu[-1],Tu0) )
			Tu1=max( Tu0,Tu[0]  )
			Tu1=min( Tu1,Tu[-1] )
		else : Tu1=Tu0
		return([f_Ta(Phi0,Tu1),f_Sl(Phi0,Tu1),f_DL(Phi0,Tu1),f_Ro(Phi0,Tu1),f_IH(Phi0,Tu1),f_MH(Phi0,Tu1)])
	elif N==1         : 
		# print('N=1 : ',LSplit)
		split=[ float(a) for a in LSplit[0].split(',')[2:] ]
		if abs(split[0]-Tu0)>10 : return([])
		return(split[1:]) # Tad,Sl,Thick,Rhou,IHrr,MHrr
	else              :
		Ph,Tu,Ta,Sl,DL,Ro,IH,MH=[],[],[],[],[],[],[],[]
		for l in LSplit : #[:-1] :
			# print('Else : ',l)
			split=[ float(a) for a in l.split(',')[1:] ]
			Ph.append( split[0] )
			Tu.append( split[1] )
			Ta.append( split[2] )
			Sl.append( split[3] )
			DL.append( split[4] )
			Ro.append( split[5] )
			IH.append( split[6] )
			MH.append( split[7] )
		if   len( set(Tu) )>1 and len( set(Ph) )==1 :
			Ph=array(Ph)
			Tu=array(Tu)
			if min(Tu)-Tu0>10 : return([])
			if min(Tu)==max(Tu) : return([ mean(Ta),mean(Sl),mean(DL),mean(Ro),mean(IH),mean(MH) ])
			Ta=array(Ta) ; f_Ta=interp1d( Tu,Ta , kind='cubic' )
			Sl=array(Sl) ; f_Sl=interp1d( Tu,Sl , kind='cubic' )
			DL=array(DL) ; f_DL=interp1d( Tu,DL , kind='cubic' )
			Ro=array(Ro) ; f_Ro=interp1d( Tu,Ro , kind='cubic' )
			IH=array(IH) ; f_IH=interp1d( Tu,IH , kind='cubic' )
			MH=array(MH) ; f_MH=interp1d( Tu,MH , kind='cubic' )
			if not Tu[0]<=Tu0<=Tu[-1] : 
				print('\033[31m=> Tu0 not in Tu => Tu : {0:.3f} , {1:.3f}      Tu0 : {2:.3f} \033[0m'.format(Tu[0],Tu[-1],Tu0) )
				Tu1=max( Tu0,Tu[0]  )
				Tu1=min( Tu1,Tu[-1] )
			else : Tu1=Tu0
			return([f_Ta(Tu1),f_Sl(Tu1),f_DL(Tu1),f_Ro(Tu1),f_IH(Tu1),f_MH(Tu1)])
		elif len( set(Tu) )==1 and len( set(Ph) )>1 :
			Ph=array(Ph)
			Tu=array(Tu)
			if min(Ph)-Phi0        >0.1 : return([]) # Out of boundaries
			if         Phi0-max(Ph)>0.1 : return([]) # Out of boundaries
			if min(Ph)==max(Ph) : return([ mean(Ta),mean(Sl),mean(DL),mean(Ro),mean(IH),mean(MH) ]) # One Phi
			Ta=array(Ta) ; f_Ta=interp1d( Ph,Ta , kind='cubic' )
			Sl=array(Sl) ; f_Sl=interp1d( Ph,Sl , kind='cubic' )
			DL=array(DL) ; f_DL=interp1d( Ph,DL , kind='cubic' )
			Ro=array(Ro) ; f_Ro=interp1d( Ph,Ro , kind='cubic' )
			IH=array(IH) ; f_IH=interp1d( Ph,IH , kind='cubic' )
			MH=array(MH) ; f_MH=interp1d( Ph,MH , kind='cubic' )
			# if not (Ph[0]<=Phi0<=Ph[-1]) :
			if not (min(Ph)<=Phi0<=max(Ph)) :
				print('\033[31m=> Phi0 not in Phi => Phi : {0:.3f} , {1:.3f}      Phi0 : {2:.3f} \033[0m'.format(min(Ph),max(Ph),Phi0) )
				Phi1=max( Phi0,Ph[0]  )
				Phi1=min( Phi1,Ph[-1] )
			else : Phi1=Phi0
			return([f_Ta(Phi1),f_Sl(Phi1),f_DL(Phi1),f_Ro(Phi1),f_IH(Phi1),f_MH(Phi1)])
		else : return( [] )
#####################################################################
def AdiaData(Dir_adia,hyb,phi,Tu) :
	# data_adia=Dir_adia+'/STe-hyb{0:03.0f}-Tu{1:.0f}.in'.format(hyb*100,Tu)
	data_adia=Dir_adia+'/STe-hyb{0:03.0f}.in'.format(hyb*100,Tu)
	#print( phi,Tu )
	if os.path.exists(data_adia) :
		
		#=====> Phi Tu
		lgrep="grep '{0:.12f},{1:.12f},{2:.12f},' {3}".format( Alpha(hyb) ,round(phi,8),Tu,data_adia) #; print('Phi Tu :',lgrep)
		Data_Adia=found_adia(data_adia,lgrep,Tu,phi)
		if Data_Adia : return( Data_Adia )
		#=====> Tu
		lgrep="grep                  ',{2:.12f},' {3}".format( Alpha(hyb) ,round(phi,8),Tu,data_adia) #; print('Tu     :',lgrep)
		Data_Adia=found_adia(data_adia,lgrep,Tu,phi)
		if Data_Adia : return( Data_Adia )
		#=====> Phi
		lgrep="grep '{0:.12f},{1:.12f},'          {3}".format( Alpha(hyb) ,round(phi,5),Tu,data_adia) #; print('Phi    :',lgrep)
		Data_Adia=found_adia(data_adia,lgrep,Tu,phi)
		if Data_Adia : return( Data_Adia )

		print('No Adia Data Found')
		return(zeros(6))
	else : print('No Adia File : '+data_adia) ; return(zeros(6))
#####################################################################
def AdiaSeparation( M ) :
    Phi0=linspace(0.2,1,81) ; NPhi=len(Phi0)
    NSol=[]
    Tad,Sl0,Dl0=[],[],[]
    for p in Phi0 :
        Selp = M[ (M[:,1]==p) , : ] ; NSol.append(len(Selp))
        Tad.append( Selp[:,3] )
        Sl0.append( Selp[:,4] )
        Dl0.append( Selp[:,5] )
        print(p,NSol[-1])
    
    NSol=array(NSol) ; MaxN=max(NSol) ; print('Max NSol : ',MaxN)
    Phi1=[ [] for n in range(MaxN) ]
    Tad1=[ [] for n in range(MaxN) ]
    Sl01=[ [] for n in range(MaxN) ]
    Dl01=[ [] for n in range(MaxN) ]
    for n in range(NPhi) :
        for q in range(NSol[n]) :
            Phi1[q].append(     Phi0[n]         )
            Tad1[q].append(     Tad [n][-(1+q)] )  #mean([ sort(Tad[n])[::-1][:-2] ]) )
            Sl01[q].append(     Sl0 [n][-(1+q)] )  #mean([ sort(Sl0[n])[::-1][:-2] ]) )
            Dl01[q].append(     Dl0 [n][-(1+q)] )  #mean([ sort(Dl0[n])      [:-2] ]) )
    
    Phi2,Tad2,Sl02,Dl02=[],[],[],[]
    for n in range(NPhi) :
        Phi2.append(     Phi0[n]       )
        Tad2.append(mean(sort( [ x for x in Tad [n][::-1] if x ] )[:-1]))
        Sl02.append(mean(sort( [ x for x in Sl0 [n][::-1] if x ] )[:-1]))
        Dl02.append(mean(sort( [ x for x in Dl0 [n]       if x ] )[:-1]))
        print( '{0:.2f}  ,  {1:.0f}  ,  {2:.3f}  ,  {3:.2f}'.format(Phi2[-1],Tad2[-1],Sl02[-1],Dl02[-1]*1e3) )

    return( [Phi0,Tad,Sl0,Dl0] , [Phi1,Tad1,Sl01,Dl01] , [Phi2,Tad2,Sl02,Dl02] )
#####################################################################
def ThLighter2(time,T) :
    Npac=78
    Ntot=len(time)
    Ndec=int(Ntot/Npac)
    t_out,T_out=[],[]
    for n in range(Ndec) :
        t_out.append(      time[n*Npac           ]  )
        T_out.append( mean(T   [n*Npac:(n+1)*Npac]) )
    if (n+1)*Npac!=Ntot : 
    	# print( '\033[31m Wrong count : Nend {0}  ,  Ntot {1} \033[0m'.format((n+1)*Npac,Ntot) )
    	t_out.append(   time[(n+1)*Npac ]  )
    	T_out.append( mean(T[(n+1)*Npac:]) )
    # t_out.append(time[])
    return(array(t_out),array(T_out))
#####################################################################
def ThLighter(time,T) :
    n=0 ; Ntot=len(time)
    t_out,T_out=[],[]
    while n<Ntot :
        tl=[] ; t=T[n]
        while n<Ntot and T[n]==t :
            tl.append( time[n] )
            n+=1
        t_out.append(mean(tl)) ; print(len(tl))
        T_out.append(     t  )
    return( array(t_out) , array(T_out) )
#####################################################################
def Read(name) :

	tc=float(os.popen('sed -n  "4p" '+name).read())
	fe=float(os.popen('sed -n  "8p" '+name).read())
	dt=float(os.popen('sed -n "10p" '+name).read())
	NL=int  (os.popen("wc -l {0} | cut -d'/' -f1".format(name)).read())-12

	Axes = len( os.popen("grep PosAxe "+name+" | cut -b1-18").read() )>0

	dt0=1/fe
	time=dt0*array( range(NL) ) #; print(dt,time[-1])
	# time=linspace(0,dt,NL) ; dt0=1/fe ; dt1=dt/(NL-1) ; edt=abs(dt0-dt1)/dt0
	# print( '=> {0} {1:.3e} \033[0m'.format( util.ColorVJR( edt , [1e-3,1e-2,1e-1] ) , edt ) )

	MV  =zeros((NL, 4))
	for n in range( 4) : MV[:,n]=[ float(v) for v in os.popen('sed -e "1,12d" {0} | cut -f{1} | tr , .'.format(name,n+1+Axes)).read().split('\n')[:-1] ]
	MT  =zeros((NL,12))
	for n in range(12) : MT[:,n]=[ float(v) for v in os.popen('sed -e "1,12d" {0} | cut -f{1} | tr , .'.format(name,n+5+Axes)).read().split('\n')[:-1] ]

	return(time,MV,MT,[tc,fe,dt,NL])
#####################################################################
def Read2(name) :

	Lines=os.popen('sed -e "12d" '+name+' | tr , . ').read().split('\n')[:-1]

	tc=float( Lines[3] )
	fe=float( Lines[7] )
	dt=float( Lines[9] )
	NL=  len( Lines    )-11

	Axes = len( os.popen("grep PosAxe "+name+" | cut -b1-18").read() )>0

	dt0=1/fe
	time=dt0*array( range(NL) )

	MV  =zeros((NL, 4))
	MT  =zeros((NL,16))
	for n in range(NL) :
		V=Lines[n+11].split('	')
		MV[n,:]=V[Axes*2 :4+Axes*2  ]
		MT[n,:]=V[        4+Axes*2: ]

	return(time,MV,MT,[tc,fe,dt,NL])
#####################################################################
def Read3(name) :

	with open(name,'rb') as fin :
		Lines=fin.readlines()
	fin.closed

	tc=float( Lines[3] ) #; print(tc)
	fe=float( Lines[7] ) #; print(fe)
	dt=float( Lines[9] ) #; print(dt)
	NL=  len( Lines    )-12
	NT=len( Lines[12].decode().split('\t') )-6

	Axes = len( os.popen("grep -a PosAxe "+name+" | cut -b1-18").read() )>0

	dt0=1/fe #; print(dt0)
	time=dt0*array( range(NL) )

	MV  =zeros((NL, 4))
	MT  =zeros((NL,NT))
	for n in range(NL) :
		V=Lines[n+12].decode().replace(',','.').split('\t')
		MV[n,:]=V[Axes*2 :4+Axes*2  ]
		MT[n,:]=V[        4+Axes*2: ]

	return(time,MV,MT,[tc,fe,dt,NL])
#####################################################################
def ReadAuto(name,dx,dt,AV,BND) :
    [N0,N1]=BND ; Nt=sum(BND)

    with open(name,'rb') as fin :
        Lines=fin.readlines()
    fin.closed
    Neuc=1e4

    NT=len( Lines[12].decode().split('\t') )-6
    fe=float( Lines[7] )   
    NL0=len( Lines    )-12 ; NL=NL0-Nt

    time=zeros( NL    )
    VZ  =zeros( NL    )
    MV  =zeros((NL, 4))
    MT  =zeros((NL,NT))
    for n in range(N0,NL0-N1) :
        V=Lines[n+12].decode().replace(',','.').split('\t')
        time[n]=V[0  ]
        VZ  [n]=V[1  ]
        MV[n,:]=V[2:6]
        MT[n,:]=V[6: ]

    VZm =[]
    MVm,MVd =[],[]
    MTm,MTd =[],[]
    if AV :
        for z in VZ :
            Euc=(z*Neuc)%(dx*Neuc)==0
            if Euc and not z in VZm :
                VZm.append(z) ; Sel=(VZ==z)
                Sel_t=time[Sel]>max(time[Sel])-dt
                MVm.append(mean(MV[Sel,:][Sel_t],axis=0)) ; MVd.append(std(MV[Sel,:][Sel_t],axis=0))
                MTm.append(mean(MT[Sel,:][Sel_t],axis=0)) ; MTd.append(std(MT[Sel,:][Sel_t],axis=0))
            #elif not Euc : print( z )
        VZm=array(VZm)
        MVm=array(MVm) ; MVd=array(MVd)
        MTm=array(MTm) ; MTd=array(MTd)

    return(time,[VZ,MV,MT],[VZm,MVm,MTm],[MVd,MTd],[fe,NL])
#####################################################################
def AvAuto(Z,V,dz) :
    Z2,Vm,Vd=[],[],[]
    N=1e4
    for z in Z :
        if ((z*N)/(dz*N))%1==0 and not z in Z2 :
            Z2.append(z) ; Sel=V[Z==z]
            Vm.append(mean( Sel ))
            Vd.append(std ( Sel ))
            if False : #z>30 :
            	fig,ax=plt.subplots() ; ax.set_title('z : {0:.3f}'.format(z))
            	ax.plot( Sel) ; plt.show()
    return(array(Z2),array(Vm),array(Vd))
#####################################################################
def Extract(data,out,time,t0,t1) :
	#out=data[:-4]+'--Extract-t0-{0}-t1-{1}.dat'.format(t0,t1)
	util.Section('Extract : '+out,0,1,'g')
	os.system('head -n12 '+data+' > '+out)
	if t1==0 : t1=time[-1]
	if t0< 0 : t0=time[-1]+t0
	# print(t0,t1)
	Dt0=abs(time-t0) ; n0=list(Dt0).index(min(Dt0))
	Dt1=abs(time-t1) ; n1=list(Dt1).index(min(Dt1))
	# n0,n1 = list(Dt0).index(min(Dt0)),list(Dt1).index(min(Dt1))
	time2=time[n0:n1+1]-t0 #; print(len(time2))
	# os.system("sed -n '{0}p' ".format(n0+13)+data+" | cut -d'\t' -f1")
	# os.system("sed -n '{0}p' ".format(n1+13)+data+" | cut -d'\t' -f1")
	# with open(out,'w') as File :
	for n in range(n1-n0+1) :
		line=os.popen("sed -n '{0}p' ".format(n+13)+data+" | cut -d'\t' -f2-").read()
		os.system( "echo '{0}\t{1}' >> {2}".format(time2[n],line[:-1],out) )
		# File.write('{0},{1}\n'.format(time2[n],line))
	# File.closed
	# os.system("sed -n  '13,{1}p' {2} | cut -f1 >  time.temp ".format(n0+13,13+n1-n0,data) )
	# os.system("sed -n '{0},{1}p' {2} | cut -f1 >  data.temp ".format(n0+13,13+n1   ,data) )
	# os.system("paste time.temp data.temp > temp.temp")
#####################################################################
def CutMT(Cut,MT,Vt,dn) :
	(NT,Nt)=shape(MT) ; Nc=len(Cut)
	OUT=zeros((Nc,NT))
	for n in range(Nc) :
		n0=list(Vt).index(Cut[n])
		OUT[n,:]=mean( MT[ : , max(n0-dn,0) : min(n0+dn+1,Nt) ] , axis=1 )
	return(OUT)
#####################################################################
def TAveraging(MT,I,Tcut,Zbd) :
	[Z0,Z1]=Zbd
	Tav=mean( MT[ I , : ] , axis=1 ) ; Nav=len(Tav) #; Tav[5]=0.5*(Tav[4]+Tav[6])
	#if 11 in I and name=='/work/fmuller/DATA-MANIP/Quartz-ULT-01/Mrange-Phi051-00-Select-01/P3000.dat' : Tav[5]=0.5*(Tav[4]+Tav[6])
	#if 11 in I : Tav[5]=0.5*(Tav[4]+Tav[6])
	Tmi=      MT[ I , : ].min(1)
	Tma=      MT[ I , : ].max(1)
	Tst=std ( MT[ I , : ] , axis=1 )
	Zav=linspace(Z0,Z1,Nav)           ; Sel=Tav<Tcut
	Zsm=linspace(Z0,Z1,Nav*10) if Nav>1 else []
	# [Tsm0,Ddum]=util.SmoothMean(Tav,ones(Nav),1) ; Tav=0.2*array(Tsm0)+0.8*Tav
	# (VP,EP,detP)=util.InterP( Tav , Zav , 5*(Nav==8)+3*(Nav==4) , Nav )
	# Tsm=array(util.EvalP( VP, Zsm ))
	#if I : Xfd=FlamePos(Zav[Sel],Tav[Sel]) ; XT=Zav[Tav==max(Tav[Sel])][0]
	if len(I)>1 : 
		if Nav>3 : f_T=interp1d( Zav,Tav, kind='cubic'  )  ; Tsm=f_T(Zsm)
		else     : f_T=interp1d( Zav,Tav, kind='linear' )  ; Tsm=f_T(Zsm)
		# Tsm,Zsm=Tav,Zav
		Xfd=FlamePos(Zsm,Tsm) ; XT=Zav[Tav==max(Tav[Sel])][0]
		Xsm=Zsm[ Tsm==max(Tsm) ][0]
	else : Xfd=0 ; Xsm=Tmi,Tma,0 ; XT=0 ; Tsm=[]
	return(Zav,Zsm,Tav,Tmi,Tma,Tsm,Tst,Xfd,XT,Xsm)
#####################################################################
def ReadConv(name,FConv,IPlot,Zbd) :
	util.Section('Read '+name,0,3,'b')
	(time,MV,MT,[tc,fe,dt,NL])=Read3(name)
	[ICool,IBiTh,ICent,ISide,IBack,ITu]=IPlot
	Tcut=3000

	#=====> Conv T
    # (nt,nT)=shape(MT) ; nT-=2 ; MT2=[] #; print('MT : ',nt,nT)
	(nt,nT)=shape(MT) ; MT2=[] #; print('MT : ',nt,nT)
	for n in range(nT) :
		(t_l,T_l)=ThLighter2(time,MT[:,n]) #; print('t_l : ',len(t_l),'T_l : ',len(T_l))
		MT2.append( T_l )
	MT2=array(MT2) #; print('MT2 : ',shape(MT2))

	#=====> Conv V
	MV2=array([ FConv[n](MV[:,n]) for n in range(4) ])

	#=====> Tu
	if ITu>=0 : Tu_av=mean( MT2[ ITu , : ] ) ; Tu_sd=std( MT2[ ITu , : ] )
	else      : Tu_av=0 ; Tu_sd=0

	#=====> Bi Th
	Ts_av=mean( MT2[ IBiTh[0] , : ] ) ; Ts_sd=std( MT2[ IBiTh[0] , : ] )
	Tb_av=mean( MT2[ IBiTh[1] , : ] ) ; Tb_sd=std( MT2[ IBiTh[1] , : ] )

	#=====> Pyro
	SelPyro=MV2[0]>Tpyro[0]
	if sum(SelPyro) : Tp_av=mean( MV2[0][SelPyro] ) ; Tp_sd=std( MV2[0] )
	else            : Tp_av,Tp_sd=0,0

	if IBack : Zbd_l,Zbd_r,Zbd_b,ILeft,IRigh=[10,70],[10,70],[0,40],ISide[::2],ISide[1::2]
	else     : Zbd_l,Zbd_r,Zbd_b,ILeft,IRigh=[ 0,60],[10,70],[0,70],ISide[ :4],ISide[4:  ]

	return( [t_l,MT2,time,MV2],TAveraging(MT2,ICent,Tcut,Zbd),TAveraging(MT2,ISide,Tcut,[0,70]),TAveraging(MT2,ISide[::2],Tcut,Zbd_l),TAveraging(MT2,ISide[1::2],Tcut,Zbd_r),TAveraging(MT2,IBack,Tcut,Zbd_b),[Tu_av,Tu_sd,Ts_av,Ts_sd,Tb_av,Tb_sd,Tp_av,Tp_sd] )
#####################################################################
def ReadTb(dirP,rev,title,qa_ref,SAVE,dirp) :
    P=[] ; Tb0,Tb1,Tb2=[],[],[] ; Td0,Td1,Td2=[],[],[] ; Tf=[] ; Qa0,Qh0=[],[] ; Qa,Qh=[],[] ; Qa1,Qh1=[],[] ; Ts=[]
    for f in os.popen('ls '+dirP+'/Z*.dat').read().split('\n')[:-1] :
        print(f)
        P.append( float(f.split('/')[-1].split('.')[0][1:]) )
        #(time,MV,MT,[tc,fe,dt,NL])=Read2(f)
        (time,MV,MT,[tc,fe,dt,NL])=Read3(f)
        
        Tb0.append( mean(MT[:,10]) )  ; Td0.append( std(MT[:,10]) )  #  0 
        Tb1.append( mean(MT[:,16]) )  ; Td1.append( std(MT[:,16]) )  # 10 
        Tb2.append( mean(MT[:,17]) )  ; Td2.append( std(MT[:,17]) )  # 11 

        Tf .append( mean(MT[:, 8]) )

        qa=ConvQa(MV[:,2])/qa_ref
        qh=ConvQh(MV[:,3])/11.206
        Qa0.append( min(qa) ) ; Qa.append( mean(qa) ) ; Qa1.append( max(qa) )
        Qh0.append( min(qh) ) ; Qh.append( mean(qh) ) ; Qh1.append( max(qh) )

        Ts .append( mean(ConvTs(MV[:,0])) )
    
    if rev : P.reverse()
    if title=='X Profile S' : P=linspace(0,100,len(P))
    P=array(P)
    root=os.popen('echo '+title+' | tr " " - ').read().split('\n')[0] ; print(root)
    
    #=========================> Fig Big
    fig,ax=plt.subplots(nrows=2) ; fig.suptitle(title)
    ax2=ax[0].twinx() ; ax3=ax[1].twinx()
    ax[1].set_xlabel('P  [mm]')
    ax[0].set_ylabel('Tb     [C]') ; ax2.set_ylabel('Tf [C]')
    ax[1].set_ylabel('Q/Qset [-]') ; ax3.set_ylabel('Ts [C]')

    ax[0].plot( P,Tb0 , 'r' )
    ax[0].plot( P,Tb1 , 'g' )
    ax[0].plot( P,Tb2 , 'b' )
    ax2  .plot( P,Tf  , 'k' )

    ax[1].plot( P,Qa  , 'g' , P,Qa0,':g' , P,Qa1,':g')
    ax[1].plot( P,Qh  , 'b' , P,Qh0,':b' , P,Qh1,':b')
    ax3  .plot( P,Ts  , 'r' )
    if qa_ref==120 : ax3.set_ylim(700,810)
    if SAVE : util.SaveFig(fig,dirp+'/TProf-'+root+'.pdf')

    #=========================> Fig Small
    if 'X' in title :
        xc0,xc1=30,70
        Sel = (P>=xc0)*(P<=xc1)
        Tb0=array(Tb0)
        Tb1=array(Tb1)
        Tb2=array(Tb2)
        Ts =array(Ts )
        MTb0=mean(Tb0[Sel])
        MTb1=mean(Tb1[Sel])
        MTb2=mean(Tb2[Sel])
        MTs =mean(Ts [Sel])
        fig0,ax0=plt.subplots(nrows=2) ; fig0.suptitle(title+' Zoom')
        ax0[0].set_ylabel('Tb/mean [-]') ; ax0[0].set_xlim(xc0,xc1)
        ax0[1].set_ylabel('Ts/mean [-]') ; ax0[1].set_xlim(xc0,xc1)
        ax0[1].set_xlabel('P  [mm]')
        ax0[0].plot( P[Sel] , Tb0[Sel]/MTb0 , 'o-r' )
        ax0[0].plot( P[Sel] , Tb1[Sel]/MTb1 , 'o-g' )
        ax0[0].plot( P[Sel] , Tb2[Sel]/MTb2 , 'o-b' )
        ax0[1].plot( P[Sel] , Ts [Sel]/MTs  , 'o-r' )
        if SAVE : util.SaveFig(fig0,dirp+'/TProf-'+root+'-Zoom.pdf')

    return( [P,[Tb0,Tb1,Tb2],[Td0,Td1,Td2],Tf,[Qa,Qh],Ts] , [fig,ax,ax2] )
#####################################################################
def ReadTrans(f,title,SAVE,fout) :
    
    (time,MV,MT,[tc,fe,dt,NL])=Read2(f)

    Tb0=MT[:, 0]
    Tb1=MT[:,10]
    Tb2=MT[:,11]
    
    # Tf =MT[:,14]
    Tf =MT[:,2]
    
    Qa =ConvQa(MV[:,2])
    Qh =ConvQh(MV[:,3])
    
    SelTs=MV[:,0]>=Vpyro[0]
    Ts =ConvTs(MV[SelTs,0]) ; timeTs=time[SelTs]
    Hr =       MV[:,1]
    
    fig,ax=plt.subplots(nrows=3) ; fig.suptitle(title)
    ax2=ax[0].twinx() ; ax3=ax[1].twinx()
    ax[0].set_ylabel('Tb     [C]')  ; ax2.set_ylabel('Ts [C]')
    ax[1].set_ylabel('Qa [nL/min]') ; ax3.set_ylabel('Tf [C]')
    ax[2].set_ylabel('PM [V]')
    ax[2].set_xlabel('time   [s]')
    ax[0].set_xlim(0,time[-1]) ; ax[0].set_xticks([])
    ax[1].set_xlim(0,time[-1]) ; ax[1].set_xticks([])
    ax[2].set_xlim(0,time[-1])

    if mean(Tb0)>100 : ax[0].plot( time  ,Tb0 , 'g' )
    ax[0].                   plot( time  ,Tb1 , 'r' )
    ax[0].                   plot( time  ,Tb2 , 'b' )
    ax2  .                   plot( timeTs,Ts  , 'k' )

    ax[1].plot( time  ,Qa  , 'g' )
    # ax[1].plot( time,Qh  , 'b' )
    ax3  .plot( time  ,Tf  , 'r' )

    ax[2].plot( time  ,Hr , 'k' )

    if SAVE : util.SaveFig(fig,fout)
#---------------------------------------------------------------------
######################################             Writing           #
#---------------------------------------------------------------------
def WriteBiTh(time,Ts,Tb,Tps,Tpb,C,Tad,name) :
    N=len(time)
    util.Section(name,0,1,'b')
    print(len(time),len(Ts),len(Tb),len(Tps),len(Tpb),len(C))
    if len(Tps )==1 : Tps =Tps [0]*ones(N)
    if len(Tpb )==1 : Tpb =Tpb [0]*ones(N)
    if len(C   )==1 : C   =C   [0]*ones(N)
    if len(time) :
        ts=Ts/mean(Ts) ; tsa=Ts/Tad
        tb=Tb/mean(Tb) ; tba=Tb/Tad
        Tps=100*array(Tps)
        Tpb=100*array(Tpb)
        C  =100*array(C  )
    else : ts=tb=zeros(N)
    with open(name,'rb','w') as out :
        out.write('time,Ts,Tb,ts,tb,tsa,tba,Tps,Tpb,C \n') # s , K , K , - , - , % , % , %
        for n in range(N) : out.write('{0:.12e},{1:.12e},{2:.12e},{3:.12e},{4:.12e},{5:.12e},{6:.12e},{7:.12e},{8:.12e},{9:.12e} \n'.
                               format(  time[n],   Ts[n],   Tb[n],   ts[n],   tb[n],  tsa[n],  tba[n],  Tps[n],  Tpb[n],    C[n]   ))
    out.closed
#####################################################################
def WritePVariations(time,Ts,Tb,Qa,Qh,md,pt,name) :
    N=len(time)
    with open(name,'rb','w') as out :
        out.write('time,Ts,Tb,Qa,Qh,md,pt \n')
        for n in range(N) : out.write('{0:.12e},{1:.12e},{2:.12e},{3:.12e},{4:.12e},{5:.12e},{6:.12e} \n'.
                               format(  time[n],   Ts[n],   Tb[n],   Qa[n],   Qh[n],   md[n],   pt[n]) )
    out.closed
#####################################################################
def WriteAll(t,M,Iout,Lout,name) : 
	util.Section('Writing : '+name,0,3,'b')
	Nt=len(t) ; NI=len(Iout)
	with open(name,'w') as out :
		out.write('time')
		for l in range(NI) : out.write( ','+Lout[l] )
		out.write('\n')
		for n in range(Nt) :
			out.write('{0:.12f}'.format(t[n]))
			for l in range(NI) : out.write(',{0:.12e}'.format( M[Iout[l],n] ))
			out.write('\n')
	out.closed
#####################################################################
def writeprofile(Zav,Tav,Tst,name) :
	util.Section('Writing : '+name,0,1,'g')
	# print(Zav)
	# print(Tav)
	# print(Tst)
	with open(name,'w') as profdat :
		profdat.write('Z,Tav,Tst\n')
		for n in range(len(Zav)) :
			profdat.write('{0:.12f},{1:.12f},{2:.12f}\n'.format(Zav[n],Tav[n],Tst[n]))
	profdat.closed
#------------------------------------
def WriteProfiles(Data,PicDat,IPlot,Arg,root) :
	util.Section('Data '+root,0,3,'g') ; Tcut=3000
	[ICool,IBiTh,ICent,ISide,IBack,ITu]=IPlot

	([t_l,MT2,time,MV2],Data_c,Data_s,Data_l,Data_r,Data_b,[Tu_av,Tu_sd,Ts_av,Ts_sd,Tb_av,Tb_sd,Tp_av,Tp_sd])=Data
	[Zav_c,Zsm_c,Tav_c,Tmi_c,Tma_c,Tsm_c,Tst_c,Xfd_c,XT_c,Xs_c]=Data_c
	[Zav_s,Zsm_s,Tav_s,Tmi_s,Tma_s,Tsm_s,Tst_s,Xfd_s,XT_s,Xs_s]=Data_s
	[Zav_l,Zsm_l,Tav_l,Tmi_l,Tma_l,Tsm_l,Tst_l,Xfd_l,XT_l,Xs_l]=Data_l
	[Zav_r,Zsm_r,Tav_r,Tmi_r,Tma_r,Tsm_r,Tst_r,Xfd_r,XT_r,Xs_r]=Data_r
	[Zav_b,Zsm_b,Tav_b,Tmi_b,Tma_b,Tsm_b,Tst_b,Xfd_b,XT_b,Xs_b]=Data_b
	writeprofile(Zav_c,Tav_c,Tst_c,root+'--ProfCent.dat') ; writeprofile(Zsm_c,Tsm_c,Tsm_c,root+'--PrsmCent.dat')
	writeprofile(Zav_s,Tav_s,Tst_s,root+'--ProfSide.dat') ; writeprofile(Zsm_s,Tsm_s,Tsm_s,root+'--PrsmSide.dat')
	writeprofile(Zav_l,Tav_l,Tst_l,root+'--ProfLeft.dat') ; writeprofile(Zsm_l,Tsm_l,Tsm_l,root+'--PrsmLeft.dat')
	writeprofile(Zav_r,Tav_r,Tst_r,root+'--ProfRigh.dat') ; writeprofile(Zsm_r,Tsm_r,Tsm_r,root+'--PrsmRigh.dat')
	writeprofile(Zav_b,Tav_b,Tst_b,root+'--ProfBack.dat') ; writeprofile(Zsm_b,Tsm_b,Tsm_b,root+'--PrsmBack.dat')

	([Zpcut,Tpcut,Xfp,XTp],[Zpic,Tpics,Tpic,Tsmooth,Xsm,Xfsm])=PicDat ; Npic=len(Tpic)
	with open(root+'--ProfPict.dat','w') as profpic :
		profpic.write('n,Z,R,G,B,Av,Tsm\n')
		for n in range(Npic) :
			profpic.write('{0},{1:.12f},{2},{3},{4},{5},{6}\n'.format( n,Zpic[n],Tpics[n,0],Tpics[n,1],Tpics[n,2],Tpic[n],Tsmooth[n] ))
	profpic.closed
######################################################################
def WritePos(Data,PicDat,root) :
	([t_l,MT2,time,MV2],Data_c,Data_s,Data_l,Data_r,Data_b,[Tu_av,Tu_sd,Ts_av,Ts_sd,Tb_av,Tb_sd,Tp_av,Tp_sd])=Data
	[Zav_c,Zsm_c,Tav_c,Tmi_c,Tma_c,Tsm_c,Tst_c,Xfd_c,XT_c,Xs_c]=Data_c
	[Zav_s,Zsm_s,Tav_s,Tmi_s,Tma_s,Tsm_s,Tst_s,Xfd_s,XT_s,Xs_s]=Data_s
	[Zav_l,Zsm_l,Tav_l,Tmi_l,Tma_l,Tsm_l,Tst_l,Xfd_l,XT_l,Xs_l]=Data_l
	[Zav_r,Zsm_r,Tav_r,Tmi_r,Tma_r,Tsm_r,Tst_r,Xfd_r,XT_r,Xs_r]=Data_r
	[Zav_b,Zsm_b,Tav_b,Tmi_b,Tma_b,Tsm_b,Tst_b,Xfd_b,XT_b,Xs_b]=Data_b
	([Zpcut,Tpcut,Xfp,XTp],[Zpic,Tpics,Tpic,Tsm,Xsm,Xfsm])=PicDat
	with open(root+'--FlamePos.dat','w') as out :
		out.               write( 'Case,Xf,XT,Xs,Tmax,Tmsm\n' ) #; print( Xfd_c,XT_c,Xs_c,Tav_c[Zav_c==XT_c][0],Tsm_c[Zsm_c==Xs_c][0] )
		if len(Zav_c): out.write( 'Cent,{0:.12f},{1:.12f},{2:.12f},{3:.12f},{4:.12f}\n'.format(Xfd_c,XT_c,Xs_c,Tav_c[Zav_c==XT_c][0],Tsm_c[Zsm_c==Xs_c][0] ) )
		if len(Zav_s): out.write( 'Side,{0:.12f},{1:.12f},{2:.12f},{3:.12f},{4:.12f}\n'.format(Xfd_s,XT_s,Xs_s,Tav_s[Zav_s==XT_s][0],Tsm_s[Zsm_s==Xs_s][0] ) )
		if len(Zav_l): out.write( 'Left,{0:.12f},{1:.12f},{2:.12f},{3:.12f},{4:.12f}\n'.format(Xfd_l,XT_l,Xs_l,Tav_l[Zav_l==XT_l][0],Tsm_l[Zsm_l==Xs_l][0] ) )
		if len(Zav_r): out.write( 'Righ,{0:.12f},{1:.12f},{2:.12f},{3:.12f},{4:.12f}\n'.format(Xfd_r,XT_r,Xs_r,Tav_r[Zav_r==XT_r][0],Tsm_r[Zsm_r==Xs_r][0] ) )
		if len(Zpic) : out.write( 'Pict,{0:.12f},{1:.12f},{2:.12f},{3:.12f},{4:.12f}\n'.format(Xfp  ,XTp ,Xsm ,Tpcut[Zpcut==XTp ][0],Tsm  [Zpic ==Xsm ][0] ) )
	out.closed
#---------------------------------------------------------------------
######################################             Ploting           #
#---------------------------------------------------------------------
def PlotV(time,MV,axV,col) :
    for i in range(2) :
        for j in range(2) :
            axV[i,j].cla()
            axV[i,j].plot(time,MV[:,2*i+j],col)
#####################################################################
def PlotT(time,MT,axT,col) :
    for i in range(4) :
        for j in range(4) :
            axT[i,j].cla()
            axT[i,j].plot(time,MT[:,4*i+j],col)
#####################################################################
def PLOTT(time,MT,AxT,col) : # Plot T versus time
    AxT[0].cla()
    AxT[1].cla()
    AxT[2].cla()
    AxT[0].plot(time,MT[:,10],col+':')
    AxT[0].plot(time,MT[:,11],col+'--')
    for i in range(7) :
        AxT[1].plot(time,MT[:,1+i],color=col,linestyle=(0,(1,1+i*4)))
        AxT[1].plot(time,MT[:,1+i],color=col,linestyle='-',linewidth=0.1)
    for i in range(2) :
        AxT[2].plot(time,MT[:,8+i],color=col,linestyle=(0,(3,1+i*4)))
        AxT[2].plot(time,MT[:,8+i],color=col,linestyle='-',linewidth=0.1)
#####################################################################
def PLOTE(time,MT,AxE,col) : # Plot T parois Vs Deb
    [Deb,Tb,Ts,Tu]=MT
    AxE[0].cla()
    AxE[1].cla()
    AxE[2].cla()
    AxE[0].plot(Deb,Tb[:,0],col)

    for i in range(7) :
        AxE[1].plot(Deb,Ts[i],color=col,linestyle=(0,(1,1+i*4)))
        AxE[1].plot(Deb,Ts[i],color=col,linestyle='-',linewidth=0.1)
    for i in range(2) :
        AxE[2].plot(Deb,Tu[i],color=col,linestyle=(0,(3,1+i*4)))
        AxE[2].plot(Deb,Tu[i],color=col,linestyle='-',linewidth=0.1)
#####################################################################
def PLOTP(Deb,Pic,AxP,col) :
    AxP.cla()
    AxP.set_aspect(1)
    AxP.axis('off')
    title='Hyb : {0} , Pow : {1}'.format(Pic[0][0].split('/')[-3],Pic[0][0].split('/')[-2])
    AxP.set_title( title , fontsize=30 ) #-15*(len(Pic)==1) )
    util.Section(title,0,1,'b')
    Deb1=[]
    Pic1=[]
    Ni,Nj,Nk=3712,5568,3
    Di1,Di2=Ni+5,Ni*0.3
    Dj=1.3*Nj
    Nd=len(Deb)
    for n in range(Nd) :
        if Deb[n] in Deb1 :
            i=Deb1 .index( Deb[n])
            Pic1[i].append(Pic[n])
        else              :
            Deb1.append( Deb[n])
            Pic1.append([Pic[n]])

    Nd1 =len(Deb1)
    Debs=sort(Deb1) ; DSel=[ Deb1.index(d) for d in Debs ]
    Pics=[ Pic1[i] for i in DSel ]

    mi,i0,j0=0,0,0
    for n in range(Nd1) :
        d=Debs[n]
        for pic in Pics[n] :
            if pic :
                cas=PicNameConverter(pic[0]).split('/')[9] ; print('cas : ',cas)
                AxP.text(j0+0.1*Nj,i0+0.3*Di2,cas,fontsize=15)
                for im in pic :
                    print('=> read : ',im)
                    IM=plt.imread(im)
                    AxP.imshow(IM,extent=(j0,j0+Nj,i0-Ni,i0))
                    i0-=Di1
                i0-=Di2
        j0+=Dj ; mi=min(mi,i0) ; i0=0
    jl,jr=0,j0+Nj-Dj         ; dj=abs(jr-jl)
    ib,it=mi-Ni+Di2+Di1,Di2  ; di=abs(it-ib)
    print('i : ',ib,it,di)
    print('j : ',jl,jr,dj)
    AxP.set_xlim((jl,jr))
    AxP.set_ylim((ib,it))
    FigP.set_size_inches(dj/2.5e3 , di/2.5e3 )
#####################################################################
def PlotSave(fig,name,Param,f,SAVE,dpi0) :
    if not ( os.path.exists(name) ) :
        ( time,M,ax,col)=Param
        f(time,M,ax,col)
        if SAVE :
            fig.savefig(name,dpi=dpi0) ; print('=> Saving',name)
    # else : print('=> Allready',name)
#####################################################################
def PlotAll(Data,PicDat,SetCut,Legend,IPlot,Arg,root) :
	util.Section('Plot '+root,0,3,'g')
	([t_l,MT2,time,MV2],Data_c,Data_s,Data_l,Data_r,Data_b,[Tu_av,Tu_sd,Ts_av,Ts_sd,Tb_av,Tb_sd,Tp_av,Tp_sd])=Data
	(                   [Zpcut,Tpcut,      Xfp  ,XTp ],[Zpic ,Tpics,Tpic,Tsm,Xsm,Xfsm])=PicDat ; Npic=len(Zpic)
	[LCool,LBiTh,LCent,LSide,LBack,LVolt]=Legend
	[ICool,IBiTh,ICent,ISide,IBack,ITu]=IPlot
	(Ncut,Dcut)=SetCut ; Zb=(Ncut-1)*Dcut

	[Zav_c,Zsm_c,Tav_c,Tmi_c,Tma_c,Tsm_c,Tst_c,Xfd_c,XT_c,Xs_c]=Data_c
	[Zav_s,Zsm_s,Tav_s,Tmi_s,Tma_s,Tsm_s,Tst_s,Xfd_s,XT_s,Xs_s]=Data_s
	[Zav_l,Zsm_l,Tav_l,Tmi_l,Tma_l,Tsm_l,Tst_l,Xfd_l,XT_l,Xs_l]=Data_l
	[Zav_r,Zsm_r,Tav_r,Tmi_r,Tma_r,Tsm_r,Tst_r,Xfd_r,XT_r,Xs_r]=Data_r
	[Zav_b,Zsm_b,Tav_b,Tmi_b,Tma_b,Tsm_b,Tst_b,Xfd_b,XT_b,Xs_b]=Data_b

	figCool,axCool=plt.subplots(nrows=2) ; bxCool=[    a.twinx() for a in axCool ]
	figBiTh,axBiTh=plt.subplots(       ) ; bxBiTh=axBiTh.twinx()
	figCent,axCent=plt.subplots(nrows=8)
	figSide,axSide=plt.subplots(nrows=8)
	
	figVolt,axVolt=plt.subplots(nrows=4)
	
	figProf,axProf=plt.subplots(nrows=2) ; bxProf=[    a.twinx() for a in axProf ]
	axProf[0].set_ylabel('Tav prof')
	axProf[1].set_ylabel('Tst prof [\%]')

	if   len(ICool)>2 : c=2 ; b=2 ; e=1 #; bxCool=[ a.twinx() for a in axCool ]
	else              : c=1 ; b=0 ; e=0 #; bxCool=[ a.twinx() for a in axCool ]
	if ITu>=0         : 
		axCool[1].plot( t_l , MT2[ ITu , : ] , 'k' ) ; axCool[1].set_ylabel('Tu [dc]')
	d=len(ISide)

	if Arg[10] :
		#=====> Temporal
		for n in range(2*e) : axCool[n].plot( t_l  , MT2[ ICool[c*n  ] , : ] , 'r' ) ; axCool[n].set_ylabel(LCool[c*n  ])
		for n in range(b  ) : bxCool[n].plot( t_l  , MT2[ ICool[2*n+1] , : ] , 'b' ) ; bxCool[n].set_ylabel(LCool[2*n+1])
		for n in range(1  ) : axBiTh   .plot( t_l  , MT2[ IBiTh[  n  ] , : ] , 'r' ) ; axBiTh   .set_ylabel(LBiTh[  n  ])
		for n in range(1  ) : bxBiTh   .plot( t_l  , MT2[ IBiTh[  n+1] , : ] , 'b' ) ; bxBiTh   .set_ylabel(LBiTh[  n+1])
		for n in range(8  ) : axCent[n].plot( t_l  , MT2[ ICent[  n  ] , : ] , 'k' ) ; axCent[n].set_ylabel(LCent[  n  ])
		for n in range(d  ) : axSide[n].plot( t_l  , MT2[ ISide[  n  ] , : ] , 'k' ) ; axSide[n].set_ylabel(LSide[  n  ])
		
		#=====> Voltage
		for n in range(4) : axVolt[n].plot( time , MV2[       n          ] , 'k' ) ; axVolt[n].set_ylabel(LVolt[n])
		for n in range(4) : axVolt[n].plot( time , util.SmoothMean(MV2[n],ones(len(time)),100)[0] , 'r' )
		# [Hr_sm,Ddum]=util.SmoothMean(MV2[1],ones(len(time)),10)
		# axVolt[1].plot( time,Hr_sm,'r' )

	SelT_c=Tav_c<2000 ; Tsd_c=Tst_c/Tav_c
	SelT_l=Tav_l<2000 ; Tsd_l=Tst_l/Tav_l
	SelT_r=Tav_r<2000 ; Tsd_r=Tst_r/Tav_r
	SelT_s=Tav_s<2000 ; Tsd_s=Tst_s/Tav_s
	#=====> Profiles
	if ISide :
		axProf[0].plot( Zav_l , Tav_l     , '-og' , Zav_r , Tav_r     , '-ob' , Zav_s , Tav_s     , '-oy' ) #; axProf[0].set_ylabel('Tav prof')
		axProf[1].plot( Zav_l , Tsd_l*100 , '-og' , Zav_r , Tsd_r*100 , '-ob' , Zav_s , Tsd_s*100 , '-oy' ) #; axProf[1].set_ylabel('Tst prof [%]')
	if ICent :
		#=====> Center
		axProf[0].plot( Zsm_c         , Tsm_c             ,  ':r' ) #; axProf[0].set_ylabel('Tav prof')
		# axProf[0].plot( Zav_c[SelT_c] , Tav_c[SelT_c]     , '-or' ) #; axProf[0].set_ylabel('Tav prof')
		axProf[0].errorbar( Zav_c[SelT_c] , Tav_c[SelT_c] , fmt='-or' , yerr=[Tav_c[SelT_c]-Tmi_c[SelT_c],Tma_c[SelT_c]-Tav_c[SelT_c]] ) #; axProf[0].set_ylabel('Tav prof')
		axProf[1].plot( Zav_c[SelT_c] , Tsd_c[SelT_c]*100 , '-or' ) #; axProf[1].set_ylabel('Tst prof [%]')
		#=====> burned gazs
		axProf[0].plot( Zav_c[-1]+10  , Tp_av             , 'oy'  , Zb    , Ts_av     , 'or'  , Zb    , Tb_av,'ob' )
		#plt.show()
	#=====> Flame pos
	def PlotPos(axProf,Xf,XT,Xs,Tav,Tsd,Sel,c) :
		# print(Tav[Sel],Tsd[Sel])
		axProf[0].plot( 2*[XT],[min(Tav[Sel]),max(Tav[Sel])    ],c , 2*[Xf],[min(Tav[Sel]),max(Tav[Sel])    ],':'+c , 2*[Xs],[min(Tav[Sel]),max(Tav[Sel])    ],'--'+c )
		axProf[1].plot( 2*[XT],[min(Tsd[Sel]),max(Tsd[Sel])*100],c , 2*[Xf],[min(Tsd[Sel]),max(Tsd[Sel])*100],':'+c , 2*[Xs],[min(Tsd[Sel]),max(Tsd[Sel])*100],'--'+c )
	if ICent : PlotPos(axProf,Xfd_c,XT_c,Xs_c,Tav_c,Tsd_c,SelT_c,'r')
	if ISide : PlotPos(axProf,Xfd_l,XT_l,Xs_l,Tav_l,Tsd_l,SelT_l,'g')
	if ISide : PlotPos(axProf,Xfd_r,XT_r,Xs_r,Tav_r,Tsd_r,SelT_r,'b')

	# axProf[0].plot( 2*[Xfd_c],[min(Tav_c[SelT_c]),max(Tav_c[SelT_c])],':r' , 2*[XT_c],[min(Tav_c[SelT_c]),max(Tav_c[SelT_c])],'r' )
	# axProf[1].plot( 2*[Xfd_c],[min(Tsd_c[SelT_c]),max(Tsd_c[SelT_c])],':r' , 2*[XT_c],[min(Tsd_c[SelT_c]),max(Tsd_c[SelT_c])],'r' )
	# print(Xfd,XTd)

	if len(Tpic) :
		for n in range(Ncut) : bxProf[0].plot( 2*[n*Dcut] , [0,255] ,  ':k' , linewidth=0.5 )
		bxProf[0].plot( Zpic  , Tpics[:,0] , 'r'  )
		bxProf[0].plot( Zpic  , Tpics[:,1] , 'g'  )
		bxProf[0].plot( Zpic  , Tpics[:,2] , 'b'  )
		bxProf[0].plot( Zpic  , Tpic       , 'k'  )
		# bxProf[0].plot( Zpic  , Tpics[:,0] , 'k'    )
		bxProf[0].plot( Zpcut , Tpcut      , 'ok'  )
		bxProf[0].plot( Zpic  , Tsm        , ':k'  )
		bxProf[0].plot( 2*[Xfp] , [0,255]  , 'o:k'  )
		bxProf[0].plot( 2*[XTp] , [0,255]  , 'o-k'  )
		bxProf[0].plot( 2*[Xsm] , [0,255]  , 'o--k' )

	#=====> Legend
	for n in range(2-1) : axCool[n].set_xticks([])
	for n in range(2-1) : bxCool[n].set_xticks([])
	for n in range(8-1) : axCent[n].set_xticks([])
	
	for n in range(4-1) : axVolt[n].set_xticks([])
	
	# if Arg[0] : plt.show()
	if Arg[1] and Arg[10] :
		util.SaveFig(figCool,root+'--Cool.pdf')
		util.SaveFig(figBiTh,root+'--BiTh.pdf')
		util.SaveFig(figCent,root+'--Cent.pdf')
		util.SaveFig(figSide,root+'--Side.pdf')
		util.SaveFig(figVolt,root+'--Volt.pdf')
	if Arg[1] :
		util.SaveFig(figProf,root+'--Prof.pdf')

	plt.close('all')
#####################################################################
def PlotScatFlat(Field,r) :
    (x,y,z)=Field
    x-=r
    y-=r

    Ray=(x-0.5*(max(x)-min(x)))**2+(y-0.5*(max(y)-min(y)))**2
    selray=0.75*r
    SelRay=Ray<selray**2
    zst=std( z[SelRay] )
    zav=mean(z[SelRay] )
    zmi=min( z[SelRay] )
    zma=max( z[SelRay] )

    Np=int(1e2)
    # Np=int(1e1)
    # edges = linspace(0, 2*r, Np)
    edges = linspace(-r, r, Np)
    dh=2*r/( Np-1 ) ; s0=(dh*1e-3)**2
    centers = edges[:-1] + diff(edges[:2])[0] / 2.
    XI, YI = meshgrid(centers, centers)
    # RI=(XI-r)**2+(YI-r)**2
    RI=XI**2+YI**2
    Out=RI>r**2    ; Nout=Out.sum()
    Int=abs(Out-1) ; Nint=Int.sum()
    print('N Out : ',Nout)

    print('=> Interpolation')
    rbf = Rbf(x, y, z, epsilon=2)
    ZI = rbf(XI, YI) ; print('ZI : ',shape(ZI),ZI[0,0])
    ZI[Out]=0

    # Cburn=Circle([r,r],r,101)
    # Cav=Circle([r,r],selray,101)
    Cburn=Circle([0,0],   r  ,101)
    Cav  =Circle([0,0],selray,101)
    X_edges, Y_edges = meshgrid(edges, edges) ; print('X_edges : ',shape(X_edges),'   ,   Y_edges : ',shape(Y_edges))
    # lims = dict(cmap='RdBu_r', vmin=0, vmax=max(z))
    lims = dict(cmap='RdBu_r', vmin=0, vmax=ZI.max())

    fig,ax=plt.subplots()
    fig.suptitle('Uav : {0:.3f} m/s   Ust : {1:.2f} , Uma : {2:.2f}  ,  Umi : {3:.2f}  ,  DU : {4:.2f} pc'.format( zav , 100*zst/zav , 100*zma/zav , 100*zmi/zav , 100*(zma-zmi)/zav ) )
    # pcmesh=ax.pcolormesh(X_edges, Y_edges, ZI, edgecolors='', shading='flat', **lims)
    print('=> Plot Interpolation')
    pcmesh=ax.pcolormesh(XI, YI, ZI, edgecolors='', shading='gouraud', **lims)
    # pcmesh=ax.pcolormesh(XI, YI, ZI, edgecolors='', shading='flat', **lims)
    ax.plot(Cburn[0],Cburn[1],'k') ; ax.plot(Cav[0],Cav[1],'w') 
    ax.plot(Cburn[0],Cburn[2],'k') ; ax.plot(Cav[0],Cav[2],'w') 
    print('=> Plot Scatter')
    # ax.scatter(x, y, 25, z, edgecolor='w', lw=0.1, **lims)
    ax.scatter(x, y, 10, z, edgecolor='w', lw=0.1, **lims)
    
    ax.set_title('Velocity Map [m/s]')
    ax.set_xlim(-r,r) #0,2*r)
    ax.set_ylim(-r,r) #0,2*r)
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')

    fig.colorbar(pcmesh)

    #==================================================> Flow rate
    Q0=sum( ZI )*( s0 ) ; Q1=Q0*6e4 ; Q2=Q1/1.1
    print('\n=> Total Flow rate : {0:.3e} m3/s  ,  {1:.2f}  L/min , {2:.2f}  nL/min'.format(Q0,Q1,Q2) )
    #S0=pi*(r*1e-3)**2 ; S1=Nint*s0
    #Umoy_out=Q0/S1
    Umoy_out=mean(z) ; Ustd_out=std(z) ; Up=Ustd_out/Umoy_out
    print('=> Umin : {0:.2f}  m/s ,  Umax : {1:.2f}  m/s'.format(min(z),max(z)))
    print('=> Umoy_out : {0:.2f}  m/s ,  Up : {1:.2f} %'.format(Umoy_out,Up*100))

    # plt.show()
    # fig.savefig(PlotDir+'VelocityMap-ScatFlat.pdf')

    return(XI,YI,ZI,dh,Out,fig)
#---------------------------------------------------------------------
######################################             Optic             #
#---------------------------------------------------------------------
def Collimation01(f,x1,h1) : x0=f**2/x1 ; h0=f*h1/x1 ; return(x0,h0)
def Collimation10(f,x0,h0) : x1=f**2/x0 ; h1=x1*h0/f ; return(x1,h1)
def PinHolRh(f,Lb,Dt,Db) : return( f*(Db-Dt)/(2*Lb) ) # Holle
def PinHolDb(f,Lb,Dt,rh) : return( Dt+2*rh*Lb/f )     # Burner
#---------------------------------------------------------------------
######################################             FlameSpeed        #
#---------------------------------------------------------------------
def LogSlFeat(Sl,Tu) :
	N=len(Sl)
	Sl0=min(Sl) ; Sl2=array(Sl)/Sl0 ; logSl=log(Sl2)
	Tu0=min(Tu) ; Tu2=array(Tu)/Tu0 ; logTu=log(Tu2)
	if any(logSl<0) : print('\033[31m=> Neg logSl \033[0m')
	if any(logTu<0) : print('\033[31m=> Neg logTu \033[0m')
	(P,E0,det)=util.InterP(logSl,logTu,2,N)

	log_pos = (logSl!=0)
	if N-sum(log_pos)!=1 : print('\033[31m=> Zero logSl not 1 \033[0m')
	logSl_out=array(util.EvalP( P,logTu )) ; E_log=abs( logSl_out[log_pos]-logSl[log_pos] )/logSl[log_pos]
	Sl_out=Sl0*exp(P[0])*Tu2**P[1]         ; E_out=abs(      Sl_out         -   Sl          )/   Sl

	# fig,ax=plt.subplots()
	# ax.plot( logTu,logSl    ,'.k' )
	# ax.plot( logTu,logSl_out,'r'  )
	# plt.show()

	return( P,E_log,E_out )
#---------------------------------------------------------------------
def EaZeldovich( Tu,Ta,Sl ) :
	R_g = 8314.4621
	Tu0=min(Tu) ; Sel_Tu=(Tu!=Tu0) ; Tu1=Tu[Sel_Tu]
	Ta0=min(Ta) ; Sel_Ta=(Ta!=Ta0) ; Ta1=Ta[Sel_Ta]
	Sl0=min(Sl) ; Sel_Sl=(Sl!=Sl0) ; Sl1=Sl[Sel_Sl]
	if sum(Sel_Tu*Sel_Ta*Sel_Sl)!=len(Tu)-1 : 
		print( 'Sel Problem => ITu : {0} , ITa : {1} , ISl : {2}'.format( list(Sel_Tu).index(0),list(Sel_Ta).index(0),list(Sel_Sl).index(0) ) )

	C=2*R_g*Ta1*Ta0/( Ta0-Ta1 )
	logTu=log( Tu1     / Tu0     )
	logSl=log( Sl1     / Sl0     )
	logTa=log( Ta1     / Ta0     )
	logDT=log((Ta0-Tu0)/(Ta1-Tu1))
	Ea_Z = C*( logTu-logSl )
	Ea_1 = C*( logTu-logSl + 0.5*3.88*logTa+1.0*logDT )
	Ea_2 = C*( logTu-logSl + 0.5*4.90*logTa+1.5*logDT )

	return( Tu1,Ea_Z,Ea_1,Ea_2 )
#---------------------------------------------------------------------
def SlZeldovich( Ea,Tu,Ta ) :
	R_g = 8314.4621
	Tu0,Tu1=min(Tu),max(Tu) ; Sel_Tu=(Tu!=Tu0) ; Tu1=Tu #[Sel_Tu]
	Ta0,Ta1=min(Ta),max(Ta) ; Sel_Ta=(Ta!=Ta0) ; Ta1=Ta #[Sel_Ta]
	Tu2=Tu1/Tu0
	Sl2=Tu2*exp( -0.5*Ea*(1/Ta1-1/Ta0)/R_g )
	# e= - 0.5*Ea*(1/Ta1-1/Ta0)/R_g #; print( min(e),max(e) )
	# e= -(Ea/(2*R_g))*(1/Ta1-1/Ta0) #; print( min(e),max(e) )
	# e= log(Tu2) - 0.5*Ea*(1/Ta1-1/Ta0)/R_g ; print( min(e),max(e) )
	# Sl2=Tu2*exp( e )
	# Sl2=exp( e )
	return( Tu2,Sl2 )
#---------------------------------------------------------------------
def DSlZeldovich( Ea,Tu,Ta ) :
	R_g = 8314.4621
	Tu0,Tu1=min(Tu),max(Tu) ; Sel_Tu=(Tu!=Tu0) ; Tu1=Tu #[Sel_Tu]
	Ta0,Ta1=min(Ta),max(Ta) ; Sel_Ta=(Ta!=Ta0) ; Ta1=Ta #[Sel_Ta]
	Tu2=Tu1/Tu0
	InvTa=(1/Ta1-1/Ta0)
	C=-0.5/R_g
	DSl2= Tu2*( C*InvTa*exp( C*Ea*InvTa ) )
	return( DSl2 )
#---------------------------------------------------------------------
def DErZeldovich( Ea,Tu,Ta,Sl ) : 
	Sl2=Sl/min(Sl)
	(Tu3,Sl3)=SlZeldovich( Ea,Tu,Ta )
	return( mean( 2*DSlZeldovich( Ea,Tu,Ta )*( Sl3-Sl2 )/Sl2 ) )
#---------------------------------------------------------------------
def ErZeldovich( Ea,Tu,Ta,Sl ) :
# def ErZeldovich( Ea,Param ) :
	# (Tu,Ta,Sl)=Param
	Sl2=Sl/min(Sl)
	(Tu3,Sl3)=SlZeldovich( Ea,Tu,Ta )
	# return( max(abs(Sl3-Sl2)/Sl2) )
	return( mean( ((Sl3-Sl2)/Sl2)**2 ) )
	# return( mean( ((Sl3-Sl2)/Sl2) ) )
#---------------------------------------------------------------------
def HanAlpha( a,Tu,Ta ) :
	R_g = 8314.4621
	Tu0,Tu1=min(Tu),max(Tu) ; Sel_Tu=(Tu!=Tu0) ; Tu1=Tu[Sel_Tu]
	Ta0,Ta1=min(Ta),max(Ta) ; Sel_Ta=(Ta!=Ta0) ; Ta1=Ta[Sel_Ta]
	if sum(Sel_Tu*Sel_Ta)!=len(Tu)-1 : 
		print( 'Sel Problem => ITu : {0} , ITa : {1}'.format( list(Sel_Tu).index(0),list(Sel_Ta).index(0) ) )
	X=(1/Ta0-1/Ta1 )/log(Tu1/Tu0) ; x=1
	Ea=2*R_g*( a-x )/X
	return( Tu1,Ea )
#---------------------------------------------------------------------
def OptiZeldovich( Ea0,Tu,Ta,Sl ) :
	tol=1e-20
	# print(Ea0)
	# Res=minimize( ErZeldovich , Ea0 , args=(Tu,Ta,Sl) , method='BFGS' , tol=1e-20 , options={'disp':True} )
	# Res=minimize( ErZeldovich , Ea0 , args=(Tu,Ta,Sl) , tol=1e-20 , method='Newton-CG' , jac=DErZeldovich , options={'disp':True,'xtol':1e-20} )
	# Res=minimize_scalar( ErZeldovich , bracket=Ea0 , args=(Tu,Ta,Sl) , tol=tol , method='brent' , options={'xtol':tol} )
	# Res=minimize_scalar( ErZeldovich , bracket=Ea0 , args=(Tu,Ta,Sl) , tol=tol , method='golden' , options={'xtol':tol} )
	Res=minimize_scalar( ErZeldovich , args=(Tu,Ta,Sl) )

	# Res=util.Found0_2( ErZeldovich , [Tu,Ta,Sl] , Ea0 , 1e-12 , 10 , True )
	# print(Res)

	# if Res : return( mean(Res) )
	# else   : return(0)
	return( Res.x )
#---------------------------------------------------------------------
def ConfigFact(rb) :

	#=========================> Size
	rp=50e-3
	ro=35e-3
	hb=20e-3
	hs=40e-3
	hh=30e-3
	ht=hs+hh
	ho=ht-hb
	db=2*rb

	#=========================> Surface
	Rh=(rp-ro)/hh
	Sp=  pi*rp**2
	Ss=2*pi*rp*hs
	Sh=  pi*(ro+rp)*hh*sqrt(1+Rh)
	So=  pi*ro**2
	Sb=4*pi*rb**2

	#============================================
	util.Section('Configuration factors',1,5,'r')
	#============================================
	
	#=========================> Sb
	Rbp=rp/hb
	Rbo=ro/ho
	Rbs=rp/hb
	Fbp=0.5*(1-1/sqrt(1+Rbp**2))
	Fbo=0.5*(1-1/sqrt(1+Rbo**2))
	Fbs=       1/sqrt(1+Rbs**2)
	Fbh=1-(Fbp+Fbo+Fbs)
	Sum_b=Fbp+Fbo+Fbs+Fbh
	print( 'Fbp : {0:.3f}  ,  Fbs : {1:.3f}  ,  Fbh : {2:.3f}  ,  Fbo : {3:.3f}  ,  Sum : {4:.3f}'.format(Fbp,Fbs,Fbh,Fbo,Sum_b) )
	
	#=========================> Sp
	Rpf=rp/hs
	Rpt=rp/ht
	Rot=ro/ht
	Xpf=1+(1+Rpf**2)/Rpf**2
	Xpo=1+(1+Rot**2)/Rpt**2
	Fpf=0.5*(Xpf-sqrt(Xpf**2-4))
	Fpb=Fbp*Sb/Sp
	Fpo=0.5*(Xpo-sqrt(Xpo**2-4*(Rot/Rpt)**2))
	Fps=1-Fpb-Fpf
	Fph=1-Fpb-Fpo-Fps
	Sum_p=Fpb+Fpo+Fps+Fph
	print( 'Fps : {0:.3f}  ,  Fph : {1:.3f}  ,  Fpo : {2:.3f}  ,  Fpb : {3:.3f}  ,  Sum : {4:.3f}'.format(Fps,Fph,Fpo,Fpb,Sum_p) )
	
	#=========================> So
	Rfh=rp/hh
	Roh=ro/hh
	Xof=1+(1+Rfh**2)/Roh**2
	Fof=0.5*(Xof-sqrt(Xof**2-4*(Rfh/Roh)**2))
	Fob=Fbo*Sb/So
	Fop=Fpo*Sp/So
	Foh=1-Fof
	Fos=1-Fob-Fop-Foh
	Sum_o=Fob+Fop+Foh+Fos
	print( 'Fop : {0:.3f}  ,  Fos : {1:.3f}  ,  Foh : {2:.3f}  ,  Fob : {3:.3f}  ,  Sum : {4:.3f}'.format(Fop,Fos,Foh,Fob,Sum_o) )
	
	#=========================> Ss
	Fsb=Fbs*Sb/Ss
	Fsp=Fps*Sp/Ss
	Fso=Fos*So/Ss
	Fss=1-Fsb-2*Fsp
	Fsh=1-Fsb-Fsp-Fso-Fss
	Sum_s=Fsb+Fsp+Fsh+Fso+Fss
	print( 'Fsp : {0:.3f}  ,  Fsh : {1:.3f}  ,  Fso : {2:.3f}  ,  Fsb : {3:.3f}  ,  Sum : {4:.3f}  ,  Fss : {5:.3f}'.format(Fsp,Fsh,Fso,Fsb,Sum_s,Fss) )
	
	#=========================> Sh
	Fhb=Fbh*Sb/Sh
	Fhp=Fph*Sp/Sh
	Fhs=Fsh*Ss/Sh
	Fho=Foh*So/Sh
	Fhh=1-Fhb-Fhp-Fhs-Fho
	Sum_h=Fhb+Fhp+Fhs+Fho+Fhh
	print( 'Fhp : {0:.3f}  ,  Fhs : {1:.3f}  ,  Fho : {2:.3f}  ,  Fhb : {3:.3f}  ,  Sum : {4:.3f}  ,  Fhh : {5:.3f}'.format(Fhp,Fhs,Fho,Fhb,Sum_h,Fhh) )
	
	return( array([
	[  0,Fps,Fph,Fpo,Fpb],
	[Fsp,Fss,Fsh,Fso,Fsb],
	[Fhp,Fhs,Fhh,Fho,Fhb],
	[Fop,Fos,Foh,  0,Fob],
	[Fbp,Fbs,Fbh,Fbo, 0 ]
	]) )
#---------------------------------------------------------------------
def BiThermocouple_hat(u1,T,E,rb,MF,gas) :
	[Tu,Tp,Ts,Th,To,Tb]=T          #===> Temperature
	[   ep,es,eh,eo,eb]=E          #===> Emissivity
	[   tp,ts,th,to,tb]=1-array(E) #===> Reflectivity
	[lamb,rocp,visc]=gas
	# tp=1
	# ts=0

	#=========================> Size
	db=2*rb

	#=========================> Gas (phi 0.3 Tu 300)
	alph=lamb/rocp
	Pr=visc/alph
	Re=u1*db/visc
	# Nu= 0.42*Pr**0.2+0.57*Pr**(1/3)*Re**0.5    #===> Cylinder
	Nu= 2 + (0.4*Re**0.5+0.06*Re**(2/3))*Pr**0.4 #===> Sphere
	h0=Nu*lamb/db
	h1=(1-eb)*h0/eb

	[[Fpp,Fps,Fph,Fpo,Fpb],
	[ Fsp,Fss,Fsh,Fso,Fsb],
	[ Fhp,Fhs,Fhh,Fho,Fhb],
	[ Fop,Fos,Foh,Foo,Fob],
	[ Fbp,Fbs,Fbh,Fbo,Fbb]]=MF

	#============================================
	# util.Section('Computation',0,5,'r')
	#============================================
	# print('Adim BiTh => Pr : {0:.3f}  ,  Re : {1:.3f}  ,  Nu : {2:.3f}  ,  h0 : {3:.3f}'.format(Pr,Re,Nu,h0) )

	#=========================> Emittence
	Mb=si*Tb**4
	Mp=si*Tp**4
	Ms=si*Ts**4
	Mh=si*Th**4
	Mo=si*To**4
	
	#=========================> Tensor
	V=array([
	ep*Mp+tp*Fpb*(h1*Tb+Mb),
	es*Ms+ts*Fsb*(h1*Tb+Mb),
	eh*Mh+th*Fhb*(h1*Tb+Mb),
	eo*Mo+to*Fob*(h1*Tb+Mb),
	eb*Mb+        h0*Tb
	])
	
	A=diag([-tp,-ts,-th,-to,eb])
	Id=diag(ones(5)) ; Id[4,4]=h0
	
	M0=array([
	[  0,Fps,Fph,Fpo,-Fpb*h1],
	[Fsp,Fss,Fsh,Fso,-Fsb*h1],
	[Fhp,Fhs,Fhh,Fho,-Fhb*h1],
	[Fop,Fos,Foh,  0,-Fob*h1],
	[Fbp,Fbs,Fbh,Fbo, 0     ]
	])
	
	M=Id+dot(A,M0)
	
	#=========================> Resolution
	det=linalg.det(  M)
	J  =linalg.solve(M,V)
	
	Tg=J[-1]

	return(Tg)
#---------------------------------------------------------------------
def RRE(u1,R,T,eb,gas,Suth) :
	(rb0,rb1)=R
	(Tb0,Tb1)=T
	[lamb,rocp,visc]=gas
	db0=2*rb0
	db1=2*rb1
	alph=lamb/rocp
	Pr=visc/alph
	mu0=Sutherland(Suth,300)
	mub=visc
	Re0=u1*db0/visc
	Re1=u1*db1/visc
	# Nu0=0.42*Pr**0.2+0.57*Pr**(1/3)*Re0**0.5    #===> Cylinder
	# Nu1=0.42*Pr**0.2+0.57*Pr**(1/3)*Re1**0.5    #===> Cylinder
	Nu0= 2 + (mub/mu0)*(0.4*Re0**0.5+0.06*Re0**(2/3))*Pr**0.4 #===> Sphere
	Nu1= 2 + (mub/mu0)*(0.4*Re1**0.5+0.06*Re1**(2/3))*Pr**0.4 #===> Sphere
	# Nu0=0.24+0.56*Re0**0.45 #===> Cylinder small Re
	# Nu1=0.24+0.56*Re1**0.45 #===> Cylinder small Re
	# Nu0=(mub/mu0)*(0.4*Re0**0.5+0.06*Re0**(2/3))*Pr**0.4
	# Nu1=(mub/mu0)*(0.4*Re1**0.5+0.06*Re1**(2/3))*Pr**0.4
	# Nu0=0.35*Re0**0.5+0.052*Re0**(2/3)
	# Nu1=0.35*Re1**0.5+0.052*Re1**(2/3)
	h0=Nu0*lamb/db0
	h1=Nu1*lamb/db1
	print('Adim RRE0 => Pr : {0:.3f}  ,  Re0 : {1:.3f}  ,  Nu0 : {2:.3f}  ,  h0 : {3:.3f}'.format(Pr,Re0,Nu0,h0) )
	print('Adim RRE1 => Pr : {0:.3f}  ,  Re1 : {1:.3f}  ,  Nu1 : {2:.3f}  ,  h1 : {3:.3f}'.format(Pr,Re1,Nu1,h1) )
	# print('eb : ',eb,'   si : ',si)

	Tg=( eb*si*(Tb0**4-Tb1**4)+h0*Tb0-h1*Tb1 )/(h0-h1)
	RRE1=( eb*si*(Tb0+Tb1)*(Tb0**2+Tb1**2) + h0 ) / (h0-h1)
	RRE0=( eb*si*(Tb0+Tb1)*(Tb0**2+Tb1**2) + h1 ) / (h1-h0)

	Tw04= Tb0**4-(Tg-Tb0)*h0/(eb*si)
	Tw14= Tb1**4-(Tg-Tb1)*h1/(eb*si)
	Tw4 = (h0*h1*(Tb1-Tb0)/(eb*si)+h0*Tb1**4-h1*Tb0**4)/(h0-h1)
	print('Tw04 : {0:.3f}  ,  Tw14 : {1:.3f}  ,  Tw4 : {2:.3f}'.format(Tw04,Tw14,Tw4) )
	h0_max=eb*si*Tb0**4/(Tg-Tb0)
	h1_max=eb*si*Tb1**4/(Tg-Tb1)
	print( 'h0_max : {0:.3f}  ,  h1_max : {1:.3f}'.format(h0_max,h1_max) )

	Tw0=( abs(Tw04) )**0.25
	Tw1=( abs(Tw14) )**0.25

	return(RRE0,RRE1,Tg,Tw0,Tw1)
#---------------------------------------------------------------------
def RRE_It(u0,R,T,eb,gas) :
	[Pr,Suth,Fresh]=gas
	[T0,ro0,cp]=Fresh

	Tg0=0
	Tg1=T[0]
	while abs(Tg1-Tg0)>1e-1 :
		mu=Sutherland(Suth,Tg1)
		ro=ro0*T0/Tg1
		ub=u0*Tg1/T0
		al=mu/Pr
		la=ro*cp*al
		Tg0=Tg1
		gas=[la,ro*cp,mu]

		(RRE0,RRE1,Tg1,Tw0,Tw1)=RRE(ub,R,T,eb,gas,Suth)
		print('Iteration => Tg0 : {0:.3f}  ,  Tg1 : {1:.3f}'.format(Tg0,Tg1) )

	print('\033[32m Convergence => Tg : {0:.3f} \033[0m'.format(Tg1) )
	return(RRE0,RRE1,Tg1,Tw0,Tw1,gas,ub)
#---------------------------------------------------------------------
def La_Ray( poro,dp,T ) : ext=2.656*sqrt(1-poro)/dp ; return( 16*si*T**3/(3*ext) )
#---------------------------------------------------------------------
def Sutherland(Set,T) :
	[mu0,T0,C]=Set
	cst=mu0*(T0+C)/T0**1.5
	return( cst*T**1.5/(T+C) )
#---------------------------------------------------------------------
def Younis_0(dp) : #,L) :
	C=-399.995758758165*dp + 0.687196666384 #fY_C(dp) #0.819*( 1- 7.33*(dp/L) )
	m= 443.739927050641*dp + 0.361220417338 #fY_m(dp) #0.36 *( 1+15.50*(dp/L) )
	return(C,m)	
#---------------------------------------------------------------------
def Younis_1(ro,U,mu,la,dp) :
	(C,m)=Younis_0(dp) ; print('C : {0:.3f}  ,  m : {1:.3f}'.format(C,m))
	Re=ro*U*dp/mu
	Nu_y=max(C*Re**m,0)
	hy=la*Nu_y/dp**2
	return(hy)
#---------------------------------------------------------------------
def Younis(U,dp,T) :
	(f_ro,f_mu,f_nu,f_cp,f_la,f_aa,f_pr)=AirProp()
	ro=f_ro(T)
	la=f_la(T)
	mu=f_mu(T)
	return(Younis_1(ro,U,mu,la,dp))
#---------------------------------------------------------------------
def ph_Jung(Phi) : return(   Phi/(2-Phi) )
def Ph_Jung(phi) : return( 2*phi/(phi+1) )
#---------------------------------------------------------------------
def Sl_Jung(Phi) : return( -4.7-245.3*Phi+926.8*Phi**2-440.3*Phi**3 )
#---------------------------------------------------------------------
def Dq_Jung(Phi) :                    return( exp( 6.49-22.83*Phi+29.79*Phi**2-18.47*Phi**3+4.59*Phi**4 ) )
def dq_Jung(phi) : Phi=Ph_Jung(phi) ; return( exp( 6.49-22.83*Phi+29.79*Phi**2-18.47*Phi**3+4.59*Phi**4 ) )
#---------------------------------------------------------------------
def Pe_Jung(phi,al) :
	Phi=Ph_Jung(phi)
	return( Dq_Jung(Phi)*Sl_Jung(Phi)/al )
#---------------------------------------------------------------------
def AirProp() :
	f_ro=interp1d( AirData[:,0] , AirData[:,1] , kind='cubic' )
	f_mu=interp1d( AirData[:,0] , AirData[:,2] , kind='cubic' )
	f_nu=interp1d( AirData[:,0] , AirData[:,3] , kind='cubic' )
	f_cp=interp1d( AirData[:,0] , AirData[:,4] , kind='cubic' )
	f_la=interp1d( AirData[:,0] , AirData[:,5] , kind='cubic' )
	f_aa=interp1d( AirData[:,0] , AirData[:,6] , kind='cubic' )
	f_pr=interp1d( AirData[:,0] , AirData[:,7] , kind='cubic' )
	return(f_ro,f_mu,f_nu,f_cp,f_la,f_aa,f_pr)
#---------------------------------------------------------------------
def GasProp(File_diff,Tu) :
	(M_diff,tit_diff)=util.ReadFile(File_diff ,',',1) ; M_diff=array(M_diff)
	Tu0=M_diff[ : , 2] ; Sel=(Tu0==Tu)
	Phi=M_diff[Sel, 1]
	Rho=M_diff[Sel, 6]
	Lau=M_diff[Sel, 9]
	Cpu=M_diff[Sel,10]
	Alp=Lau/(Rho*Cpu)
	f_ro=interp1d( Phi,Rho , kind='cubic' )
	f_la=interp1d( Phi,Lau , kind='cubic' )
	f_cp=interp1d( Phi,Cpu , kind='cubic' )
	f_al=interp1d( Phi,Alp , kind='cubic' )
	return(f_ro,f_la,f_cp,f_al)
#---------------------------------------------------------------------
def ReadTemp(name) : return([ float(x) for x in os.popen("sed -n '2p' "+name).read().split(',') ])
#---------------------------------------------------------------------
def Um_Pois(Q,D,T) : return( 2*Vit(Q,D,T) )
# def Um_Pois(Q,D) : return( 2*Q/( pi*(0.5*D)**2 ) )
# def Um_Pois(Q,D) : return( 2*Q/( pi*(0.5*D)**4 ) )
#---------------------------------------------------------------------
def Ur_Pois(r,D,Q) :
	Um=Um_Pois(Q,D,T)
	# return( Um*( (0.5*D)**2 - r**2 ) )
	return( Um*( 1-(2*r/D)**2 ) )
#---------------------------------------------------------------------
def rU_Pois(u,D,Q,T) :
	Um=Um_Pois(Q,D,T)
	if u>Um : 
		print('u : {0:.3f}  ,  Um : {1:.3e}'.format(u,Um))
		return(0)
	# return( sqrt( (0.5*D)**2-u/Um ) )
	return( sqrt(1-u/Um)*(0.5*D) )
#---------------------------------------------------------------------
def G_Pois(Q,D,Tg,Tu) : return( (Tg/Tu)*4*(Q/6e4)/(pi*(0.5*D)**3) )
#---------------------------------------------------------------------
# def DrG_Pois( Qa,QH2,QCH4,D ) : return( 4*( DaQtot(Qa,QH2,QCH4)/(Qa+QH2+QCH4) + 6*Ea_r/D )/pi )
def DrG_Pois( Qa,QH2,QCH4,D ) : return( 4*( DaQtot(Qa,QH2,QCH4)/(Qa+QH2+QCH4) + 3*1e-2 )/pi )
#---------------------------------------------------------------------
def Ze(Ea,Tu,Tad,R) : return( Ea*(Tad-Tu)/(R*Tad**2) )
#---------------------------------------------------------------------
def Exp(Tu,Tb) : return( Tb/Tu )
#---------------------------------------------------------------------
def Cb_Le(Ze,phi) : return( 1 + Ze*(1/phi-1) )
#---------------------------------------------------------------------
def Le_eff(LeE,LeD,phi,Ze) : 
	Cb=Cb_Le(Ze,phi)
	return( 1+(LeE-1+Cb*(LeD-1))/(1+Cb) )
#---------------------------------------------------------------------
def Bet_M(Exp,Ze,Le) : return( Exp+0.5*Ze*(Le-1) )
#---------------------------------------------------------------------
def LM(LeE,LeD,phi,Ea,Tu,Tad,Tb,Lf,R) :
	ze =Ze(Ea,Tu,Tad,R)
	exp=Exp(Tu,Tb)
	Le =Le_eff(LeE,LeD,phi,ze)
	Bet=Bet_M(exp,ze,Le)
	return( Lf*(Bet-(exp-1)) )
#---------------------------------------------------------------------
def Hoferichter(dq,dp,LM) : 
	# print('dq : {0:.3f} mm  ,  dp : {1:.3f} mm  ,  dp/dq : {2:.3f}'.format(dq*1e3,dp*1e3,dp/dq))
	# print('LM : {0:.3f} mm'.format(LM*1e3) )
	y2=1-2*dq/dp
	u=0.5/( 1 + (4*LM/dp)*y2 - y2**2 )
	print('\033[31m U/Sl Hofe Output : {0:.3f} \033[0m'.format(u) )
	return( u )
#---------------------------------------------------------------------
def Hoferichter2(dp,LM,U,Sl) :
	print('\033[32m U/Sl Hofe Input : {0:.3f} \033[0m'.format(U/Sl) )
	b=-4*LM/dp
	c=0.5*Sl/U-1
	D=b**2-4*c #; print('Delta : ',D)
	x1=(-b-sqrt(D))/2
	x2=(-b+sqrt(D))/2
	# print('x1  :',x1)
	# print('x2  :',x2)
	# dq1=dp*(x1+1)/2
	# dq2=dp*(x2+1)/2
	dq1=0.5*dp*(1-x1)
	dq2=0.5*dp*(1-x2)
	# print('dq1 : ',dq1)
	# print('dq2 : ',dq2)
	return( dq1,dq2 )
#---------------------------------------------------------------------
def Clavin_MaD(Tu,Tb,Ea,Lef) : 
	bet=(Tb-Tu)*Ea/(R*Tb**2)
	x0=linspace(1e-8,(Tb-Tu)/Tu,int(1e3))
	y0=log(1+x0)/x0 
	return( Tb*log(Tb/Tu)/(Tb-Tu) + 0.5*bet*(Lef-1)*Tu*trapz(x0,y0)/(Tb-Tu) )
def Clavin_MaC(Tu,Tb,Ea,Lef) : 
	bet=(Tb-Tu)*Ea/(R*Tb**2)
	x0=linspace(1e-8,(Tb-Tu)/Tu,int(1e3))
	y0=log(1+x0)/x0 
	return(                         0.5*bet*(Lef-1)*Tu*trapz(x0,y0)/(Tb-Tu) )
#---------------------------------------------------------------------
def SlH2_Hof(phi)  : return( 2.1600*phi**4-9.1537*phi**3+12.4930*phi**2-3.7952*phi+0.3972 )
def SlCH4_Hof(phi) : return( 2.1494*phi**4-9.3553*phi**3+13.3410*phi**2-7.0035*phi+1.2279 )
#---------------------------------------------------------------------
def LewisVonElbe( Pe,Sl,dh,LM,dp ) :
	dq=Pe*dh ; dq2=1-2*dq/dp
	gc0=Sl/dq
	gc1=Sl/(dq+LM)
	gc2=Sl/( 0.25*dp*(1-dq2**2)+LM*dq2 )
	return( gc0,gc1,gc2 )
#---------------------------------------------------------------------
def Omega(Y,T,B,Ta) : return( B*Y*exp(-Ta/T) )
# def Omega(Y,T,Q,Ta) :
    # C=1-Y/Y[0]
    # om=Y*exp(-Ta/T) ; Iom=trapz(C,om) ; B=Q/Iom ; print(B)
    # return( B*om )
#---------------------------------------------------------------------
def Check_Sl( Phi,Sl ) :
    IM_Sl=plt.imread('/work/fmuller/Python/Manip/Jung/FlameSpeed.png') ; [Ni_Sl,Nj_Sl,Nk_Sl]=shape(IM_Sl)
    x0,x1=0,2
    y0_Sl,y1_Sl=-9,400
    fig,ax=plt.subplots()
    im=ax.matshow( IM_Sl , extent=(x0,x1,y0_Sl,y1_Sl) , cmap=cm.binary_r )
    ax.set_xticks([0,0.5,1,1.5,2]) ; ax[0].set_aspect( (Ni_Sl/Nj_Sl)*(x1-x0)/(y1_Sl-y0_Sl) )
    plt.show()
