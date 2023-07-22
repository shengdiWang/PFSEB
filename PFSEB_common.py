#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 14:16:41 2023

@author: shengdiwang
"""

import numpy as np
import pandas as pd
from os import path
import matplotlib.pyplot as plt
import math
import cmath
import sympy
from sympy import *
from sympy import Function, Symbol
from scipy import integrate




def snowThermalCap(roSnow):
    # snow thermal capacity from HTESSEL

    roIce = 920
    cIce  = 2.05
    snCap = roSnow * (cIce/roIce)     
    
    return snCap


def snowThermalCon(roSnow):
    # snow thermal conductivity from HTESSEL

    roIce = 920
    kIce  = 2.29
    snTC = kIce * (roSnow / roIce)**1.88 # snow conductivity
    
    return snTC


def snowAlbedo(roSnow):
    # snow albedo from Eq. 31

    snalbedo = 1 - 0.247 * math.sqrt(0.16 + 110 * (roSnow/1000)**4)

    return snalbedo


def atmosphericVaporPressure(t):
    # Atmospheric vapor pressure as Eq. 3
    # avp = 10 ** (11.40 - 2353.0 / t)
    avp = 611 * sympy.exp(17.62 * (t - 273.15) / (243.12 + (t - 273.15)))
    
    return avp


def rho_a(Tz, p):

    R_a = 287.058 # specific gas constant of air [ J/(kg K) ]
    res = p/(R_a*Tz) # air density [kg m^(-3)]

    return res


# Solve the equation by Newton-Raphson method in order to estimate TS
def newton(function, derivative, x0, tolerance, max_iterations):
    x1 = 0
    if abs(x0-x1)<= tolerance and abs((x0-x1)/x0)<= tolerance:
        return x0
    
    k = 1
    while k <= max_iterations:
        x1 = x0 - (function(x0)/derivative(x0))
        if abs(x0-x1)<= tolerance and abs((x0-x1)/x0)<= tolerance:
            return x1
        x0 = x1
        k = k + 1
        # Stops the method if it hits the number of maximum iterations
        if k > max_iterations:
            print("ERROR: Exceeded max number of iterations")
            return x1 # Returns the value


def longWaveOut(emiV, Ts0):
    # Emitted longwave radiation as Eq. 4

    SIGMA = 5.670374E-8 # Stefan-Boltzmann constant [W m-2 K-4]
    LWout = -emiV * SIGMA * Ts0 ** 4
    
    return LWout


# Sensible heat flux    
def sensible_heat(roAir, CAir, ra, TA, TS):
    
    QH = - roAir * CAir / ra * (TA - TS)
   
    return QH


def latent_heat(roAir, latSub, ra, EA, ES0, PA):
    
    QE = - roAir * latSub / (ra + 50) * (0.622 * (EA - ES0) / PA)
    
    return QE


# stable
def psi_H(zeta1, zeta2):

    if zeta1 < 0:
        
        # def funH(x):
        #     return 0.95 * (1 - 11.6 * x) ** (-0.5)
        # res = integrate.quad(funH, zeta2, zeta1)[0]
        
        res = 1.9 * cmath.atanh((1 - 11.6 * zeta1) ** 0.5) + cmath.log(zeta1) \
              - (1.9 * cmath.atanh((1 - 11.6 * zeta2) ** 0.5) + cmath.log(zeta2))
    
    elif zeta1 > 0:
        
        def funH(x):
            return 0.95 + 7.8 * x
        res = integrate.quad(funH, zeta2, zeta1)[0]
        # res = 0
    else:
        res = 0
    
    return res
    # return 0


# stable
def psi_M(zeta1, zeta2):

    if zeta1 < 0:
        
        # def funM(x):
        #     return (1 - 19 * x) ** (-0.25)
        # res = integrate.quad(funM, zeta2, zeta1)[0]
        
        res = -2 * atan((1 - 19 * zeta1) ** (1/4)) + 2 * math.log(1 + (1 - 19 * zeta1) ** (1/4)) + math.log(1 + (1 - 19 * zeta1) ** 0.5) \
              - (-2 * atan((1 - 19 * zeta2) ** (1/4)) + 2 * math.log(1 + (1 - 19 * zeta2) ** (1/4)) + math.log(1 + (1 - 19 * zeta2) ** 0.5))
    
    elif zeta1 > 0:
        
        def funM(x):
            return 1 + 6 * x
        res = integrate.quad(funM, zeta2, zeta1)[0]
        # res = 0

    else:
        res = 0
    
    return res
    # return 0


# Turbulent coefficients       
def atmospheric_correction(L, z, d, zom): # Lstar, wsHeight, 0, RZ
# stable               
    if L>0 :
         sai_m = sai_h = -5 * ((z-d)/L)
         return sai_m , sai_h
# unstable          
    if L<0 :
         x = y = float(pow((1 - (16 * (z-d)/L)),0.25))
         x0 = float(pow((1 - (16 * ((z-d)/L) * (zom/z))),0.25))
         y0 = float(pow((1 - (16 * ((z-d)/L) * ((0.1 * zom)/z))),0.25))
         sai_m = np.log ((1 + pow(x,2))/(1 + pow(x0,2))) + 2 * np.log ((1 + x)/(1 + x0)) - 2 * np.arctan(x) + 2 * np.arctan(x0)
         sai_h = 2 * np.log ((1 + y)/(1 + y0))
         return sai_m , sai_h        
# neutral         
    if L==0:
         sai_m = sai_h = 0
         return sai_m, sai_h
     
    return 0, 0


# Aerodynamic resistance for soil                                                                            
def ra_air(wsHeight, RZ, Lstar, wsi):
    
    if Lstar == 0:
        return 1000000000000000
    else:
        zeta1 = wsHeight / Lstar
        zeta2 = RZ / Lstar
        Ustar = wsi * VONK / (math.log(wsHeight / RZ) - psi_M(zeta1, zeta2))
        ra = (math.log(wsHeight / RZ) - psi_H(zeta1, zeta2)) / (VONK * Ustar)
        return ra
    
    # if Lstar == 0:
    #     return 1000000000000000
    # else:
    #     zeta1 = wsHeight / Lstar
    #     zeta2 = RZ / Lstar
    #     Ustar = wsi * VONK / (math.log(wsHeight / RZ) - atmospheric_correction(Lstar, wsHeight, 0, RZ)[0])
    #     ra = (math.log(wsHeight / RZ) - atmospheric_correction(Lstar, wsHeight, 0, RZ)[1]) / (VONK * Ustar)
    #     return ra


def Monin_obokhov(roAir, Ustar, TA, QH, QE):
    
    Lstar = roAir * CAir * Ustar ** 3 * TA / (VONK * GRAVIT * (QH + 0.61*CAir*TA*QE/latSub))
    
    if Lstar > 1e6:
        Lstar = 1e6
    elif Lstar < -1e6:
        Lstar = -1e6
    elif 0 < Lstar < 1e-6:
        Lstar = 1e-6
    elif -1e-6 < Lstar < 0:
        Lstar = -1e-6
    else:
        Lstar = Lstar
    # elif 0 < Lstar < 10:
    #     Lstar = 10
    # elif -10 < Lstar < 0:
    #     Lstar = -10
    # else:
    #     Lstar = Lstar

    return Lstar


def ustar(wsi, wsHeight, RZ, Lstar): 
    
    zeta1 = wsHeight / Lstar
    zeta2 = RZ / Lstar
    
    return (wsi * VONK)/(math.log(wsHeight / RZ) - psi_M(zeta1, zeta2))

    # return (wsi * VONK)/(math.log(wsHeight / RZ) - atmospheric_correction(Lstar, wsHeight, 0, RZ)[0])





# === PROFILE ===

# 5-layer snow, 16 layers for 0-1 m, 190 layers for 1-20 m [0.1 m interval]


soil1m = [0, 3, 7, 11, 15, 22, 29, 35, 44, 50, 58, 66, 75, 84, 93, 100]
soil1m = [i / 100 for i in soil1m]
soilTck = 0.1 # thickness for deep soil (> 1 m)


snowN = 5   # maximum snow layer
SLTck1 = 1 # soil layer 1 thickness
SLTck2 = 9 # soil layer 2 thickness
soilDepth = 10  # total soil depth
NODE = int(snowN + len(soil1m) + (soilDepth - 1) / soilTck) # layer number
intiST = -3 # initial soil temperature [C]


# K soil
RKKK = np.zeros([7, 1])
RKKK[6, 0] = 2.92 * 24.0 * 3600.0


# C soil
RCCC = np.zeros([7, 1])


# snow & soil node depth
X = np.zeros([NODE, 1])
for i in range(snowN, NODE):
    if i < snowN + len(soil1m): # 5-layer snow & 0-1 m soil
        X[i, 0] = soil1m[i - snowN] 
    else:                       # soil profile > 1 m
        X[i, 0] = 1.0 + (i - snowN - len(soil1m) + 1) * soilTck;
        

# snow & soil node thickness
DX = np.zeros([NODE - 1, 1])
for i in range(snowN, NODE):
    DX[i - 1, 0] = X[i, 0] - X[i - 1, 0]
    

# === CONSTANT ===

Tzero = 0 # frozen temperature
QQ    = 0 # the lower boundary condition

# AIR
roAir = 1.225  # air density [kg m-3]
CAir  = 1003.5  # air thermal capacity [J m-3 C-1]

# WATER
roWater = 1000.0 # water density

# ICE
roIce = 920   # ice density
cIce  = 2.05  # ice volumetric heat capacity [MJ m-3 K-1]
kIce  = 2.29  # thermal conductivity of ice  [W m-1 K-1]

GRAVIT   = 9.807   # gravitational acceleration
VONK     = 0.4     # Von Karman's constant
TF       = 273.15  # unit C to K
latSub   = 2.5E6 # latent heat of sublimation [J kg-1]
wsHeight = 2.0     # REFERENCE HEIGHT for wind speed measure
albedoG  = 0.2     # snow-free albedo
SIGMA    = 5.670374E-8





class SEB(object):


    def __init__(self, metaf, site = 1):
        
        self.meta = pd.read_csv(metaf)
        self.meta = self.meta[self.meta['site'] == site]
        
        # simulation period
        self.begYr      = int(self.meta.beg)
        self.endYr      = int(self.meta.end)
        
        self.begDate1   = pd.to_datetime(str(self.begYr) + '0101')     # start date for run
        self.endDate1   = pd.to_datetime(str(self.begYr + 5) + '1231') # end date for run
        self.dateRange1 = pd.date_range(self.begDate1, self.endDate1)
        
        self.begDate2   = pd.to_datetime(str(self.begYr) + '0101')     # start date for run
        self.endDate2   = pd.to_datetime(str(self.endYr) + '1231')     # end date for run
        self.dateRange2 = pd.date_range(self.begDate2, self.endDate2)
        
        
        
        # ---- import forcing ----
        
        self.focf        = str(site) + '.csv'
        self.focf        = path.join(dir_foc, self.focf);
        self.foc         = pd.read_csv(self.focf)

        self.DATES       = pd.to_datetime(self.foc.date)
        self.AIRT        = self.foc.satFinal         # air temperature
        self.SND         = self.foc.snd              # snow depth
        self.dewP        = self.foc.tempD            # dew-point temperature
        self.SRin        = self.foc.sw               # solar radiation
        self.TRin        = self.foc.lwin             # theraml radiation
        self.WS          = self.foc.ws               # wind speed
        self.PRE         = self.foc.Pa               # pressure
        self.snowRo      = self.foc.sdn              # snow density
        self.emi         = self.foc.emi              # emissivity
        self.foc['z']    = 0.015 # 0.001 2E-5
        self.foc.loc[(self.foc.snd >0) ,'z'] = 0.005
        self.RZ          = self.foc.z                # roughness length
        self.RH          = self.foc.rh               # relative humidity
        self.albe        = self.foc.albedo           # albedo
        self.obsGST      = self.foc.tempS            # observed soil temp at 0.01 m

        self.EES     = 0.0001;
        self.TDAYS   = 3600.0 * 24.0 # 1 day in second
        
        
        
        
    def initial(self, intiST):
        
        RTT = np.zeros([NODE, 1])
        for i in range(snowN, NODE):
            RTT[i, 0] = intiST
         
        return RTT
        
        
        
        
    def SEB_RTS0(self):
        
         # RTT   = self.initial(-3)
         NTB1  = int(self.DATES[self.DATES == self.begDate1].index.tolist()[0])
         
         # Define symbolic variables
         TS    = Symbol('TS')
         
         root  = TF + self.AIRT[0] # initial root
         Lstar = - 10000             # initial Lstar
         
         
         
         RTS0 = []
         LE   = []
         H    = []
         L    = []
         RA   = []
         
         for NDAYI in range(0, len(self.dateRange2)):
                    
            dayi  = NTB1 + NDAYI
            print(self.DATES[dayi])
            
            
            
            # SNOW LAYER NUMBER AND THICKNESS
            SNOWH  = self.SND[dayi]             # snow depth
            ROSNOW = self.snowRo[dayi]          # snow density
            
            KSNOW = snowThermalCon(ROSNOW)      # snow conductivity
            CSNOW = snowThermalCap(ROSNOW)      # snow capacity
            KSNOW = KSNOW * 3600.0 * 24.0
            
            # albedo dayi
            if SNOWH > 0:
                albei = snowAlbedo(ROSNOW)
            else:
                albei = albedoG
             
            
            # snow layer [1-5]
            for i in range(0, 5):
                  
                  RKKK[i, 0] = KSNOW
                  RCCC[i, 0] = CSNOW
            
            
            TA     = TF + self.AIRT[dayi]
            tDew   = TF + self.dewP[dayi]
            QSI    = self.SRin[dayi]
            AQSI   = (1 - albei) * QSI
            wsi    = self.WS[dayi]
            PA     = self.PRE[dayi]
            roAir  = rho_a(TA, PA)
            EA     = atmosphericVaporPressure(tDew)
            QLI    = self.TRin[dayi]
            EMI    = self.emi[dayi]
            RZ     = self.RZ[dayi]
            
            # Surface energy balance
            
            # Up longwave
            Qle   = longWaveOut(EMI, TS)
            
            # ra : aerodynamic resistances
            ra    = ra_air(wsHeight, RZ, Lstar, wsi)
            
            # Sensible heat
            QH    = sensible_heat(roAir, CAir, ra, TA, TS)
            
            # Latent heat
            ES0   = atmosphericVaporPressure(TS)
            QE    = latent_heat(roAir, latSub, ra, EA, ES0, PA)
            
            # balance
            equ   = AQSI + QLI + Qle - QH - QE
            
            TS      = sympy.symbols('TS')
            f       = sympy.simplify(equ)
            f_prime = f.diff(TS)
             
            
            if __name__ == "__main__":
                def function(x):
                    return f.subs({'TS':x}) # The main function
                def derivative(x):
                    return f_prime.subs({'TS':x}) # The derivative of the main function
            
            
            root  = float(newton(function, derivative, root, 0.0001, 500))
            
            
            # Monin-Obokhov length calculation
            Ustar = ustar(wsi, wsHeight, RZ, Lstar)
            Lstar = Monin_obokhov(roAir, Ustar, TA, 
                                  QH.subs({'TS':root}), QE.subs({'TS':root}))
            
            # print(ra)            
            # print(Ustar)
            # print(Lstar)
            # print(QE.subs({'TS':root}))
            # print(QH.subs({'TS':root}))
            # print(root)
            # print(TA)
            
            RTS0.append(root - 273.15)
            LE.append(QE.subs({'TS':root}))
            H.append(QH.subs({'TS':root}))
            L.append(Lstar)
            RA.append(ra)


         df = pd.DataFrame({'date':self.dateRange2, 
                            'QE':LE, 'QH':H, 'RTS0':RTS0, 'Lstar':L,
                            'RA': RA})
         
         return df



# === SET PATH ===


dir_mod = '/Users/shengdiwang/Library/CloudStorage/OneDrive-个人/桌面/JRA/SIM'
dir_foc = '/Users/shengdiwang/Library/CloudStorage/OneDrive-个人/桌面/JRA/FOC'
dir_var = '/Users/shengdiwang/Library/CloudStorage/OneDrive-个人/桌面/JRA/VAR'


# === META FILE ===

metaf = 'meta.csv'
metaf = path.join(dir_var, metaf)


# === RUN SEB ===
# metaf:meta path, site:site name
test = SEB(metaf, site = 1)
df = test.SEB_RTS0()
df.to_csv(path.join(dir_mod, 'SEB_ch04_ori.csv'), index = False)






# df_quad = pd.read_csv(path.join(dir_mod, 'SEB_psy_CG_quad.csv'))
# df_quad['date'] = pd.to_datetime(df_quad['date'])
# mask = df_quad['date'] >= pd.to_datetime('20000101')
# mask *= (df_quad['date'] <= pd.to_datetime('20001231'))
# df_quad = df_quad[mask]


df_ori = pd.read_csv(path.join(dir_mod, 'SEB_ch04_ori.csv'))
df_ori['date'] = pd.to_datetime(df_ori['date'])
mask = df_ori['date'] >= pd.to_datetime('20001001')
mask *= (df_ori['date'] <= pd.to_datetime('20010930'))
df_ori = df_ori[mask]


# df_cg = pd.read_csv(path.join(dir_mod, 'SEB_ch04.csv'))
# df_cg['date'] = pd.date_range('20001001', '20011001', freq='3H')[:-1]
# df_cg = pd.DataFrame(df_cg).set_index('date')
# df_cg = df_cg.resample('D').mean()


test.foc['date'] = pd.to_datetime(test.foc['date']) 
mask = test.foc['date'] >= pd.to_datetime('20001001')
mask *= (test.foc['date'] <= pd.to_datetime('20011001'))
test.foc = test.foc[mask]



def gcp():
    
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.tick_params(width=0.5,length=2)
    plt.legend(loc=2)




cm = 1/2.54
fig = plt.figure(figsize=(12*cm, 8*cm))
font = {'family': 'Arial',
        'size': 4}
import matplotlib
matplotlib.rc('font', **font)



ax1 = fig.add_subplot(2,2,1)
# plt.plot(df_cg.index, df_cg.QE, 'r', lw = 0.5, label='ctrl')
# plt.plot(df_cg.date, df_quad.QE, 'k', lw = 0.5)
plt.plot(df_ori.date, df_ori.QE, 'g', lw = 0.5, label='sim')
gcp()
plt.ylabel("QE")



ax2 = fig.add_subplot(2,2,2)
# plt.plot(df_cg.index, df_cg.QH, 'r', lw = 0.5, label='ctrl')
# plt.plot(df_quad.date, df_quad.QH, lw = 0.5)
plt.plot(df_ori.date, df_ori.QH, 'g', lw = 0.5, label='sim')
gcp()
plt.ylabel("QH")



ax3 = fig.add_subplot(2,2,3)
# plt.plot(df_cg.index, df_cg.RTS0, 'r', lw = 0.5, label='ctrl')
# plt.plot(df_quad.date, df_quad.RTS0, lw = 0.5)
plt.plot(df_ori.date, df_ori.RTS0, 'g', lw = 0.5, label='sim')
gcp()
plt.ylabel("Ts")
# plt.ylim(-25, 33)



# ax4 = fig.add_subplot(2,2,4)
plt.plot(test.foc.date, test.foc.satFinal, 'k', lw = 0.5, label='Ta')
plt.plot(test.foc.date, test.foc.snd * 100, 'grey', lw = 0.5, label='snd')
plt.plot(test.foc.date, test.foc.sw, 'r', lw = 0.5, label='sw')
plt.plot(test.foc.date, test.foc.lwin, 'y', lw = 0.5, label='sw')
gcp()
plt.ylabel("Ta")
# plt.ylim(-25, 33)





plt.savefig('ch04.png', dpi = 400, bbox_inches='tight')






# plt.plot(df_ori.date, df_ori.Lstar, 'g', lw = 0.5, label='sim')


# plt.ylim(-1, 1)













