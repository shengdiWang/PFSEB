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




def latentWater(T):
    
    # latent heat of sublimation [J kg-1]
    L = 1E3 * (2500.8 - 2.36 * (T - 273.15))
    
    return L


def latentIce():
        
    return 2835000


def snowThermalCap(roSnow):
    # snow thermal capacity from HTESSEL

    roIce = 920
    cIce  = 2.05
    snCap = roSnow * (cIce / roIce)     
    
    return snCap


def snowThermalCon(roSnow):
    # snow thermal conductivity from HTESSEL

    roIce = 920
    kIce  = 2.29
    snTC  = kIce * (roSnow / roIce) ** 1.88 # snow conductivity
    
    return snTC


def snowAlbedo(roSnow):
    # snow albedo from Eq. 31

    snalbedo = 1 - 0.247 * math.sqrt(0.16 + 110 * (roSnow / 1000) ** 4)

    return snalbedo


def atmosphericVaporPressure(t):
    # Atmospheric vapor pressure as Eq. 3
    # avp = 10 ** (11.40 - 2353.0 / t)
    
    avp = 611.2 * sympy.exp(17.62 * (t - 273.15) / (243.12 + (t - 273.15)))
    
    return avp


def surfaceVaporPressure(t):
    # Atmospheric vapor pressure as Eq. 3
    # avp = 10 ** (11.40 - 2353.0 / t)
    
    if t > TF:
        avp = 611.2 * sympy.exp(17.62 * (t - 273.15) / (243.12 + (t - 273.15)))
    else:
        avp = 611.2 * sympy.exp(22.46 * (t - 273.15) / (272.61 + (t - 273.15)))
        
    return avp


def rho_a(Tz, p):

    R_a = 287.058    # specific gas constant of air [ J/(kg K) ]
    res = p / (R_a * Tz) # air density [kg m^(-3)]

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


def longWaveOut(emiV, Ts0, Lin):
    # Emitted longwave radiation as Eq. 4

    SIGMA = 5.67e-8 # Stefan-Boltzmann constant [W m-2 K-4]
    LWout = emiV * SIGMA * Ts0 ** 4 + (1 - emiV) * Lin
    
    return LWout


# stable
def psi_H(zeta1, zeta2):

    if zeta1 <= 0:
        
        # def funH(x):
        #     return (1 - 0.95 * (1 - 11.6 * x) ** (-0.5)) / x
        # res = integrate.quad(funH, zeta2, zeta1)[0]
        
        res = 1.9 * cmath.atanh((1 - 11.6 * zeta1) ** 0.5) + cmath.log(zeta1) \
              - (1.9 * cmath.atanh((1 - 11.6 * zeta2) ** 0.5) + cmath.log(zeta2))
    
    else:
        
        # def funH(x):
        #     return (1 - 0.95 - 7.8 * x) / x
        # res = integrate.quad(funH, zeta2, zeta1)[0]

        res = ((-5 + 5 ** 0.5) * cmath.log(-3 + 5 ** 0.5 - 2 * zeta1) - (5 + 5 ** 0.5) * cmath.log(3 + 5 ** 0.5 + 2 * zeta1))/2 \
              - (((-5 + 5 ** 0.5) * cmath.log(-3 + 5 ** 0.5 - 2 * zeta2) - (5 + 5 ** 0.5) * cmath.log(3 + 5 ** 0.5 + 2 * zeta2))/2)
    
    return res


# stable
def psi_M(zeta1, zeta2):

    if zeta1 <= 0:
        
        # def funM(x):
        #     return (1 - (1 - 19 * x) ** (-0.25)) / x
        # res = integrate.quad(funM, zeta2, zeta1)[0]
        
        res = -2 * atan((1 - 19 * zeta1) ** (1/4)) + 2 * math.log(1 + (1 - 19 * zeta1) ** (1/4)) + math.log(1 + (1 - 19 * zeta1) ** 0.5) \
              - (-2 * atan((1 - 19 * zeta2) ** (1/4)) + 2 * math.log(1 + (1 - 19 * zeta2) ** (1/4)) + math.log(1 + (1 - 19 * zeta2) ** 0.5))
    
    else:
        
        # def funM(x):
        #     return (1 - (1 + 6 * x)) / x
        # res = integrate.quad(funM, zeta2, zeta1)[0]

        res = -19.5 * (1 + zeta1) ** (1/3) - 7.5367 * atan(0.57735 - 1.72489 * (1 + zeta1) ** (1/3)) + 4.35131 * log(3+4.4814 * (1+zeta1) ** (1/3)) - 2.17566 * log(3 - 4.4814 * (1 + zeta1) ** (1/3) + 6.69433 * (1 + zeta1) ** (2/3)) \
              - (-19.5 * (1 + zeta2) ** (1/3) - 7.5367 * atan(0.57735 - 1.72489 * (1 + zeta2) ** (1/3)) + 4.35131 * log(3+4.4814 * (1+zeta2) ** (1/3)) - 2.17566 * log(3 - 4.4814 * (1 + zeta2) ** (1/3) + 6.69433 * (1 + zeta2) ** (2/3)))
    
    return res


# Sensible heat flux    
def sensible_heat(roAir, CAir, TA, TS, wsi, RZ, Lstar):
    
    # QH = - roAir * CAir / ra * (TA - TS)
    QH = - roAir * CAir * VONK ** 2 * wsi / (math.log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar)) * (TA - TS) / (math.log(wsHeight / RZ)- psi_H(wsHeight / Lstar, RZ / Lstar))
    
    return QH


def latent_heat_w(roAir, latSub, SH, ES0, PA, wsi, RZ, Lstar):
    
    # QE = - roAir * latSub / (ra + 50) * (0.622 * (EA - ES0) / PA)
    QE = - roAir * latSub * VONK ** 2 * wsi / (math.log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar)) * (SH - 0.622 * ES0 / PA) \
        / (math.log(wsHeight / RZ)- psi_H(wsHeight / Lstar, RZ / Lstar)) + 50 * VONK ** 2 * wsi / (math.log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar))
         
    return QE


def latent_heat_i(roAir, latSub, SH, ES0, PA, wsi, RZ, Lstar):
    
    # QE = - roAir * latSub / (ra + 50) * (0.622 * (EA - ES0) / PA)
    QE = - roAir * latentIce() * VONK ** 2 * wsi / (math.log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar)) * (SH - 0.622 * ES0 / PA) \
        / (math.log(wsHeight / RZ)- psi_H(wsHeight / Lstar, RZ / Lstar))
    
    return QE


def Monin_obokhov(roAir, Ustar, TA, t, QH, QE):
    
    
    # if t > TF:
    #     latSub = latentWater(TA)
    #     Lstar  = - roAir * CAir * Ustar ** 3 * TA / (VONK * GRAVIT * (QH + 0.61*CAir*TA*QE/latSub))
    # else:
    #     latSub = latentIce()
    #     Lstar  = - roAir * CAir * Ustar ** 3 * TA / (VONK * GRAVIT * (QH + 0.61*CAir*TA*QE/latSub))
    
    latSub = latentWater(TA)
    Lstar  = - roAir * CAir * Ustar ** 3 * TA / (VONK * GRAVIT * (QH + 0.61*CAir*TA*QE/latSub))
    
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

    return Lstar


def ustar(wsi, wsHeight, RZ, Lstar): 
    
    zeta1 = wsHeight / Lstar
    zeta2 = RZ / Lstar
    
    return (wsi * VONK)/(math.log(wsHeight / RZ) - psi_M(zeta1, zeta2))






# === PROFILE ===

# 5-layer snow, 16 layers for 0-1 m, 190 layers for 1-20 m [0.1 m interval]


soil1m    = [0, 3, 7, 11, 15, 22, 29, 35, 44, 50, 58, 66, 75, 84, 93, 100]
soil1m    = [i / 100 for i in soil1m]
soilTck   = 0.1 # thickness for deep soil (> 1 m)


snowN     = 5   # maximum snow layer
SLTck1    = 1 # soil layer 1 thickness
SLTck2    = 9 # soil layer 2 thickness
soilDepth = 10  # total soil depth
NODE      = int(snowN + len(soil1m) + (soilDepth - 1) / soilTck) # layer number
intiST    = -3 # initial soil temperature [C]


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
# roAir = 1.225  # air density [kg m-3]
CAir  = 1005.7  # air thermal capacity [J m-3 C-1]

# WATER
roWater = 1000.0 # water density

# ICE
roIce = 920   # ice density
cIce  = 2.05  # ice volumetric heat capacity [MJ m-3 K-1]
kIce  = 2.29  # thermal conductivity of ice  [W m-1 K-1]

GRAVIT   = 9.81   # gravitational acceleration
VONK     = 0.4     # Von Karman's constant
TF       = 273.15  # unit C to K
wsHeight = 2.0     # REFERENCE HEIGHT for wind speed measure
albedoG  = 0.2     # snow-free albedo
SIGMA    = 5.67e-8




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
        
        self.begDate2   = pd.to_datetime(str(self.begYr) + '1001')     # start date for run
        self.endDate2   = pd.to_datetime(str(self.endYr) + '1001')     # end date for run
        self.dateRange2 = pd.date_range(self.begDate2, self.endDate2, freq='3H')[:-1]
        
        
        
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
        self.foc['emi']  = 0.97
        self.foc.loc[(self.foc.snd > 0) ,'emi'] = 0.99
        self.emi         = self.foc.emi              # emissivity
        self.foc['z']    = 0.001
        self.foc.loc[(self.foc.snd > 0) ,'z'] = 0.0005
        self.RZ          = self.foc.z                # roughness length
        self.RH          = self.foc.rh               # relative humidity
        self.albe        = self.foc.albedo           # albedo
        self.obsGST      = self.foc.tempS            # observed soil temp at 0.01 m
        self.SH          = self.foc.q                # q
        self.TS          = self.foc.tempS            # tempS
        
        self.EES     = 0.0001
        self.TDAYS   = 3600.0 * 24.0 # 1 day in second
        
        
        
        
    def initial(self, intiST):
        
        RTT = np.zeros([NODE, 1])
        for i in range(snowN, NODE):
            RTT[i, 0] = intiST
         
        return RTT
        
        
        
        
    def SEB_RTS0(self):
        
         # RTT   = self.initial(-3)
         NTB1  = int(self.DATES[self.DATES == self.begDate2].index.tolist()[0])
         
         # Define symbolic variables
         # TS    = Symbol('TS')
         
         # root  = TF + self.AIRT[0]   # initial root
         # root  = TF + -0.01
         
         Lstar = -100000             # initial Lstar
         QE    = 0
         QH    = 0
         Ustar = 10
         
         
         RTS0 = []
         LE   = []
         H    = []
         L    = []
         # RA   = []
         
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
            SH     = self.SH[dayi]
            root   = TF + self.TS[dayi]
            
            # Surface energy balance
            
            # Up longwave
            Qle   = longWaveOut(EMI, root, QLI)
            
            # ra : aerodynamic resistances
            # ra    = ra_air(wsHeight, RZ, Lstar, wsi)
            
            # Sensible heat
            QH    = sensible_heat(roAir, CAir, TA, root, wsi, RZ, Lstar)
            
            # Latent heat
            ES0   = surfaceVaporPressure(root)
            
            
            if root > TF:
                latSub = latentWater(root)
                QE = latent_heat_w(roAir, latSub, SH, ES0, PA, wsi, RZ, Lstar)
            else:
                latSub = latentIce()
                QE = latent_heat_i(roAir, latSub, SH, ES0, PA, wsi, RZ, Lstar)
            
            
            # balance
            equ     = AQSI + QLI - Qle - QH - QE
            
            # TS      = sympy.symbols('TS')
            # f       = sympy.simplify(equ)
            # f_prime = f.diff(TS)
             
            
            # if __name__ == "__main__":
            #     def function(x):
            #         return f.subs({'TS':x}) # The main function
            #     def derivative(x):
            #         return f_prime.subs({'TS':x}) # The derivative of the main function
            
            
            # root  = float(newton(function, derivative, root, self.EES, 500))
            
            
            # Monin-Obokhov length calculation
            Ustar = ustar(wsi, wsHeight, RZ, Lstar)
            # Lstar = Monin_obokhov(roAir, Ustar, TA, root,
            #                       QH.subs({'TS':root}), QE.subs({'TS':root}))
            Lstar = Monin_obokhov(roAir, Ustar, TA, root, QH, QE)
            
            # print(ra)            
            # print(Ustar)
            # print(Lstar)
            # print(QE.subs({'TS':root}))
            # print(QH.subs({'TS':root}))
            # print(root)

            
            RTS0.append(root - 273.15)
            # LE.append(QE.subs({'TS':root}))
            # H.append(QH.subs({'TS':root}))
            LE.append(QE)
            H.append(QH)
            L.append(Lstar)
            # RA.append(ra)


         df = pd.DataFrame({'date':self.dateRange2, 
                            'QE':LE, 'QH':H, 'RTS0':RTS0, 'Lstar':L})
         
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


test = SEB(metaf, site = 6)
df = test.SEB_RTS0()
df.to_csv(path.join(dir_mod, 'SEB_cg3_ori.csv'), index = False)













df_ori = pd.read_csv(path.join(dir_mod, 'SEB_cg3_ori.csv'))
df_ori['date'] = pd.to_datetime(df_ori['date'])
mask = df_ori['date'] >= pd.to_datetime('20001001')
mask *= (df_ori['date'] <= pd.to_datetime('20011001'))
df_ori = df_ori[mask]


df_cg = pd.read_csv(path.join(dir_mod, 'SEB_CG.csv'))
df_cg['date'] = pd.date_range('20001001', '20011001', freq='3H')[:-1]
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

plt.plot(df_ori.date, df_ori.QE, 'g', lw = 0.5, label='sim')
plt.plot(df_cg.date, df_cg.QE, 'r', lw = 0.5, label='ctrl')

gcp()
plt.ylabel("QE")
plt.ylim(-50, 350)

df_cg.QE.mean()
df_ori.QE.mean()


ax2 = fig.add_subplot(2,2,2)

plt.plot(df_ori.date, df_ori.QH, 'g', lw = 0.5, label='sim')
plt.plot(df_cg.date, df_cg.QH, 'r', lw = 0.5, label='ctrl')

gcp()
plt.ylabel("QH")
plt.ylim(-50, 350)

df_cg.QH.mean()
df_ori.QH.mean()


ax3 = fig.add_subplot(2,2,3)

plt.plot(df_cg.date, df_cg.RTS0, 'r', lw = 0.5, label='ctrl')
plt.plot(df_ori.date, df_ori.RTS0, 'g', lw = 0.5, label='sim')

gcp()
plt.ylabel("Ts")
plt.ylim(-50, 33)

df_cg.RTS0.mean()
df_ori.RTS0.mean()


# ax4 = fig.add_subplot(2,2,4)
# plt.plot(test.foc.date, test.foc.satFinal, 'k', lw = 0.5, label='Ta')
# plt.plot(test.foc.date, test.foc.snd * 100, 'grey', lw = 0.5, label='snd')
# plt.plot(test.foc.date, test.foc.sw, 'r', lw = 0.5, label='sw')
# plt.plot(test.foc.date, test.foc.lwin, 'y', lw = 0.5, label='sw')
gcp()
plt.ylabel("Ta")
plt.ylim(-50, 33)



# ax4 = fig.add_subplot(2,2,4)
# # plt.plot(df_ori.date, df_ori.Lstar, 'grey', lw = 0.5, label='sim')
# # plt.ylim(-40, 40)
# gcp()
# plt.ylabel("Lstar")


# plt.plot(df_ori.date, np.array(df_ori.RTS0) - np.array(df_cg.RTS0), 'g', lw = 0.5)
# plt.axhline(0, color='k', lw = 0.2, ls='--')



plt.savefig('cg3.png', dpi = 400, bbox_inches='tight')










