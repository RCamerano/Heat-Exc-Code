# -*- coding: utf-8 -*-
"""
Created on Fri May 27 16:24:16 2022

@author: Utente
"""
import numpy as np
import Thermeprop
import CoolProp.CoolProp as CP

fluid = 10068
mixture = 1
fraction = 1
p = 2
T_inlet = Thermeprop.getThermProp(mixture, 'T', 'P', p, 'Q', 0, fluid, fraction)
T_outlet = Thermeprop.getThermProp(mixture, 'T', 'P', p, 'Q', 0, fluid, fraction)
P_inlet = p
P_outlet = p
x_inlet = 0.99
x_outlet = 0.01
D_hydraulic = 0.05
M = 0.05

x = np.linspace(x_inlet,x_outlet,10)
alfa_discrete = np.zeros(len(x))

cp_l = Thermeprop.getThermProp(mixture, 'C', 'P', np.mean([P_inlet,P_outlet]), 'Q', 0, fluid, fraction)
k_l = Thermeprop.getThermProp(mixture, 'L', 'P', np.mean([P_inlet,P_outlet]), 'Q', 0, fluid, fraction)
viscosity_l = Thermeprop.getThermProp(mixture, 'V', 'P', np.mean([P_inlet,P_outlet]), 'Q', 0, fluid, fraction)

Pr_l = viscosity_l*cp_l/k_l
rho_l = Thermeprop.getThermProp(mixture, 'D', 'P', np.mean([P_inlet,P_outlet]), 'Q', 0, fluid, fraction)
velocity_l = M/(rho_l*np.pi*0.25*D_hydraulic**2)
Re_l = rho_l*velocity_l*D_hydraulic/viscosity_l

alfa_l = 0.023*(Re_l**0.4)*(Pr_l**0.4)*k_l/D_hydraulic
p_crit = CP.PropsSI("Pcrit","",0,"",0,"R1233zd(E)")/100000
pr = p/p_crit

for i in range(len(x)):
    alfa_discrete[i] = alfa_l*((1-x[i])**0.8 + (((1-x[i])**0.04)*3.8*x[i]**0.76)/pr**0.38)
    
alfa = np.mean(alfa_discrete)