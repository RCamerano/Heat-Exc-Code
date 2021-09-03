# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 09:50:44 2021

@author: Paolo
"""
from CoolProp.CoolProp import PropsSI
import numpy as np

def heat_transfer_factor_tube(Re,HE_ratio):
    if Re <= 2000:
        # Il fattore j_h si ottiene con una relazione del tipo j_h = A*Re**B, con B costante, ed A funzione di HE_ratio
        A = 1.9*HE_ratio**-0.361
        B = -0.64
        j_h = A*Re**B
    elif Re >= 20000:
        j_h = 0.0281*Re**-0.2
    else:
        A = 1.9*HE_ratio**-0.361
        B = -0.64
        j_h_1 = A*2000**B
        j_h_2 = 0.0281*20000**-0.2
        
        # Interpolazione lineare
        m = (j_h_1 - j_h_2)/(2000 - 20000)
        b = (2000*j_h_2 - 20000*j_h_1)/(2000 - 20000)
        j_h = m*Re + b
        
    return j_h
    

def heat_transfer_factor_shell(Re,baffle_cut):
    if baffle_cut == 15:
        j_h = 0.48*Re**0.464
    elif baffle_cut == 20:
        j_h = 0.49*Re**0.464
    elif baffle_cut == 25:
        j_h = 0.4181*Re**0.464
    elif baffle_cut == 30:
        j_h = 0.4058*Re**0.4675
    elif baffle_cut == 35:
        j_h = 0.3936*Re**0.471
    elif baffle_cut == 40:
        j_h = 0.3848*Re**0.473
    elif baffle_cut == 45:
        j_h = 0.376*Re**0.475
        
    return j_h
    
# La funzione 'heat_exchange_coefficient_condensation' calcola il coefficiente di scambio per fluido che riscalda/raffredda     
def heat_transfer_coefficient_1ph(T_inlet,T_outlet,P_inlet,P_outlet,fluid,Re_inlet,Re_outlet,j_h,D_hydraulic,idx):
    
    # Dittus-Boelter correlation - Tube side
    viscosity_mean = PropsSI('V', 'T', np.mean([T_inlet[idx],T_outlet[idx]]), 'Q', 0, fluid[idx])
    T_wall = np.mean([np.mean([T_inlet[0],T_outlet[0]]),np.mean([T_inlet[1],T_outlet[1]])])
    viscosity_wall = PropsSI('V', 'T', T_wall, 'Q', 0, fluid[idx])
    
    k = PropsSI('L', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'H', np.mean([T_inlet[idx],T_outlet[idx]]), fluid[idx])
    C_p = PropsSI('Cpmass', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'H', np.mean([T_inlet[idx],T_outlet[idx]]), fluid[idx])
    
    Pr = viscosity_mean*C_p/k
    
    Nu = j_h[idx]*(np.mean([Re_inlet[idx],Re_outlet[idx]])**0.8)*(Pr[idx]**0.33)*(viscosity_mean/viscosity_wall)**0.14
    alfa = Nu*k/D_hydraulic
    
    return alfa
    
# La funzione 'heat_exchange_coefficient_condensation' calcola il coefficiente di scambio per un fluido bifase che condensa   
def heat_exchange_coefficient_condensation(T_inlet,T_outlet,P_inlet,P_outlet,x_inlet,x_outlet,D_hydraulic,fluid,idx,M):                                         
    
    x = np.linspace(x_inlet,x_outlet,10)
    alfa_discrete = np.zeros(len(x))
    
    cp_l = PropsSI('Cpmass', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx])
    k_l = PropsSI('L', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx])
    viscosity_l = PropsSI('V', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx])
    
    Pr_l = viscosity_l*cp_l/k_l
    rho_l = PropsSI('D', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx])
    velocity_l = M[idx]/(rho_l*np.pi*0.25*D_hydraulic**2)
    Re_l = rho_l*velocity_l*D_hydraulic/viscosity_l
    
    alfa_l = 0.023*(Re_l**0.4)*(Pr_l**0.4)*k_l/D_hydraulic
    p_star = PropsSI('D', 'T', np.mean([T_inlet[idx],T_outlet[idx]]), 'Q', 0, fluid[idx])/PropsSI('PCRIT', fluid[idx])
    
    for i in range(len(x)):
        alfa_discrete[i] = alfa_l*((1-x[i])**0.8 + (((1-x[i])**0.04)*3.8*x[i]**0.76)/p_star**0.38)
        
    alfa = np.mean(alfa_discrete)
    
    return alfa

    