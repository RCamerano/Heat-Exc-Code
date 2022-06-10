# -*- coding: utf-8 -*-
"""
Created on Fri May 27 14:50:11 2022

@author: Utente
"""
import numpy as np
from scipy.optimize import fsolve
import Thermeprop
import CoolProp.CoolProp as CP
import pandas as pd

h = []
x = []

mixture = 1
fraction = 1
fluid = 10068
p = 1 # bar
M = 0.005 # kg/s
ID = 0.0088 # m
A_cross = np.pi*(ID**2) / 4
G = M/A_cross

with open('Cond_tubes' + '.txt', 'w') as f:
    for i in range (50,55,5):
        Q = i/100
        T_sat = Thermeprop.getThermProp(mixture, 'T', 'P', p, 'Q', Q, fluid, fraction)
        rho_flow = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', Q, fluid, fraction)
        viscosity_flow = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', Q, fluid, fraction)
        velocity = M / rho_flow / A_cross
        # Re = rho_flow * velocity * ID / viscosity_flow
        
        # Calcolo delle proprietà del fluido
        rho_L = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 0, fluid, fraction)
        rho_G = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 1, fluid, fraction)
        viscosity_L = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 0, fluid, fraction)
        viscosity_G = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 1, fluid, fraction)
        visc_ratio = viscosity_L / viscosity_G
        cp_L = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 0, fluid, fraction)
        k_L = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 0, fluid, fraction)
        Pr_L = viscosity_L*cp_L/k_L
        
        # Calcolo slip ratio 'S' e spessore strato liquido 'sigma'
        
        X_tt = ((1 - Q)/Q)**0.875 * (rho_G/rho_L)**0.5*(viscosity_L/viscosity_G)**0.125
        Ft = ( (Q**3 * G**2) / (rho_G**2 * 9.81 * ID * (1-Q)) )**0.5
        ## Calcolo dell'alfa: Void Fraction
       
        # Modello di Graham | Migliori predizioni sperimentali per refrigeranti con tubi di ID > 2 mm
        alfa = ( 1 + X_tt + Ft**-1)**-0.321
        
        def theta_stratified(vars, *data):
            A_L, R = data
            theta_strat = vars
            return A_L - (0.5*R**2)*(2*np.pi - theta_strat - np.sin(2*np.pi - theta_strat))
        
        A_L = A_cross*(1 - alfa)
        data =  (A_L, ID/2)
        theta_strat = fsolve(theta_stratified, np.pi, args = data)[0]
        beta = np.pi - theta_strat/2
        
        # Calcolo sigma con totale stratificazione
        # sigma = ID/2 * (1 - np.cos(beta))
        
        # Calcolo sigma con fluido disposto come settore di corona circolare
        r = np.sqrt(ID**2 / 4 - A_L / beta)
        sigma = ID/2 - r
        
        # Calcolo slip ratio S
        NUM = 0.333*alfa**2 + 1.333*visc_ratio**(-1) * alfa*(1-alfa) + visc_ratio**(-1) * (1-alfa)**2
        DEN = visc_ratio * alfa**2 + 1.333 * visc_ratio * alfa*(1-alfa) + 0.333*(1-alfa)**2
        S = visc_ratio**2 *alfa/(1-alfa) * NUM / DEN
        
        # Se il flusso di massa specifico è inferiore a 200 kg/m^2s viene utilizzata la relazione di
        # Akers per il coefficente di scambio, altrimenti la relazione di Shah. 
        if G < 200:
            G_e = G*((1 - Q) + Q*(rho_L/rho_G)**0.5)
            Re_e = ID*G_e/viscosity_L
            if Re_e > 50000:
                C = 0.0265
                n = 0.8
            else:
                C = 5.03
                n = 1/3
            h.append( C*Re_e**n*Pr_L**0.333*k_L/ID )
        else:
            p_crit = CP.PropsSI("Pcrit","",0,"",0,"R1233zd(E)")/100000
            pr = p/p_crit
            
            # Conti con h_LT
            Re_LT = G * ID / viscosity_L
            h_LT = 0.023*(Re_LT**0.8)*(Pr_L**0.4)
            # h.append(h_LT*( (1 - Q)**0.8 + 3.8*Q**0.76*(1 - Q)**0.04/(pr)**0.38 ) )
            
            # Conti separati: uso h_LS
            # Re_LS = G * (1-Q) * ID / viscosity_L
            # h_LS = 0.023*(Re_LS**0.8)*(Pr_L**0.4)
            # Z = (1/Q - 1)**0.8 * pr**0.4
            # h.append(h_LS * (1 + 3.8/(Z**0.95)))
            
            
            # 06.06.2022
            
            # Calcolo h_LT tramite versione di Shah 2009
            n = 0.0058 + 0.557*pr
            
            # h_LT = 0.023*(Re_LT**0.8)*(Pr_L**0.4)*k_L/sigma*(rho_L/rho_G)**0.333
            h_I = h_LT * (viscosity_L/14/viscosity_G)**n * ( (1 - Q)**0.8 + (3.8* Q**0.76 * (1 - Q)**0.04 / (pr)**0.38) )
            
            Re_LS = G * (1-Q) * ID / viscosity_L
            h_NU = 1.32 * Re_LS**-0.333 * ( ( rho_L *(rho_L - rho_G) * 9.81 * k_L**3) / viscosity_L**2 )**0.333
            
            h.append(h_I + h_NU)
                    
        x.append(Q)
            
Cond_tube = pd.DataFrame({'Quality': x,
                   'Alfa': h})

writer = pd.ExcelWriter(r'C:\\Users\\Utente\\Desktop\\Cond_tube.xlsx', engine='xlsxwriter')
Cond_tube.to_excel(writer, sheet_name='Sheet1', index=False)
writer.save()