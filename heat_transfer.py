# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 09:50:44 2021

@author: Paolo
"""

import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import Thermeprop
from CoolProp.CoolProp import PropsSI

#%%  HEAT TRANSFER FACTOR heat_transfer_coefficient_1ph

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
        j_h = 0.48*Re**-0.464
    elif baffle_cut == 20:
        j_h = 0.49*Re**-0.464
    elif baffle_cut == 25:
        j_h = 0.4181*Re**-0.464
    elif baffle_cut == 30:
        j_h = 0.4058*Re**-0.4675
    elif baffle_cut == 35:
        j_h = 0.3936*Re**-0.471
    elif baffle_cut == 40:
        j_h = 0.3848*Re**-0.473
    elif baffle_cut == 45:
        j_h = 0.376*Re**-0.475
        
    return j_h

#%% HEAT TRANSFER 1ph PRIMA VERSIONE

# La funzione 'heat_exchange_coefficient_condensation' calcola il coefficiente di scambio per fluido che riscalda/raffredda     
def heat_transfer_coefficient_1ph(mixture,fraction,T,p,fluid,Re,j_h,D_hydraulic,idx,L):
    
    if Re[idx] > 10000:
        # Dittus-Boelter correlation - Tube side
        viscosity_mean = Thermeprop.getThermProp(mixture, 'V', 'T', T[idx], 'Q', 0, fluid[idx], fraction)
        T_wall = np.mean([T[0],T[1]])
        viscosity_wall = Thermeprop.getThermProp(mixture, 'V', 'T', T_wall, 'Q', 0, fluid[idx], fraction)
        
        k = Thermeprop.getThermProp(mixture, 'L', 'T', T[idx], 'P', p[idx], fluid[idx], fraction)
        C_p = Thermeprop.getThermProp(mixture, 'C', 'T', T[idx], 'P', p[idx], fluid[idx], fraction)
        
        Pr = viscosity_mean*C_p/k
       
        # Nu = j_h[idx]*(Re[idx]**0.8)*(Pr**0.33)*(viscosity_mean/viscosity_wall)**0.14
        Nu = 0.023*(Re[idx]**0.8)*(Pr**0.365)*(viscosity_mean/viscosity_wall)**0.14
     
        alfa = Nu*k/D_hydraulic
        
    elif Re[idx] < 2100:
        Nu = 3.66
        k = Thermeprop.getThermProp(mixture, 'L', 'T', T[idx], 'P', p[idx], fluid[idx], fraction)
        # aggiunta k, verifica se corretta
        alfa = Nu*k/D_hydraulic
    
    else:
        viscosity_mean = Thermeprop.getThermProp(mixture, 'V', 'T', T[idx], 'Q', 0, fluid[idx], fraction)
        T_wall = np.mean([T[0],T[1]])
        viscosity_wall = Thermeprop.getThermProp(mixture, 'V', 'T', T_wall, 'Q', 0, fluid[idx], fraction)
        
        k = Thermeprop.getThermProp(mixture, 'L', 'T', T[idx], 'P', p[idx], fluid[idx], fraction)
        C_p = Thermeprop.getThermProp(mixture, 'C', 'T', T[idx], 'P', p[idx], fluid[idx], fraction)
        
        Pr = viscosity_mean*C_p/k
        
        Nu_10000 = j_h[idx]*(10000**0.8)*(Pr**0.33)*(viscosity_mean/viscosity_wall)**0.14
        
        Nu_2100 = 3.66
        
        m = (Nu_10000 - Nu_2100) / (10000 - 2100)
        q = Nu_2100 - m*2100
        Nu = m*Re[idx] + q
        
        alfa = Nu*k/D_hydraulic
    
    return alfa

#%% SHELL

# La funzione 'heat_transfer_coefficient_1ph_Shell' calcola il coefficiente di scambio per un fluido monofase
# che scalda/raffredda nello shell.
def heat_transfer_coefficient_1ph_Shell(mixture,fraction,T,T_tube,P,velocity,viscosity,fluid,rho,Deq,OD,n_tubi,D_shell_int,baffle_cut,B,D_bundle,baffle_hole_d,pitch,tube_disposition):
    
    print('Calcolo coefficente di scambio termico mono-fase lato shell in corso... \n')
    # Parametri termodinamici per il calcolo dello scambio termico sono calcolati
    T_wall = np.mean([T,T_tube])
    C_p_w = Thermeprop.getThermProp(mixture, 'C', 'T', T_wall, 'P', P, fluid, fraction)
    k_w = Thermeprop.getThermProp(mixture, 'L', 'T', T_wall, 'P', P, fluid, fraction)
    viscosity_w = Thermeprop.getThermProp(mixture, 'V', 'T', T_wall, 'P', P, fluid, fraction)
    Pr_w = viscosity_w*C_p_w/k_w
    Cp = Thermeprop.getThermProp(mixture, 'C', 'T', T, 'P', P, fluid, fraction)
    k = Thermeprop.getThermProp(mixture, 'L', 'T', T, 'P', P, fluid, fraction)
    Pr = viscosity*Cp/k
    gamma = Pr/Pr_w
    Re_d = rho*velocity*OD/viscosity
    
    # Calcolo del numero di Nusselt sulla base del regime di moto all'interno dello shell.
    if Re_d < 100:
        Nu_d = 0.9*Re_d**0.4*Pr**0.36*gamma**0.25
    elif Re_d > 100 and Re_d < 1000:
        Nu_d = 0.52*Re_d**0.5*Pr**0.36*gamma**0.25
    elif Re_d > 1000 and Re_d < 2E5:
        Nu_d = 0.27*Re_d**0.63*Pr**0.36*gamma**0.25
    elif Re_d > 2E5 and Re_d < 2E6:
        Nu_d = 0.033*Re_d**0.8*Pr**0.36*gamma**0.25
    else:
        print('Reynolds number Beyond the range of validity of the heat exchanger model \n')
    
    if n_tubi < 16:
        print('WARNING: numero di tubi inferiore a 16. Il modello utilizzato potrebbe non essere valido. \n')
    
    # # Funzione per risolvere numericamente l'equazione dell'altezza del segmento circolare, nota la sua area.
    # def segmento_circolare(vars, *data):
    #     A, r = data
    #     l_c = vars
    #     return A - r**2*np.arccos(1 - l_c/r) + (r - l_c)*(r**2 - (r - l_c)**2)**0.5
    
    # # the distance from the baffle tip to the shell inside diameter [m]
    # data =  (0.01*baffle_cut*np.pi*D_shell_int**2/4, D_shell_int/2)
    # l_c = fsolve(segmento_circolare, D_shell_int/2, args = data)[0]
    
    l_c = D_shell_int*0.01*baffle_cut
    
    # JC is the correction factor taking into consideration the baffle cut and baffle spacing and 
    # is the average for the entire exchanger:
    phi = (D_shell_int - 2*l_c)/D_bundle
    F_c = (np.pi + phi*np.sin(np.arccos(phi)) - 2*np.arccos(phi))/np.pi
    J_C = 0.55 + 0.72*F_c
    
    # JL is the correction factor taking into consideration the tube-to-baffle leakages and the
    # baffle-to-shell leakages:
    D_baffle = 0.97*D_shell_int #(D_shell_int + D_bundle)/2
    delta_sb = (D_shell_int - D_baffle)/2
    theta_1 = np.arcsin(1 - 2*l_c/D_shell_int)
    A_sb = 0.5*(np.pi + theta_1)*D_shell_int*delta_sb
    
    C_1 = D_shell_int - D_bundle
    delta_tb = baffle_hole_d - OD
    theta_3 = 2*np.arccos((D_shell_int - 2*l_c)/(D_shell_int - C_1))
    F_w = (theta_3 - np.sin(theta_3))/(2*np.pi)
    A_tb = (np.pi*OD*(1 - F_w)*n_tubi*delta_tb)/2
    
    theta_2 = 2*np.arccos((D_shell_int/2 - l_c)/(D_shell_int/2))
    A_wg = D_shell_int**2*(theta_2 - np.sin(theta_2))/8
    n_tw = F_w*n_tubi
    A_wt = 0.25*np.pi*n_tw*OD**2
    A_w = A_wg - A_wt
    
    r_a = A_sb/(A_sb + A_tb)
    r_b = (A_sb + A_tb)/A_w
    J_L = 0.44*(1 - r_a) + (1 - 0.044*(1 - r_a))*np.exp(-2.2*r_b)
    
    # JB is the correction factor taking into consideration the bundle and partition bypass effects,
    # i.e. the portion of the shell side fluid that can not take partto the heat transfer with the tube bundle (bypass):
    S_b =B*(D_bundle - OD)
    S_m = B*(D_shell_int - D_bundle + (D_bundle - OD)/pitch*(pitch - OD) )
    if Re_d < 100:
        C_j = 1.35
    else:
        C_j = 1.25
    if tube_disposition == 't':
        pitch = 0.5*pitch*3**0.5
    N_c = D_shell_int*(1 - 2*0.01*baffle_cut)/pitch
    N_ss = 1 # Number of sealing strips
    r_ss = N_ss/N_c
    J_B = 1
    if r_ss < 0.5:
        J_B = np.exp( (-C_j*S_b/S_m)*(1 - (2*r_ss)**0.33) )
    
    # JS is the correction factor that accounts for variations in baffle spacing at the inlet and 
    # outlet sections:
    # JS è uguale ad 1 se la distanza tra l'ingresso e il primo baffle e la distanza tra l'uscita
    # e l'ultimo baffle sono uguali alla distanza tra tutti i baffle.
    J_S = 1
    
    # JR is the correction factor that accounts for the temperature gradient when the shell-side 
    # fluid is in laminar flow:
    Nr_cc = (D_shell_int - 2*l_c)/pitch
    Re_s = rho*velocity*D_shell_int/viscosity
    if Re_s > 100:
        J_R = 1
    elif Re_s < 20:
        J_R = (10/Nr_cc)**0.18
    else:
        J_R_100 = 1
        J_R_20 = (10/Nr_cc)**0.18
        m = (J_R_100 - J_R_20) / (100 - 20)
        q = J_R_20 - m*20
        J_R = m*Re_s + q
        
    h = Nu_d*k/OD
    h = J_C*J_L*J_B*J_S*J_R*h
    
    return h

# La funzione 'heat_transfer_coefficient_condensation_Shell' calcola il coefficiente di scambio
# per un fluido bifase che evapora nello shell utilizando la relazione di Shekriladze and Gomelauri  
def heat_transfer_coefficient_condensation_Shell(mixture, fraction, p,Q,T,fluid,velocity,T_wall,Re,D_eq):
    k_L = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 0, fluid, fraction)
    rho_L = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 0, fluid, fraction)
    rho_G = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 1, fluid, fraction)
    viscosity_L = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 0, fluid, fraction)
    cp_L = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 0, fluid, fraction)
    if Re < 1E6:
        h = 0.9*k_L/D_eq*(rho_L*velocity*D_eq/viscosity_L)**0.5
    else:
        h_shear = 0.9*k_L/D_eq*(rho_L*velocity*D_eq/viscosity_L)**0.5
        if T < 373.15:
            A = 8.071
            B = 1730.63
            C = 233.426
        else:
            A = 8.140
            B = 1810.94
            C =244.485
        T_Antoine = T - 273.15
        p_Antoine = 10**(A - B/(C + T_Antoine))
        T_gi = Thermeprop.getThermProp(mixture, 'T', 'P', p_Antoine, 'Q', 0, fluid, fraction)
        delta_T = T_gi - T_wall
        i_fg = Thermeprop.getThermProp(mixture, 'H','P',p,'Q',1,fluid, fraction) - Thermeprop.getThermProp(mixture, 'H','P',p,'Q',0,fluid, fraction)
        i_fg_prime = i_fg * ( 1 + 0.68*(cp_L*delta_T/i_fg) )
        h_grav = 0.725*( rho_L*(rho_L - rho_G)*9.81*i_fg_prime*k_L**3/(D_eq*viscosity_L*delta_T) )**0.25
        h = (0.5*h_shear**2 + (0.25*h_shear**0.25 + h_grav**4)**0.5)**0.5
        
    return h

# La funzione 'heat_transfer_coefficient_evaporation_Shell' calcola il coefficiente di scambio
# per un fluido bifase che evapora nello shell utilizando la relazione di Thome and Robinson 
def heat_transfer_coefficient_evaporation_Shell(mixture,fraction,p,Q,fluid,pitch,A_cross,OD,q):
    
    print('Calcolo coefficente di scambio termico evaporazione lato shell in corso... \n')
    
    # Calcolo proprietà del fluido
    rho_L = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 0, fluid, fraction)
    rho_G = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 1, fluid, fraction)
    viscosity_L = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 0, fluid, fraction)
    sigma = Thermeprop.getThermProp(mixture, 'surface_tension', 'P', p, 'Q', 1, fluid, fraction)
    cp_L = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 0, fluid, fraction)
    k_L = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 0, fluid, fraction)
    Pr_L = viscosity_L*cp_L/k_L
    pr = p/Thermeprop.getThermProp(mixture, 'Pcrit', fluid, fraction)
    M =  Thermeprop.getThermProp(mixture, 'M', 'P', p, 'Q', 1, fluid, fraction)
    
    # Surface roughness
    R_p = 0.4E-6
    
    # Il calcolo di alcuni parametri richiede un procedimento iterativo dal momento che questi parametri dipendono
    # l'uno dall'altro. Si itera sull'errore percentuale del parametro 'epsilon'. Le iterazioni si fermano quando
    # l'errore è inferiore all'1%
    err = 1
    tollerance = 0.001
    epsilon = 0.1
    while abs(err) > tollerance:
        u_G = Q*A_cross/epsilon/rho_G
        Cap = viscosity_L*u_G/sigma
        Ri = (rho_L - rho_G)**2*9.81*pitch/(A_cross)**2
        S = 1 + 25.7*(Ri*Cap)**0.5*(OD/pitch)
        epsilon_new = ( 1 + (S*(1-Q)*rho_G/Q/rho_L) )**-1
        
        err = (epsilon_new - epsilon)/epsilon
        epsilon = epsilon_new
    
    # Calcolo dei parametri che portano al coefficiente di scambio
    u_L = A_cross*(1 - Q)/rho_L/(1 - epsilon)
    A_hex = 6*(pitch/3)*(pitch/2)
    A_cfl = A_hex - np.pi*OD**2/4
    A_L = A_cfl*(1 - epsilon)
    D_delta = ( 4*A_L/np.pi + OD**2 )**0.5
    delta = (D_delta - OD)/2
    Re_delta = 4*rho_L*u_L*delta/viscosity_L
    alfa_cb = 4.032*Re_delta**0.236*Pr_L**0.4*(k_L/delta)
    alfa_nb = 55*( pr**(0.12 - 0.2*np.log10(R_p))*(-np.log10(pr))**-0.55*M**-0.5*q**0.67 )
    h = (alfa_nb**2 + alfa_cb**2)**0.5
    
    return h

# # La funzione 'heat_transfer_coefficient_evaporation_Shell' calcola il coefficiente di scambio
# # per un fluido bifase che evapora nello shell utilizando la relazione di Honda
# def heat_transfer_coefficient_condensation_Shell(p,Q,fluid,pitch,A_cross,OD,idx,q):
    
    
    


#%% TUBE

def heat_transfer_coefficient_1ph_Tube(mixture,fraction,T,T_out,P,viscosity,Re,fluid,L,ID):
    
    print('Calcolo coefficente di scambio termico mono-fase lato tubi in corso... \n')
    
    # Calcolo del numero di Prandtl
    Cp = Thermeprop.getThermProp(mixture, 'C', 'T', T, 'P', P, fluid, fraction)
    k = Thermeprop.getThermProp(mixture,'L','T', T, 'P', P, fluid, fraction)
    Pr = viscosity*Cp/k
    
    # Laminare
    if Re < 2300:
        Nu = 3.66 + (0.668*Re*Pr*(ID/L))/(1 + 0.40*(Re*Pr*(ID/L))**0.666)
    # Turbolento
    elif Re > 10000:
        m = 0.88 - 0.24/(Pr + 4)
        n = 1/3 + 0.5*np.exp(-0.6*Pr)
        Nu = 5 + 0.015*(Re**m)*(Pr**n)
    # Transizione
    else:
        T_wall = np.mean([T,T_out])
        viscosity_wall = Thermeprop.getThermProp(mixture, 'V', 'T', T_wall, 'Q', 0, fluid, fraction)
        phi = viscosity/viscosity_wall
        Nu = 0.116*(Re**0.666 - 125)*(Pr**0.333)*(phi**0.14)*(1 + (ID/L)**0.666)
        
    h = Nu*k/ID
    
    return h


# La funzione 'heat_exchange_coefficient_evaporation_vertical' calcola il coefficiente di scambio
# per un fluido bifase che evapora in tubi verticali utilizando la relazione di Steiner-Tabore
def heat_exchange_coefficient_evaporation_vertical(mixture,fraction,P,T,Q,Re,viscosity,ID,phi,fluid):
    
    print('Calcolo coefficente di scambio termico evaporazione tubi verticali \n')
    
    # Local normalized nucleate pool boiling coefficient at standard conditions (hNcB,o)
    data = pd.read_excel(r'heat_transfer_input.xlsx')
    data = data.loc[data.Fluid == fluid]
    h_0 = data.h_0
    P_cr = Thermeprop.getThermProp(mixture, 'Pcrit', 'P', P, 'Q', Q, fluid, fraction)
    phi_0 = data.phi_0
    P_n = P/P_cr
    F_pf = 1.2*P_n**0.27 + 2.5*P_n + P_n/(1 - P_n)
    n_f1 = 0.9 - 0.3*P_n**0.3
    R_p = 0.4*1E-6
    R_p0 = 0.4*1E-6
    h_ncb = h_0*F_pf*(phi*phi_0)**n_f1*(R_p*R_p0)**0.133
    
    # Nucleate boiling correction factor (FNcB )
    F_pf = 2.816*P_n**0.45 + (3.4 + 1.7/(1 - P_n**7))*P_n**3.7
    n_f2 = 0.8 - 0.1*np.exp(1.75*P_n)
    M = Thermeprop.getThermProp(mixture, 'M', fluid, fraction)
    F = 0.377 + 0.199*np.log(M) + 0.000028427*M**2
    ID_0 = 0.01
    F_ncb = F_pf**(phi*phi_0)**n_f2*(ID*ID_0)**-0.4*(R_p*R_p0)**0.133*F
    
    # Local liquid-phase forced convection coefficient based on the total flow as liquid (hfo)
    f = (0.7904*np.log(Re) - 1.64)**-2
    cp = Thermeprop.getThermProp(mixture, 'C', 'P', P, 'Q', Q, fluid, fraction)
    k = Thermeprop.getThermProp(mixture, 'L', 'P', P, 'Q', Q, fluid, fraction)
    Pr = viscosity*cp/k
    Nu = f/8*(Re - 1000)*Pr/(1 + 12.7*(f/8)**0.5*(Pr**0.666 - 1))
    h_fo = Nu*k/ID
    
    # Two-phase multiplier (FTP])
    rho_l = Thermeprop.getThermProp(mixture, 'D', 'P', P, 'Q', 0, fluid, fraction)
    rho_g = Thermeprop.getThermProp(mixture, 'D', 'P', P, 'Q', 1, fluid, fraction)
    F_tp = ((1 - Q)**1.5 + 1.9*Q**0.6*(rho_l/rho_g)**0.35)**1.1
    
    h = ((h_ncb*F_ncb)**3 + (h_fo*F_tp)**3)**(1/3)
    
    return h

# La funzione 'heat_exchange_coefficient_evaporation_horizontal' calcola il coefficiente di scambio
# per un fluido bifase che evapora in tubi orizzontali utilizando la relazione di Kattan-Thome-Favrat 
def heat_exchange_coefficient_evaporation_horizontal(mixture,fraction,p,Q,A_pipe,ID,fluid,q,G):
    
    print('Calcolo coefficente di scambio termico evaporazione tubi orizzontali in corso... \n')
    
    # Si assegna il nome al fluido per calcolo Surface Tension con Coolprop. Thermeprop non ancora capace di farlo
    
    if fluid == 10068:
        fluid_coolprop = 'R1233zd(E)'
    elif fluid == 1921:
        fluid_coolprop = 'water'
    else:
        print('Surface tension and molar mass not calculated for the selected fluid')
    
    # Calcolo delle proprietà del flusso.
    rho_L = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 0, fluid, fraction)
    rho_G = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 1, fluid, fraction)
    viscosity_L = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 0, fluid, fraction)
    viscosity_G = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 1, fluid, fraction)
    sigma_L = PropsSI('surface_tension', 'P', p*100000, 'Q', 0, fluid_coolprop) 
    cp_L = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 0, fluid, fraction)
    k_L = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 0, fluid, fraction)
    Pr_L = viscosity_L*cp_L/k_L
    cp_G = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 1, fluid, fraction)
    k_G = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 1, fluid, fraction)
    Pr_G = viscosity_G*cp_G/k_G

    # pr = p/Thermeprop.getThermProp(mixture, 'P', 'C', 0, ' ', 0, fluid, fraction)
    pr = p/(PropsSI('Pcrit', fluid_coolprop)/100000)
    Mol =  PropsSI('M', fluid_coolprop)*1000

    X_tt = ((1 - Q)/Q)**0.875 * (rho_G/rho_L)**0.5*(viscosity_L/viscosity_G)**0.125
    We_g = (G*Q)**2 * ID/rho_G/sigma_L
    Ft = ( (Q**3 * G**2) / (rho_G**2 * 9.81 * ID * (1-Q)) )**0.5
    ## Calcolo dell'alfa: Void Fraction
    
    # Metodo di Lockhart-Martinelli-Weber | Con i parametri m=1.33, n=0.9, o=-0.06, predice bene alfa in microcanali (pochi mm di diametro)
    # per refrigeranti tipo R134a e R410a e miscele aria-acqua con 100 < G < 300 kg/m2.s. Valori ragionevoli con Q > 0.2
    # alfa = ( 1 + ( X_tt + (1/We_g)**1.33 * (rho_L/rho_G)**0.9 ) )**-0.06
    
    # Modello di Rouhani-Axelsson | Inaccurato, cause non ancora capite
    # alfa = (Q/rho_G) * ( (1 + 0.12*(1 - Q)*( Q/rho_G + (1 - Q)/rho_L )) + 1.18*(1 - Q)*( (9.81*sigma_L*(rho_L - rho_G) )**0.25) / (G * rho_L**0.5) )**-1
    
    # Modello di Graham | Migliori predizioni sperimentali per refrigeranti con tubi di ID > 2 mm
    alfa = ( 1 + X_tt + Ft**-1)**-0.321
    
    def theta_stratified(vars, *data):
        A_L, R = data
        theta_strat = vars
        return A_L - (0.5*R**2)*(2*np.pi - theta_strat - np.sin(2*np.pi - theta_strat))
    
    A_L = A_pipe*(1 - alfa)
    data =  (A_L, ID/2)
    theta_strat = fsolve(theta_stratified, np.pi, args = data)[0] 
        
    # G_L = G*A_cross*(1 - Q)
    # Provo a correggere espressione di G_L --> G è già in kg/s. Per come è importato, G non è mass flux, ma mass flow
    G_L = G * A_pipe * (1 - Q)
    
    # Il calcolo di 'theta_dry' richiede la conosconeza dello spessore 'delta' che a sua volta dipende da 
    # 'theta_dry', per questa ragione il calcolo viene effettuato con un procedimento itrativo.
    err = 1
    tollerance = 0.01
    delta = ID/50
    while abs(err) > tollerance:
        A_low = (np.pi - theta_strat/2)*((ID/2)**2 - (ID/2 - delta)**2 )
        A_high = np.pi*((ID/2)**2 - (ID/2 - delta)**2 )
        G_low = G_L/A_low
        G_high = G_L/A_high
        theta_dry = theta_strat*(G_high - G)/(G_high - G_low)
        delta_new = np.pi*ID*(1 - alfa)/2/(2*np.pi - theta_dry)
        
        err = (delta_new - delta)/delta
        delta = delta_new
    
    # Costanti del modello.
    C = 0.0133
    m = 0.69
    
    # Calcolo dei coefficenti di scambio per fase liquida e gassosa e totale.
    h_cb = C*((4*G*(1 - Q)*delta/(1 - alfa)/viscosity_L)**m)*((Pr_L)**0.4)*k_L/delta
    h_nb = 55*( pr**0.12 * (-np.log10(pr))**-0.55 * Mol**-0.5 * q**0.67 )
    h_wet = ( h_nb**3 + h_cb**3)**(1/3)
    h_G = 0.023*(G*Q*ID/alfa/viscosity_G)**0.8*Pr_G**0.4*k_G/ID
    
    h = ( theta_dry*h_G + (2*np.pi - theta_dry)*h_wet ) / (2 * np.pi)
    
    return h


# La funzione 'heat_exchange_coefficient_condensation_vertical' calcola il coefficiente di scambio
# per un fluido bifase che condensa in tubi verticali utilizando la relazione di Carpenter e Colburn 
def heat_exchange_coefficient_condensation_vertical(mixture,fraction,p,Q1,Q2,Re,G,fluid,ID,A_cross):
    
    print('Calcolo coefficente di scambio termico condensazione tubi verticali in corso... \n')
    
    # Calcolo delle proprietà del flusso.
    rho_L = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 0, fluid, fraction)
    rho_G = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 1, fluid, fraction)    
    cp_L = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 0, fluid, fraction)
    k_L = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 0, fluid, fraction)
    viscosity_L = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 0, fluid, fraction)
    Pr_L = viscosity_L*cp_L/k_L
    
    # Surface roughness
    R_p = 0.4E-6
    
    # Calcolo del apparent interfacial friction factor basandosi sul regime di flusso nello scambiatore.
    if Re < 2300:
        f = 64/Re
    elif Re > 2300 and Re < 4000:
        f = 1.2063/Re**0.416
    else:
        def theta_stratified(vars, *data):
            Re, ID, R_p = data
            f = vars
            return 2**-0.5 - 2*np.log(R_p/ID/3.71 + 2.52/Re/f**0.5)
    
        data =  (Re, ID, R_p)
        f = fsolve(theta_stratified, 1.2063/4000**0.416, args = data)[0] 
    
    # Calcolo del flusso di massa di vapore per unità di superficie all'ingresso (1) ed uscita (2) dello scambiatore.
    alfa_1 = (1 + (1 - Q1)/Q1*(rho_G/rho_L)**0.667)**-1
    alfa_2 = (1 + (1 - Q2)/Q2*(rho_G/rho_L)**0.667)**-1
    G1 = Q1*G/(A_cross*alfa_1)
    G2 = Q2*G/(A_cross*alfa_2)
    G_bar = ((G1**2 + G1*G2 + G2**2)/3)**0.5
    
    tau = f*G_bar**2/2/rho_G
    
    h = 0.065*Pr_L**0.5*tau**0.5*rho_L**0.5*k_L*rho_L**0.5/viscosity_L
    
    return h

# La funzione 'heat_exchange_coefficient_condensation_horizontal' calcola il coefficiente di scambio
# per un fluido bifase che condensa in tubi verticali utilizando la relazione di Akers o Shah. 
def heat_exchange_coefficient_condensation_horizontal(mixture,fraction,p,Q,Re,G,fluid,ID,A_cross):

    print('Calcolo coefficente di scambio termico condensazione tubi orizzontali in corso... \n')    

    # Calcolo delle proprietà del fluido
    rho_L = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 0, fluid, fraction)
    rho_G = Thermeprop.getThermProp(mixture, 'D', 'P', p, 'Q', 1, fluid, fraction)
    viscosity_L = Thermeprop.getThermProp(mixture, 'V', 'P', p, 'Q', 0, fluid, fraction)
    cp_L = Thermeprop.getThermProp(mixture, 'C', 'P', p, 'Q', 0, fluid, fraction)
    k_L = Thermeprop.getThermProp(mixture, 'L', 'P', p, 'Q', 0, fluid, fraction)
    Pr_L = viscosity_L*cp_L/k_L
    
    # Se il flusso di massa specifico è inferiore a 200Kg/m^2s viene utilizzata la relazione di
    # Akers per il coefficente di scambio, altrimenti la relazione di Shah. 
    if G/A_cross < 200:
        G_e = G*((1 - Q) + Q*(rho_L/rho_G)**0.5)
        Re_e = ID*G_e/viscosity_L
        if Re > 50000:
            C = 0.0265
            n = 0.8
        else:
            C = 5.03
            n = 1/3
        h = C*Re_e**n*Pr_L**0.333*k_L/ID
    else:
        pr = p/Thermeprop.getThermProp(mixture, 'P', 'C', 0,' ', 0, fluid, fraction)
        Re_L = G * Q * ID / viscosity_L
        h = 0.023*Re_L**0.8*Pr_L**0.4*( (1 - Q)**0.8 + 3.8*Q**0.76*(1 - Q)**0.04/(pr)**0.38 )
        
    return h
    

#%% CORRELAZIONI D3.2 REGEN

# La funzione 'heat_exchange_coefficient_condensation' calcola il coefficiente di scambio per un fluido bifase che condensa   
def heat_exchange_coefficient_condensation(mixture,fraction,T_inlet,T_outlet,P_inlet,P_outlet,x_inlet,x_outlet,D_hydraulic,fluid,idx,M):                                         
    
    x = np.linspace(x_inlet,x_outlet,10)
    alfa_discrete = np.zeros(len(x))
    
    cp_l = Thermeprop.getThermProp(mixture, 'C', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx], fraction)
    k_l = Thermeprop.getThermProp(mixture, 'L', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx], fraction)
    viscosity_l = Thermeprop.getThermProp(mixture, 'V', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx], fraction)
    
    Pr_l = viscosity_l*cp_l/k_l
    rho_l = Thermeprop.getThermProp(mixture, 'D', 'P', np.mean([P_inlet[idx],P_outlet[idx]]), 'Q', 0, fluid[idx], fraction)
    velocity_l = M[idx]/(rho_l*np.pi*0.25*D_hydraulic**2)
    Re_l = rho_l*velocity_l*D_hydraulic/viscosity_l
    
    alfa_l = 0.023*(Re_l**0.4)*(Pr_l**0.4)*k_l/D_hydraulic
    p_star = Thermeprop.getThermProp(mixture, 'D', 'T', np.mean([T_inlet[idx],T_outlet[idx]]), 'Q', 0, fluid[idx], fraction)/Thermeprop.getThermProp(mixture, 'PCRIT', fluid[idx], fraction)
    
    for i in range(len(x)):
        alfa_discrete[i] = alfa_l*((1-x[i])**0.8 + (((1-x[i])**0.04)*3.8*x[i]**0.76)/p_star**0.38)
        
    alfa = np.mean(alfa_discrete)
    
    return alfa

# La funzione 'heat_exchange_evaporation' calcola il coefficiente di scambio e il flusso termico
# per un fluido bifase che evapora.
# N.B.: In questa funzione il significato degli indici è DIVERSO rispetto a quello utilizzato in tutto il resto
# del codice. In questa funzione il primo indice indica le proprietà del fluido che evapora, il secondo indice
# indica le proprietà dell'altro fluido (che in questa caso può raffreddare o condensare).
def heat_exchange_coefficient_evaporation(mixture,fraction,T,p,x,h,rho,D_hydraulic,velocity,fluid,idx,alfa_2,OD,ID,k,h_fouling,LMTD):
    # Il calcolo dello scambio termico avviene tramite un processo iterativo
    tollerance = 0.01
    err = 1
    # Viene assegnato un valor di primo tentativo al 'Boiling number', tramite questo valore calcoliamo i parametri
    # da cui dipende il coefficiente di scambio, calcolato questo valutiamo il flusso termico, aggiorniamo il valore
    # del 'Boiling number' e lo confrontiamo con il valore precedentemente usato. Questa procedura continua fino
    # a che lo scostamento tra due 'Boiling number' consecutivi è inferiore all' 1%.
    Bo = 0.00003
    while err > tollerance:
        viscosity_l = Thermeprop.getThermProp(mixture, 'V', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
        cp_l = Thermeprop.getThermProp(mixture, 'C', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
        k_l = Thermeprop.getThermProp(mixture, 'L', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
        Pr_l = viscosity_l*cp_l/k_l
        rho_l = Thermeprop.getThermProp(mixture, 'D', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
        rho_g = Thermeprop.getThermProp(mixture, 'D', 'P', p[idx], 'Q', 1, fluid[idx], fraction)
        alfa_l = 0.023*((rho[idx]*velocity[idx]*(1-x[idx])*D_hydraulic[idx]/viscosity_l)**0.8)*(Pr_l**0.4)*k_l/D_hydraulic[idx]
        C_0 = ((1/x[idx]-1)**0.8)*(rho_l/rho_g)**0.5
        Fr_l = ((rho[idx]*velocity[idx])**2)/((rho_l**2)*9.8*D_hydraulic[idx])
        heat_flux = Bo*velocity[idx]*rho[idx]*(Thermeprop.getThermProp(mixture, 'H', 'P', p[idx], 'Q', 1, fluid[idx], fraction) - Thermeprop.getThermProp(mixture, 'H', 'P', p[idx], 'Q', 0, fluid[idx], fraction))
        #
        if Bo > 0.0011:
            F = 14.70
        else:
            F = 15.43
        #
        if Fr_l >= 0.04:
            N = C_0
        else:
            N = 0.38*C_0*Fr_l**-0.3
        #
        psi_cb = 1.8*N**-0.8
        #
        if N > 0.1 and Bo >= 0.00003:
            psi_nb = 230*Bo**0.5
        elif N > 0.1 and Bo < 0.00003:
            psi_nb = 1 + 46*Bo**0.5
        #
        if N >= 0.1 and N <= 1:
            psi_bs = np.exp(2.74*N**-0.15)*F*Bo**0.5
        elif N < 0.1:
            psi_bs = np.exp(2.74*N**-0.1)*F*Bo**0.5
        #
        if N <= 0.1:
            psi = np.max([psi_bs,psi_cb])
        else:
            psi = np.max([psi_nb,psi_cb])
        #
        alfa_ev = psi*alfa_l
        # Noti gli 'alfa' possiamo calcolare il coefficiente di scambio globale
        if idx == 0:    
            U = (1/alfa_2 + 1/h_fouling[1] + OD*np.log(OD/ID)/(2*k) + (1/alfa_ev)*(OD/ID) + (1/h_fouling[0])*(OD/ID))**-1
        else:
            U = (1/alfa_ev + 1/h_fouling[1] + OD*np.log(OD/ID)/(2*k) + (1/alfa_2)*(OD/ID) + (1/h_fouling[0])*(OD/ID))**-1
        # Utilizziamo il coefficiente di scambio globale per valutare il flusso termico specifico
        heat_flux = U*LMTD
        # Aggiorniamo il valore del 'Boiling number'
        Bo_next = abs(heat_flux/(velocity[0]*rho[0]*(Thermeprop.getThermProp(mixture, 'H', 'P', p[0], 'Q', 1, fluid[0], fraction) - Thermeprop.getThermProp(mixture, 'H', 'P', p[0], 'Q', 0, fluid[0], fraction))))
        # Valutiamo lo scostamento
        err = abs((Bo_next-Bo)/Bo)
        Bo = Bo_next
    
    # L'ultimo 'alfa_ev' trovato corrisponde al valore cercato del coefficiente di scambio
    alfa = alfa_ev
    
    return alfa

# Metodo iterativo per il calcolo del fattore di scambio in evaporazione
def heat_exchange_coefficient_evaporation_Liu(mixture,fraction,T,p,x,h,rho,D_hydraulic,velocity,fluid,idx,alfa_2,OD,ID,k,h_fouling,LMTD,q_0,M,data,Ft):
     
    err = 1
    tollerance = 0.01
    
    viscosity_l = Thermeprop.getThermProp(mixture, 'V', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
    cp_l = Thermeprop.getThermProp(mixture, 'C', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
    k_l = Thermeprop.getThermProp(mixture, 'L', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
    Pr_l = viscosity_l*cp_l/k_l
    rho_l = Thermeprop.getThermProp(mixture, 'D', 'P', p[idx], 'Q', 0, fluid[idx], fraction)
    rho_g = Thermeprop.getThermProp(mixture, 'D', 'P', p[idx], 'Q', 1, fluid[idx], fraction)
    velocity_l = M[idx]/(rho_l*np.pi*0.25*D_hydraulic[idx]**2)
    Re_l = rho_l * velocity_l * D_hydraulic[idx] / viscosity_l
    F = (1 + x[idx]*Pr_l*(rho_l/rho_g - 1))**0.35
    S = (1 + 0.55*F**0.1 * Re_l**0.16)**-1
    alfa_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D_hydraulic[idx]
    
    viscosity = Thermeprop.getThermProp(mixture, 'V', 'P', p[idx], 'H', h[idx], fluid[idx], fraction)
    cp = Thermeprop.getThermProp(mixture, 'C', 'P', p[idx], 'H', h[idx], fluid[idx], fraction)
    k = Thermeprop.getThermProp(mixture, 'L', 'P', p[idx], 'H', h[idx], fluid[idx], fraction)
    Pr = viscosity*cp/k
    
    q = q_0
    while err > tollerance:
        M_mol = Thermeprop.getThermProp(mixture, 'D', 'P', p[idx], 'H', h[idx], fluid[idx], fraction)/1000/Thermeprop.getThermProp(mixture, 'DMOLAR', 'P', p[idx], 'H', h[idx], fluid[idx], fraction)
        alfa_nb = 55*Pr**(0.12 - 0.434*np.log(data.Selected_Value[31])) * -np.log10(Pr)**-0.55 * M_mol**-0.5 * q**0.67
        alfa_ev = ((F*alfa_l)**2 + (S*alfa_nb)**2)**0.5
        if idx == 0:    
            U = (1/alfa_2 + 1/h_fouling[1] + OD*np.log(OD/ID)/(2*k) + (1/alfa_ev)*(OD/ID) + (1/h_fouling[0])*(OD/ID))**-1
        else:
            U = (1/alfa_ev + 1/h_fouling[1] + OD*np.log(OD/ID)/(2*k) + (1/alfa_2)*(OD/ID) + (1/h_fouling[0])*(OD/ID))**-1
        
        q = U * LMTD * Ft
        
        err = abs((q-q_0)/q_0)
        
        q_0 = q
        
    return alfa_ev
        
        
    
    
    
    
    
    
    