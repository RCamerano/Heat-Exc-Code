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
    
# La funzione 'heat_exchange_coefficient_condensation' calcola il coefficiente di scambio per fluido che riscalda/raffredda     
def heat_transfer_coefficient_1ph(T,p,fluid,Re,j_h,D_hydraulic,idx,L):
    
    if Re[idx] > 10000:
        # Dittus-Boelter correlation - Tube side
        viscosity_mean = PropsSI('V', 'T', T[idx], 'Q', 0, fluid[idx])
        T_wall = np.mean([T[0],T[1]])
        viscosity_wall = PropsSI('V', 'T', T_wall, 'Q', 0, fluid[idx])
        
        k = PropsSI('L', 'P', p[idx], 'T', T[idx], fluid[idx])
        C_p = PropsSI('Cpmass', 'P', p[idx], 'T', T[idx], fluid[idx])
        
        Pr = viscosity_mean*C_p/k
        
        Nu = j_h[idx]*(Re[idx]**0.8)*(Pr**0.33)*(viscosity_mean/viscosity_wall)**0.14
            
        alfa = Nu*k/D_hydraulic
        
    elif Re[idx] < 2100:
        Nu = 3.66
        alfa = Nu*k/D_hydraulic
    
    else:
        viscosity_mean = PropsSI('V', 'T', T[idx], 'Q', 0, fluid[idx])
        T_wall = np.mean([T[0],T[1]])
        viscosity_wall = PropsSI('V', 'T', T_wall, 'Q', 0, fluid[idx])
        
        k = PropsSI('L', 'P', p[idx], 'T', T[idx], fluid[idx])
        C_p = PropsSI('Cpmass', 'P', p[idx], 'T', T[idx], fluid[idx])
        
        Pr = viscosity_mean*C_p/k
        
        Nu_10000 = j_h[idx]*(10000**0.8)*(Pr**0.33)*(viscosity_mean/viscosity_wall)**0.14
        
        Nu_2100 = 3.66
        
        m = (Nu_10000 - Nu_2100) / (10000 - 2100)
        q = Nu_2100 - m*2100
        Nu = m*Re[idx] + q
        
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

# La funzione 'heat_exchange_evaporation' calcola il coefficiente di scambio e il flusso termico
# per un fluido bifase che evapora.
# N.B.: In questa funzione il significato degli indici è DIVERSO rispetto a quello utilizzato in tutto il resto
# del codice. In questa funzione il primo indice indica le proprietà del fluido che evapora, il secondo indice
# indica le proprietà dell'altro fluido (che in questa caso può raffreddare o condensare).
def heat_exchange_coefficient_evaporation(T,p,x,h,rho,D_hydraulic,velocity,fluid,idx,alfa_2,OD,ID,k,h_fouling,LMTD):
    # Il calcolo dello scambio termico avviene tramite un processo iterativo
    tollerance = 0.01
    err = 1
    # Viene assegnato un valor di primo tentativo al 'Boiling number', tramite questo valore calcoliamo i parametri
    # da cui dipende il coefficiente di scambio, calcolato questo valutiamo il flusso termico, aggiorniamo il valore
    # del 'Boiling number' e lo confrontiamo con il valore precedentemente usato. Questa procedura continua fino
    # a che lo scostamento tra due 'Boiling number' consecutivi è inferiore all' 1%.
    Bo = 0.00003
    while err > tollerance:
        viscosity_l = PropsSI('V', 'P', p[idx], 'Q', 0, fluid[idx])
        cp_l = PropsSI('Cpmass', 'P', p[idx], 'Q', 0, fluid[idx])
        k_l = PropsSI('L', 'P', p[idx], 'Q', 0, fluid[idx])
        Pr_l = viscosity_l*cp_l/k_l
        rho_l = PropsSI('D', 'P', p[idx], 'Q', 0, fluid[idx])
        rho_g = PropsSI('D', 'P', p[idx], 'Q', 1, fluid[idx])
        alfa_l = 0.023*((rho[idx]*velocity[idx]*(1-x[idx])*D_hydraulic[idx]/viscosity_l)**0.8)*(Pr_l**0.4)*k_l/D_hydraulic[idx]
        C_0 = ((1/x[idx]-1)**0.8)*(rho_l/rho_g)**0.5
        Fr_l = ((rho[idx]*velocity[idx])**2)/((rho_l**2)*9.8*D_hydraulic[idx])
        heat_flux = Bo*velocity[idx]*rho[idx]*(PropsSI('H', 'P', p[idx], 'Q', 1, fluid[idx]) - PropsSI('H', 'P', p[idx], 'Q', 0, fluid[idx]))
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
        Bo_next = abs(heat_flux/(velocity[0]*rho[0]*(PropsSI('H', 'P', p[0], 'Q', 1, fluid[0]) - PropsSI('H', 'P', p[0], 'Q', 0, fluid[0]))))
        # Valutiamo lo scostamento
        err = abs((Bo_next-Bo)/Bo)
        Bo = Bo_next
    
    # L'ultimo 'alfa_ev' trovato corrisponde al valore cercato del coefficiente di scambio
    alfa = alfa_ev
        
    return alfa

def heat_exchange_coefficient_evaporation_Liu(T,p,x,h,rho,D_hydraulic,velocity,fluid,idx,alfa_2,OD,ID,k,h_fouling,LMTD,q_0,M,data,Ft):
     
    err = 1
    tollerance = 0.01
    
    viscosity_l = PropsSI('V', 'P', p[idx], 'Q', 0, fluid[idx])
    cp_l = PropsSI('Cpmass', 'P', p[idx], 'Q', 0, fluid[idx])
    k_l = PropsSI('L', 'P', p[idx], 'Q', 0, fluid[idx])
    Pr_l = viscosity_l*cp_l/k_l
    rho_l = PropsSI('D', 'P', p[idx], 'Q', 0, fluid[idx])
    rho_g = PropsSI('D', 'P', p[idx], 'Q', 1, fluid[idx])
    velocity_l = M[idx]/(rho_l*np.pi*0.25*D_hydraulic[idx]**2)
    Re_l = rho_l * velocity_l * D_hydraulic[idx] / viscosity_l
    F = (1 + x[idx]*Pr_l*(rho_l/rho_g - 1))**0.35
    S = (1 + 0.55*F**0.1 * Re_l**0.16)**-1
    alfa_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D_hydraulic[idx]
    
    viscosity = PropsSI('V', 'P', p[idx], 'H', h[idx], fluid[idx])
    cp = PropsSI('Cpmass', 'P', p[idx], 'H', h[idx], fluid[idx])
    k = PropsSI('L', 'P', p[idx], 'H', h[idx], fluid[idx])
    Pr = viscosity*cp/k
    
    q = q_0
    while err > tollerance:
        M_mol = PropsSI('D', 'P', p[idx], 'H', h[idx], fluid[idx])/1000/PropsSI('DMOLAR', 'P', p[idx], 'H', h[idx], fluid[idx])
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
        
        
    
    
    
    
    
    
    