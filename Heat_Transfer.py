# -*- coding: utf-8 -*-
"""
Created on Fri May 21 14:45:45 2021

@author: Paolo
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

# Definizione del coefficiente di sporcamento. Attualmente è imposto costante.
h_fouling = [1E6,1E6]

# La funzione 'heat_flux_wrapper' prende in input tutto il contenuto di un volume di controllo dello scambiatore
# e fornisce in output il valore dello scambio termico specifico e del coefficiente di scambio.
# Lo scopo della funzione è quello di stabilire quale meccanismo di scambio termico va richiamato per il calcolo
# e di fornire alle funzioni di calcolo i valori in input corretti per effettuare la modellazione.
def heat_flux_wrapper(HE_CV):
    # Identifichiamo all'interno della funzione tutti i possibili casi che possono verificarsi all'interno dello
    # scambiatore in termini di meccanismi di scambio termico e posizione dei fluidi. 
    
    # Shell side evaporation
    if HE_CV.temperature_in[0] > HE_CV.temperature_in[1]:
        if HE_CV.quality_in[1] != -1.0:  # HE_CV.quality_in[1] >= 0 and HE_CV.quality_in[1] <= 1:
            # Definizione degli input da fonire alla funzione 'heat_exchange_evaporation'
            x = [HE_CV.quality_in[1],HE_CV.quality_in[0]]
            p = [HE_CV.pressure_in[1],HE_CV.pressure_in[0]]
            T = [HE_CV.temperature_in[1],HE_CV.temperature_in[0]]
            h = [HE_CV.enthalpy_in[1],HE_CV.enthalpy_in[0]]
            rho = [HE_CV.density_in[1],HE_CV.density_in[0]]
            D_hydr = [HE_CV.D_hydraulic[1],HE_CV.D_hydraulic[0]]
            velocity = [HE_CV.velocity[1],HE_CV.velocity[0]]
            fluid = [HE_CV.fluid[1],HE_CV.fluid[0]]
            [heat_flux,alfa_ev,alfa_2] = heat_exchange_evaporation(HE_CV,T,p,x,h,rho,D_hydr,velocity,fluid)
            alfa = [alfa_2, alfa_ev]
            
    # Pipe side evaporation
    if HE_CV.temperature_in[0] < HE_CV.temperature_in[1]:
        if HE_CV.quality_in[0] != -1.0: # HE_CV.quality_in[0] >= 0 and HE_CV.quality_in[0] <= 1: 
            # Definizione degli input da fonire alla funzione 'heat_exchange_evaporation'
            x = HE_CV.quality_in
            p = HE_CV.pressure_in
            T = HE_CV.temperature_in
            h = HE_CV.enthalpy_in
            rho = HE_CV.density_in
            D_hydr = HE_CV.D_hydraulic
            velocity = HE_CV.velocity
            fluid = HE_CV.fluid
            [heat_flux,alfa_ev,alfa_2] = heat_exchange_evaporation(HE_CV,T,p,x,h,rho,D_hydr,velocity,fluid)
            alfa = [alfa_ev,alfa_2]
    
    # Pipe side hot - Shell side cold
    if HE_CV.temperature_in[0] > HE_CV.temperature_in[1]:
        #
        # 1-phase heat transfer on both side
        if HE_CV.quality_in[0] == -1.0 and HE_CV.quality_in[1] == -1.0: # (HE_CV.quality_in[1] < 0 and HE_CV.quality_in[1] > 1) and (HE_CV.quality_in[0] < 0 and HE_CV.quality_in[0] > 1):
            # Definizione degli input da fonire alla funzione 'heat_exchange_coefficient_1phase'
            x = HE_CV.quality_in
            p = HE_CV.pressure_in
            T = HE_CV.temperature_in
            h = HE_CV.enthalpy_in
            rho = HE_CV.density_in
            velocity = HE_CV.velocity
            D = HE_CV.D_hydraulic
            fluid = HE_CV.fluid
            alfa = [0,0]
            alfa[0] = heat_exchange_coefficient_1phase(T[0],p[0],x[0],h[0],rho[0],D[0],velocity[0],fluid[0])
            alfa[1] = heat_exchange_coefficient_1phase(T[1],p[1],x[1],h[1],rho[1],D[1],velocity[1],fluid[1])
            # Noti gli 'alfa' possiamo calcolare il coefficiente di scambio globale
            U = (1/alfa[1] + 1/h_fouling[1] + HE_CV.D_pipe[1]*np.log(HE_CV.D_pipe[1]/HE_CV.D_pipe[0])/(2*HE_CV.k) + (1/alfa[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]) + (1/h_fouling[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]))**-1
            # Utilizziamo il coefficiente di scambio globale per valutare il flusso termico specifico
            heat_flux = U*(T[0] - T[1])
        #
        # Condensation in pipe side 
        elif HE_CV.quality_in[0] != -1.0 and HE_CV.quality_in[1] == -1.0:
            # Definizione degli input da fonire alle funzioni 'heat_exchange_coefficient_1phase' e 'heat_exchange_coefficient_condensation'
            x = HE_CV.quality_in
            p = HE_CV.pressure_in
            T = HE_CV.temperature_in
            h = HE_CV.enthalpy_in
            rho = HE_CV.density_in
            velocity = HE_CV.velocity
            D = HE_CV.D_hydraulic
            fluid = HE_CV.fluid
            alfa = [0,0]
            alfa[0] = heat_exchange_coefficient_condensation(T[0],p[0],x[0],D[0],velocity[0],fluid[0])
            alfa[1] = heat_exchange_coefficient_1phase(T[1],p[1],x[1],h[1],rho[1],D[1],velocity[1],fluid[1])
            # Noti gli 'alfa' possiamo calcolare il coefficiente di scambio globale
            U = (1/alfa[1] + 1/h_fouling[1] + HE_CV.D_pipe[1]*np.log(HE_CV.D_pipe[1]/HE_CV.D_pipe[0])/(2*HE_CV.k) + (1/alfa[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]) + (1/h_fouling[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]))**-1
            # Utilizziamo il coefficiente di scambio globale per valutare il flusso termico specifico
            heat_flux = U*(T[0] - T[1])
            
    # Pipe side cold - Shell side hot
    if HE_CV.temperature_in[1] > HE_CV.temperature_in[0]:
        # 1-phase heat transfer on both side
        if HE_CV.quality_in[0] == -1.0 and HE_CV.quality_in[1] == -1.0: # (HE_CV.quality_in[1] < 0 and HE_CV.quality_in[1] > 1) and (HE_CV.quality_in[0] < 0 and HE_CV.quality_in[0] > 1):
            # Definizione degli input da fonire alla funzione 'heat_exchange_coefficient_1phase'
            x = HE_CV.quality_in
            p = HE_CV.pressure_in
            T = HE_CV.temperature_in
            h = HE_CV.enthalpy_in
            rho = HE_CV.density_in
            velocity = HE_CV.velocity
            D = HE_CV.D_hydraulic
            fluid = HE_CV.fluid
            alfa = [0,0]
            alfa[0] = heat_exchange_coefficient_1phase(T[0],p[0],x[0],h[0],rho[0],D[0],velocity[0],fluid[0])
            alfa[1] = heat_exchange_coefficient_1phase(T[1],p[1],x[1],h[1],rho[1],D[1],velocity[1],fluid[1])
            # Noti gli 'alfa' possiamo calcolare il coefficiente di scambio globale
            U = (1/alfa[1] + 1/h_fouling[1] + HE_CV.D_pipe[1]*np.log(HE_CV.D_pipe[1]/HE_CV.D_pipe[0])/(2*HE_CV.k) + (1/alfa[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]) + (1/h_fouling[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]))**-1
            # Utilizziamo il coefficiente di scambio globale per valutare il flusso termico specifico
            heat_flux = U*(T[0] - T[1])
            
        # Condensation in shell side
        elif HE_CV.quality_in[0] == -1.0 and HE_CV.quality_in[1] != -1.0: # (HE_CV.quality_in[0] < 0 and HE_CV.quality_in[0] > 1) and (HE_CV.quality_in[1] >= 0 and HE_CV.quality_in[1] <= 1):
            # Definizione degli input da fonire alle funzioni 'heat_exchange_coefficient_1phase' e 'heat_exchange_coefficient_condensation'
            x = HE_CV.quality_in
            p = HE_CV.pressure_in
            T = HE_CV.temperature_in
            h = HE_CV.enthalpy_in
            rho = HE_CV.density_in
            velocity = HE_CV.velocity
            D = HE_CV.D_hydraulic
            fluid = HE_CV.fluid
            alfa = [0,0]
            alfa[0] = heat_exchange_coefficient_1phase(T[0],p[0],x[0],h[0],rho[0],D[0],velocity[0],fluid[0])
            alfa[1] = heat_exchange_coefficient_condensation(T[1],p[1],x[1],D[1],velocity[1],fluid[1])
            # Noti gli 'alfa' possiamo calcolare il coefficiente di scambio globale
            U = (1/alfa[1] + 1/h_fouling[1] + HE_CV.D_pipe[1]*np.log(HE_CV.D_pipe[1]/HE_CV.D_pipe[0])/(2*HE_CV.k) + (1/alfa[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]) + (1/h_fouling[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]))**-1
            # Utilizziamo il coefficiente di scambio globale per valutare il flusso termico specifico
            heat_flux = U*(T[0] - T[1])
            
    return heat_flux, alfa

# La funzione 'heat_exchange_coefficient_1phase' calcola il coefficiente di scambio per un fluido monofase
def heat_exchange_coefficient_1phase(T,p,x,h,rho,D_hydraulic,velocity,fluid):
    k = PropsSI('L', 'P', p, 'H', h, fluid)
    Cp = PropsSI('Cpmass', 'P', p, 'H', h, fluid)
    viscosity = PropsSI('V', 'P', p, 'H', h, fluid)
    Re = rho*velocity*D_hydraulic/viscosity
    Pr = viscosity*Cp/k
    # Condizione di flusso turbolento
    if Re > 10000:
        alfa = 0.023*(Re**0.8)*(Pr**0.33)*k/D_hydraulic
    # Condizione di flusso laminare
    elif Re < 2100:
        Nu = 3.66
        alfa = Nu*k/D_hydraulic
    # Nella condizione di transitorio calcoliamo 'alfa' linearizzando tra i valori limiti per flusso
    # turbolento (i.e. Re = 10000) e flusso laminare (i.e. Re = 2100)
    else:
        alfa_tur = 0.023*(10000**0.8)*(Pr**0.33)*k/D_hydraulic
        alfa_lam = 3.66*k/D_hydraulic
        alfa = alfa_lam + (alfa_tur - alfa_lam)/(10000 - 2100)*(Re - 2100)
    
    return alfa

# La funzione 'heat_exchange_coefficient_condensation' calcola il coefficiente di scambio per un fluido bifase che condensa   
def heat_exchange_coefficient_condensation(T,p,x,D_hydraulic,velocity,fluid):                                         
    cp_l = PropsSI('Cpmass', 'P', p, 'Q', 0, fluid)
    k_l = PropsSI('L', 'P', p, 'Q', 0, fluid)
    viscosity_l = PropsSI('V', 'P', p, 'Q', 0, fluid)
    Pr_l = viscosity_l*cp_l/k_l
    rho_l = PropsSI('D', 'P', p, 'Q', 0, fluid)
    Re_l = rho_l*velocity*D_hydraulic/viscosity_l
    alfa_l = 0.023*(Re_l**0.4)*(Pr_l**0.4)*k_l/D_hydraulic
    p_star = PropsSI('D', 'T', T, 'Q', 0, fluid)/PropsSI('PCRIT', fluid)
    alfa = alfa_l*((1-x)**0.8 + (((1-x)**0.04)*3.8*x**0.76)/p_star**0.38)
    
    return alfa
 
# La funzione 'heat_exchange_evaporation' calcola il coefficiente di scambio e il flusso termico
# per un fluido bifase che evapora.
# N.B.: In questa funzione il significato degli indici è DIVERSO rispetto a quello utilizzato in tutto il resto
# del codice. In questa funzione il primo indice indica le proprietà del fluido che evapora, il secondo indice
# indica le proprietà dell'altro fluido (che in questa caso può raffreddare o condensare).
def heat_exchange_evaporation(HE_CV,T,p,x,h,rho,D_hydraulic,velocity,fluid):
    # Il calcolo dello scambio termico avviene tramite un processo iterativo
    tollerance = 0.01
    err = 1
    # Viene assegnato un valor di primo tentativo al 'Boiling number', tramite questo valore calcoliamo i parametri
    # da cui dipende il coefficiente di scambio, calcolato questo valutiamo il flusso termico, aggiorniamo il valore
    # del 'Boiling number' e lo confrontiamo con il valore precedentemente usato. Questa procedura continua fino
    # a che lo scostamento tra due 'Boiling number' consecutivi è inferiore all' 1%.
    # Il valore di primo tentativo per il 'Boiling number' dovrà essere assegnato sulla base del Boiling number 
    # al volume di controllo precedente (****DA IMPLEMENTARE****)
    Bo = 0.0011
    while err > tollerance:
        viscosity_l = PropsSI('V', 'P', p[0], 'Q', 0, fluid[0])
        cp_l = PropsSI('Cpmass', 'P', p[0], 'Q', 0, fluid[0])
        k_l = PropsSI('L', 'P', p[0], 'Q', 0, fluid[0])
        Pr_l = viscosity_l*cp_l/k_l
        rho_l = PropsSI('D', 'P', p[0], 'Q', 0, fluid[0])
        rho_g = PropsSI('D', 'P', p[0], 'Q', 1, fluid[0])
        alfa_l = 0.023*((rho[0]*velocity[0]*(1-x[0])*D_hydraulic[0]/viscosity_l)**0.8)*(Pr_l**0.4)*k_l/D_hydraulic[0]
        C_0 = ((1/x[0]-1)**0.8)*(rho_l/rho_g)**0.5
        Fr_l = ((rho[0]*velocity[0])**2)/((rho_l**2)*9.8*D_hydraulic[0])
        heat_flux = Bo*velocity[0]*rho[0]*(PropsSI('H', 'P', p[0], 'Q', 1, fluid[0]) - PropsSI('H', 'P', p[0], 'Q', 0, fluid[0]))
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
        if N > 0.1 and Bo > 0.00003:
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
        # Calcolo del coefficiente di scambio del fluido dall'altra parte della parete. Viene chiamata la funzione
        # 'heat_exchange_coefficient_1phase' se il fluido è monofase o 'heat_exchange_coefficient_condensation' se
        # il fluido è bifase.
        if PropsSI('Q', 'P', p[1], 'H', h[1], fluid[1]) == -1.0:
            alfa_2 = heat_exchange_coefficient_1phase(T[1],p[1],x[1],h[1],rho[1],D_hydraulic[1],velocity[1],fluid[1])
        else:
            alfa_2 = heat_exchange_coefficient_condensation(T[1],p[1],x[1],D_hydraulic[1],velocity[1],fluid[1])
        # Noti gli 'alfa' possiamo calcolare il coefficiente di scambio globale
        U = (1/alfa_2 + 1/h_fouling[1] + HE_CV.D_pipe[1]*np.log(HE_CV.D_pipe[1]/HE_CV.D_pipe[0])/(2*HE_CV.k) + (1/alfa_ev)*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]) + (1/h_fouling[0])*(HE_CV.D_pipe[1]/HE_CV.D_pipe[0]))**-1
        # Utilizziamo il coefficiente di scambio globale per valutare il flusso termico specifico
        heat_flux = U*(T[0] - T[1])
        # Aggiorniamo il valore del 'Boiling number'
        Bo_next = abs(heat_flux/(velocity[0]*rho[0]*(PropsSI('H', 'P', p[0], 'Q', 1, fluid[0]) - PropsSI('H', 'P', p[0], 'Q', 0, fluid[0]))))
        # Valutiamo lo scostamento
        err = abs((Bo_next-Bo)/Bo)
        Bo = Bo_next
        
    return heat_flux, alfa_ev, alfa_2
