# -*- coding: utf-8 -*-

"""
Created on Wed Jun 23 12:18:51 2021

@author: Paolo
"""

import numpy as np
import Thermeprop

def InletOutlet(mixture, fraction, fluid, P_loss, P_inlet = 'None', T_inlet = 'None', Q_inlet = 'None', T_outlet = 'None', Q_outlet = 'None'):
    
    # INLET
    if Q_inlet == 'None':
        # 
        # Pressure and temperature are given
        pressure_inlet = P_inlet # [Pa]
        temperature_inlet = T_inlet # [K]
        quality_inlet = ''
        # check su input
        temperature_sat = Thermeprop.getThermProp(mixture, 'T', 'P', pressure_inlet, 'Q', 1, fluid, fraction)
        t_sat_sup = temperature_sat + 1
        t_sat_inf = temperature_sat - 1
        if temperature_inlet > t_sat_inf and temperature_inlet < t_sat_sup:
            print('input errato: sono stati inseriti input di T e P, ma risulta essere sotto campana')
        # Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
        # se il flusso è bifase utilizziamo T e titolo di vapore.
        # print (temperature_inlet, pressure_inlet)
        density_inlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction) 
        enthalpy_inlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction)
        viscosity_inlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction)
        #
    elif P_inlet == 'None':
        #
        # Temperature and quality are given
        quality_inlet = Q_inlet # [Pa]
        temperature_inlet = T_inlet # [K]
        pressure_inlet = Thermeprop.getThermProp(mixture, 'P', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction)
        # Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
        # se il flusso è bifase utilizziamo T e titolo di vapore.
        if quality_inlet == -1.0:
            density_inlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction) 
            enthalpy_inlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction)
            viscosity_inlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction)
        else:  
            density_inlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction) 
            enthalpy_inlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction)
            viscosity_inlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction)
    #
    else:
        #
        # Pressure and quality are given
        quality_inlet = Q_inlet # [Pa]
        pressure_inlet = P_inlet # [K]
        temperature_inlet = Thermeprop.getThermProp(mixture, 'T', 'P', pressure_inlet, 'Q', quality_inlet, fluid, fraction)
        # Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
        # se il flusso è bifase utilizziamo T e titolo di vapore.
        if quality_inlet == -1.0:
            density_inlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction) 
            enthalpy_inlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction)
            viscosity_inlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_inlet, 'P', pressure_inlet, fluid, fraction)
        else:  
            density_inlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction) 
            enthalpy_inlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction)
            viscosity_inlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_inlet, 'Q', quality_inlet, fluid, fraction)
    
        
    # OUTLET
    #
    if Q_outlet == 'None':
        #
        # Per il calcolo delle consizioni di uscita, si presuppone che il FLUIDO SECONDARIO subisca perdite di carico pari a quelle massime accettabili
        pressure_outlet = pressure_inlet - P_loss   # CAPIRE COME DEFINIRE LO STATO DEL FLUSSO
        temperature_outlet = T_outlet # [K]
        quality_outlet = ''
        # check su input
        temperature_sat = Thermeprop.getThermProp(mixture, 'T', 'P', pressure_outlet, 'Q', 1, fluid, fraction)
        t_sat_sup = temperature_sat + 1
        t_sat_inf = temperature_sat - 1
        if temperature_outlet > t_sat_inf and temperature_outlet < t_sat_sup:
            print('input errato: sono stati inseriti input di T e P, ma risulta essere sotto campana')
        density_outlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_outlet, 'P', pressure_outlet, fluid, fraction) 
        enthalpy_outlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_outlet, 'P', pressure_outlet, fluid, fraction)
        viscosity_outlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_outlet, 'P', pressure_outlet, fluid, fraction)
    
    elif T_outlet == 'None':
        #
        # Per il calcolo delle consizioni di uscita, si presuppone che il FLUIDO SECONDARIO subisca perdite di carico pari a quelle massime accettabili
        pressure_outlet = pressure_inlet - P_loss   # CAPIRE COME DEFINIRE LO STATO DEL FLUSSO
        quality_outlet = Q_outlet
        temperature_outlet = Thermeprop.getThermProp(mixture, 'T', 'P', pressure_outlet, 'Q', quality_outlet, fluid, fraction)
        if quality_outlet == -1.0:
            density_outlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_outlet, 'P', pressure_outlet, fluid, fraction) 
            enthalpy_outlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_outlet, 'P', pressure_outlet, fluid, fraction)
            viscosity_outlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_outlet, 'P', pressure_outlet, fluid, fraction)
        else:
            density_outlet = Thermeprop.getThermProp(mixture, 'D', 'T', temperature_outlet, 'Q', quality_outlet, fluid, fraction) 
            enthalpy_outlet = Thermeprop.getThermProp(mixture, 'H', 'T', temperature_outlet, 'Q', quality_outlet, fluid, fraction)
            viscosity_outlet = Thermeprop.getThermProp(mixture, 'V', 'T', temperature_outlet, 'Q', quality_outlet, fluid, fraction)


    return pressure_inlet, temperature_inlet, quality_inlet, density_inlet, enthalpy_inlet, viscosity_inlet, pressure_outlet, temperature_outlet, quality_outlet, density_outlet, enthalpy_outlet, viscosity_outlet

# La temperatura di uscita dell'aria dall'air cooler è valutata sulla base del valore di primo tentativo di Uo
def air_cooler_outlet(Uo,T_inlet,T_outlet):
    # Le unità di misure sono convertite da SI a sistema imperiale
    T_inlet_wf = (T_inlet[0] - 273.15)*9/5 + 32
    T_inlet_sf = (T_inlet[1] - 273.15)*9/5 + 32
    T_outlet_wf = (T_outlet[0] - 273.15)*9/5 + 32
    
    Uo = 0.1761*Uo
    
    T_air_out = T_inlet_sf + (Uo + 1)*(np.mean([T_inlet_wf,T_outlet_wf]) - T_inlet_sf)/10
    T_air_out = (T_air_out - 32)*5/9 + 273.15
    
    return T_air_out