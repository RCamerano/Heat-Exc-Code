# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 12:18:51 2021

@author: Paolo
"""

from CoolProp.CoolProp import PropsSI

def InletOutlet(fluid,P_loss, P_inlet = None, T_inlet = None, Q_inlet = None, T_outlet = None, Q_outlet = None):
    
    # INLET
    if Q_inlet == None:
        # 
        # Pressure and temperature are given
        pressure_inlet = P_inlet #[Pa]
        temperature_inlet = T_inlet #[K]
        quality_inlet = PropsSI('Q', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
        # Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
        # se il flusso è bifase utilizziamo T e titolo di vapore.
        if quality_inlet == -1.0:
            density_inlet = PropsSI('D', 'T', temperature_inlet, 'P', pressure_inlet, fluid) 
            enthalpy_inlet = PropsSI('H', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
            viscosity_inlet = PropsSI('V', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
        else:  
            density_inlet = PropsSI('D', 'T', temperature_inlet, 'Q', quality_inlet, fluid) 
            enthalpy_inlet = PropsSI('H', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
            viscosity_inlet = PropsSI('V', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
        #
    elif P_inlet == None:
        #
        # Temperature and quality are given
        quality_inlet = Q_inlet #[Pa]
        temperature_inlet = T_inlet #[K]
        pressure_inlet = PropsSI('P', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
        # Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
        # se il flusso è bifase utilizziamo T e titolo di vapore.
        if quality_inlet == -1.0:
            density_inlet = PropsSI('D', 'T', temperature_inlet, 'P', pressure_inlet, fluid) 
            enthalpy_inlet = PropsSI('H', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
            viscosity_inlet = PropsSI('V', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
        else:  
            density_inlet = PropsSI('D', 'T', temperature_inlet, 'Q', quality_inlet, fluid) 
            enthalpy_inlet = PropsSI('H', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
            viscosity_inlet = PropsSI('V', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
    #
    else:
        #
        # Pressure and quality are given
        quality_inlet = Q_inlet #[Pa]
        pressure_inlet = P_inlet #[K]
        temperature_inlet = PropsSI('T', 'P', pressure_inlet, 'Q', quality_inlet, fluid)
        # Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
        # se il flusso è bifase utilizziamo T e titolo di vapore.
        if quality_inlet == -1.0:
            density_inlet = PropsSI('D', 'T', temperature_inlet, 'P', pressure_inlet, fluid) 
            enthalpy_inlet = PropsSI('H', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
            viscosity_inlet = PropsSI('V', 'T', temperature_inlet, 'P', pressure_inlet, fluid)
        else:  
            density_inlet = PropsSI('D', 'T', temperature_inlet, 'Q', quality_inlet, fluid) 
            enthalpy_inlet = PropsSI('H', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
            viscosity_inlet = PropsSI('V', 'T', temperature_inlet, 'Q', quality_inlet, fluid)
    
        
    # OUTLET
    #
    if Q_outlet == None:
        #
        # Per il calcolo delle consizioni di uscita, si presuppone che il FLUIDO SECONDARIO subisca perdite di carico pari a quelle massime accettabili
        pressure_outlet = P_inlet - P_loss   # CAPIRE COME DEFINIRE LO STATO DEL FLUSSO
        temperature_outlet = T_inlet #[k]
        quality_outlet = PropsSI('Q', 'T', temperature_outlet, 'P', pressure_outlet, fluid)
        if quality_outlet == -1.0:
            density_outlet = PropsSI('D', 'T', temperature_outlet, 'P', pressure_outlet, fluid) 
            enthalpy_outlet = PropsSI('H', 'T', temperature_outlet, 'P', pressure_outlet, fluid)
            viscosity_outlet = PropsSI('V', 'T', temperature_outlet, 'P', pressure_outlet, fluid)
        else:
            density_outlet = PropsSI('D', 'T', temperature_outlet, 'Q', quality_outlet, fluid) 
            enthalpy_outlet = PropsSI('H', 'T', temperature_outlet, 'Q', quality_outlet, fluid)
            viscosity_outlet = PropsSI('V', 'T', temperature_outlet, 'Q', quality_outlet, fluid)
    
    elif T_outlet == None:
        #
        # Per il calcolo delle consizioni di uscita, si presuppone che il FLUIDO SECONDARIO subisca perdite di carico pari a quelle massime accettabili
        pressure_outlet = P_inlet - P_loss   # CAPIRE COME DEFINIRE LO STATO DEL FLUSSO
        quality_outlet = Q_outlet
        temperature_outlet = PropsSI('T', 'P', pressure_outlet, 'Q', quality_outlet, fluid)
        if quality_outlet == -1.0:
            density_outlet = PropsSI('D', 'T', temperature_outlet, 'P', pressure_outlet, fluid) 
            enthalpy_outlet = PropsSI('H', 'T', temperature_outlet, 'P', pressure_outlet, fluid)
            viscosity_outlet = PropsSI('V', 'T', temperature_outlet, 'P', pressure_outlet, fluid)
        else:
            density_outlet = PropsSI('D', 'T', temperature_outlet, 'Q', quality_outlet, fluid) 
            enthalpy_outlet = PropsSI('H', 'T', temperature_outlet, 'Q', quality_outlet, fluid)
            viscosity_outlet = PropsSI('V', 'T', temperature_outlet, 'Q', quality_outlet, fluid)


    return pressure_inlet, temperature_inlet, quality_inlet, density_inlet, enthalpy_inlet, viscosity_inlet, pressure_outlet, temperature_outlet, quality_outlet, density_outlet, enthalpy_outlet, viscosity_outlet