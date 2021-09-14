# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 18:02:38 2021

@author: Riccardo
"""

import numpy as np

def Ft_TEMA_E(T1,T2,t1,t2,np_shell,np_tubi):
    if np_shell == 1 and np_tubi == 1:
        Ft = 1
    else:
        P1 = (T2 - T1) / (t1 - T1)
        R1 = (t1 - t2) / (T2 - T1)
        S = (R1**2 + 1)**0.5 / (R1 - 1)
        W = ((1 - P1*R1) / (1 - P1))**(1/np_shell)
        if R1 == 1:
            W1 = (np_shell - np_shell * P1) /  (np_shell - np_shell * P1 + P1)
            Ft = 2**0.5 * ((1 - W1) / W1) / np.log((W1/(1-W1) + 2**(-0.5)) / (W1/(1-W1) - 2**(-0.5))) #!! Controllare formula !!
        else:
            Ft = S * np.log(W) / np.log((1 + W - S + S * W) / (1 + W + S - S * W))
    return Ft

def Ft_TEMA_F():
    Ft = 1
    return Ft

def Ft_TEMA_K():
    Ft = 1   # Approssimazione accurata se il kettle svolge principalmente evaporazione (liquido in arrivo poco sottoraffreddato)
    return Ft

def Ft_TEMA_X(flow_mix,T1,T2,t1,t2):
    # Procedura tratta da Kern "Process Heat Transfer"
    K = (T1 - T2) / (T1 - t1)
    S = (t2 - t1) / (T1 - t1)
    R = K / S
    r_counter = (R - 1) / (np.log(1 - S) / (1 - R*S))
    r = S / (np.log(1 / (1 - S/K * np.log(1/(1-K)) )))
    if flow_mix[0] == 'N' and flow_mix[1] == 'Y':
       Ft = r / r_counter
    else:
        print('per ora è gestibile solo la condizione tubi-unmixed and shell-mixed. Modificare flow_mix coerentemente')
    return Ft


def Ft_Air_cooler(temperature_inlet,temperature_outlet,wf_position,sf_position,LMTD,np_tubi):
    # Calcolo del fattore di correzione Ft per non perfetta controcorrente
    A = ((temperature_inlet[wf_position] - temperature_outlet[wf_position])**2 + (temperature_inlet[sf_position] - temperature_outlet[sf_position])**2)**0.5
    B = ((temperature_inlet[wf_position] - temperature_outlet[sf_position])**0.5 + (temperature_outlet[wf_position] - temperature_inlet[sf_position])**(1/1.7))**1.7
    
    CLMTD = A/(1.7*np.log((A+B)/(B-A)))
    if np_tubi == 1:
        Ft = CLMTD/LMTD
    elif np_tubi == 2:
        CLMTD = 0.6*CLMTD + 0.4*LMTD
        Ft = CLMTD/LMTD
    elif np_tubi >= 3:
        Ft = 1
        
    return Ft