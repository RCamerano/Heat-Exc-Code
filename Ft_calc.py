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
       condition = 'true'
    else:
        print('per ora è gestibile solo la condizione tubi-unmixed and shell-mixed. Modificare flow_mix coerentemente')
    return Ft

    