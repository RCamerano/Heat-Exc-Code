# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 15:43:12 2021

@author: Paolo
"""

def shell_sizing(tube_disposition,np):
    # Triangular disposition
    if tube_disposition == 't':
        if np == 1:
            K1 = 0.319
            n = 2.142
        elif np == 2:
            K1 = 0.249
            n = 2.207
        elif np == 4:
            K1 = 0.175
            n = 2.285
        elif np == 6:
            K1 = 0.0743
            n = 2.499
        elif np == 8:
            K1 = 0.0365
            n = 2.675
    # Square disposition
    if tube_disposition == 's':
        if np == 1:
            K1 = 0.215
            n = 2.207
        elif np == 2:
            K1 = 0.156
            n = 2.291
        elif np == 4:
            K1 = 0.158
            n = 2.263
        elif np == 6:
            K1 = 0.0402
            n = 2.617
        elif np == 8:
            K1 = 0.0331
            n = 2.643
    
    return K1,n

def clearance(shell_type,D_bundle):
    # Le variabili considerate in questa funzione sono:
    # - x: diametro del bundle
    # - y: shell clearance
    # Effettuiamo la linearizzazione dei grafici in cui sono tabellati i valori dello shell clearance
    x1 = 1.2 # [m]
    x2 = 0.2 # [m]
    if shell_type == 'PT':
        y1 = 87 # [mm]
        y2 = 96 # [mm]
    elif shell_type == 'SR':
        y1 = 50 # [mm]
        y2 = 78 # [mm]
    elif shell_type == 'P':
        y1 = 38 # [mm]
        y2 = 38 # [mm]
    elif shell_type == 'U':
        y1 = 10 # [mm]
        y2 = 20 # [mm]
        
    m = (y1-y2)/(x1-x2)
    b = (x1*y2 - x2*y1)/(x1-x2)
    
    shell_clearance = 0.001*(m*D_bundle + b) # [m]
    
    return shell_clearance