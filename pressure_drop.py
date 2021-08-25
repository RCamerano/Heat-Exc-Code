# -*- coding: utf-8 -*-
"""
Created on Thu May 27 15:03:47 2021

@author: Paolo
"""
import numpy as np
import scipy.optimize
import pandas as pd

def pressure_loss(HE_CV,delta_l):
    # Calcolo delle perdite di carico distribuite all'interno del volume di controllo dovute a variazione di densita
    # e attrito. Vengono richiamate le due funzioni 'accelleration_loss' e 'friction_loss' che forniscono la perdita di carico
    # sotto forma di metri di colonna liquida persa.
    delta_p = [0,0]
    delta_p[0] = (0.5*(HE_CV.velocity[0]*HE_CV.density_in[0])**2)*(accelleration_loss(HE_CV)[0] + friction_loss(HE_CV,delta_l)[0])
    delta_p[1] = (0.5*(HE_CV.velocity[1]*HE_CV.density_in[1])**2)*(accelleration_loss(HE_CV)[1] + friction_loss(HE_CV,delta_l)[1])
    
    return delta_p
    
def concentrated_loss(HE_CV,A_fitting):
    concentrated_delta_p = [0,0]
    # Inizializzazione della variabile sigma che contiene il rapporto tra area di ingresso e area di attraversamento
    # nello scambiatore.
    sigma_inlet = [0,0]
    for i in [0,1]:
        # La grandezza sigma è definita minore di 1, è dunque opportuno calcolarla opportunamente.
        if HE_CV.A_cross[i] <= A_fitting[i]:
            sigma_inlet[i] = HE_CV.A_cross[i]/A_fitting[i]
        else:
            sigma_inlet[i] = A_fitting[i]/HE_CV.A_cross[i]
    
    # Lettura del file contenente i valori tabulati di K_c
    df = pd.read_csv('concentrated_velocity_heads.csv', index_col=0)
    # Valutazione parte costante di K. Il tipo di ingresso/uscita alo scambiatoe è specificato nel secondo indice della
    # variabile 'df'
    K_c = [(df['K_c']['elbow_45_standard']),(df['K_c']['elbow_45_standard'])]
    # Valutazione K dovuto a strettoie. Questo valore dipende dal rapporto tra l'area di attraversamento allo scambiatore
    # e l'area di ingresso/uscita
    for i in [0,1]:
        if HE_CV.A_cross[i] <= A_fitting[i]:
            K_c[i] += 0.5*(1-sigma_inlet[i])
    
    concentrated_delta_p[0] = (K_c[0] + 1 - sigma_inlet[0]**2)/HE_CV.density_in[0]
    concentrated_delta_p[1] = (K_c[1] + 1 - sigma_inlet[1]**2)/HE_CV.density_in[1]
    # Check sul moto del fluido in ingresso allo scambiatore
    if HE_CV.Reynolds[0] < 2300 or HE_CV.Reynolds[1] < 2300:
        print('Laminar flow at heat exchanger inlet/outlet. Pressure losses model not reliable')
            
    return concentrated_delta_p

def accelleration_loss(HE_CV):
    # Calcolo delle perdite di carico dovute alla variazione di densità all'interno dello scambiatore
    accelleration_delta_p = [0,0]
    accelleration_delta_p[0] = 2*(1/HE_CV.density_out[0]-1/HE_CV.density_in[0])
    accelleration_delta_p[1] = 2*(1/HE_CV.density_out[1]-1/HE_CV.density_in[1])
    
    return accelleration_delta_p
    
def friction_loss(HE_CV,delta_l):
    # calcolo delle perdite di carico dovute all'attrito con le pareti dello scambiatore
    friction_delta_p = [0,0]
    # il valore del friction coefficient è calcolato all'interno della funzione 'friction_coeff'
    f = friction_coeff(HE_CV)
    # Equazione di Fanning
    friction_delta_p[0] = 4*f[0]*delta_l/(np.mean([HE_CV.density_in[0],HE_CV.density_out[0]]))/HE_CV.D_hydraulic[0]
    friction_delta_p[1] = 4*f[1]*delta_l/(np.mean([HE_CV.density_in[1],HE_CV.density_out[1]]))/HE_CV.D_hydraulic[1]
    
    return friction_delta_p

def friction_coeff(HE_CV):
    f = [0,0]
    epsilon = 0.00005 # [m] reference Moody
    for i in [0,1]:
        if HE_CV.Reynolds[i] > 4000: # Colebrook
            # Formula semplificata per il calcolo del friction coefficient proposta da Churchill.
            f[i] = (-4*np.log10(0.27*epsilon/HE_CV.D_hydraulic[i] + (7/HE_CV.Reynolds[i])**0.9))**-2
            # Per il calcolo del friction coefficient è necessario risolvere un'equazione non-lineare. Per farlo richiamiamo
            # il risolutore 'fsolve'. La variabile incognita 'A' è uguale alla radice quafrata di 'f'.
            # Unknown_A_start = np.sqrt(0.007) # reference Riccardo
            # Objective_function_A = lambda A : 1/A + 4*np.log10(epsilon/(3.7*HE_CV.D_hydraulic[i]) + 1.256/(HE_CV.Reynolds[i]*A))
            # A = scipy.optimize.fsolve(Objective_function_A,Unknown_A_start,xtol = 1.0e-8,maxfev = 10^5)[0]
            # f[i] = A**2
        elif HE_CV.Reynolds[i] < 2100: # Hagen-Poiseuille
            # Formula di Hagen-Poiseuille per il calcolo del friction coefficient in flussi laminari.
            f[i] = 16/HE_CV.Reynolds[i]
        else: # Churchill
            # Equazione di Churchill per il calcolo del friction coefficient in regime di moto transitorio.
            A = (2.457*np.log(1/((7/HE_CV.Reynolds[i])**0.9)+0.27*epsilon/HE_CV.D_hydraulic[i]))**16
            B = (37530/HE_CV.Reynolds[i])**16
            f[i] = 2*(((8/HE_CV.Reynolds[i])**12) + 1/(A+B)**1.5)**(1/12)
            
    return f