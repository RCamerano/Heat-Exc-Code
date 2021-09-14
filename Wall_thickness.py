# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 09:25:17 2021

@author: Riccardo
"""

import pandas as pd
from sklearn.ensemble import RandomForestRegressor

# Questa funzione serve a definire lo spessore dei tubi da 0.25, 0.5, 0.75 e 1 pollice.
# Dal momento che non tutti gli spessori sono disponibili andiamo a prendere lo spessore immediatamente
# successivo a quello calcolato.
def wall_thickness(Do,wt_min):
    # Conversione del diametro nominale da pollici a metri e successiva conversione in formato deciamle a 2 cifre significative
    Do = Do/0.0254
    wt_min = wt_min*1000
    Do = float("{0:.2f}".format(Do))
    wt_min = float("{0:.2f}".format(wt_min))
    if Do == 0.25:
        if wt_min < 1.65:
            wt = 1.65
        elif wt_min >= 1.65 and wt_min < 2.24:
            wt = 2.24
        elif wt_min >= 2.24 and wt_min < 3.02:
            wt = 3.02
        else:
            print('spessore tubi non coerente')
    elif Do == 0.5:
        if wt_min < 1.65:
            wt = 1.65
        elif wt_min >= 1.65 and wt_min < 2.11:
            wt = 2.11
        elif wt_min >= 2.11 and wt_min < 2.77:
            wt = 2.77
        elif wt_min >= 2.77 and wt_min < 3.73:
            wt = 3.73
        elif wt_min >= 3.73 and wt_min < 4.78:
            wt = 4.78
        elif wt_min >= 4.78 and wt_min < 7.47:
            wt = 7.47
        else:
            print('spessore tubi non coerente')
    elif Do == 0.75:
        if wt_min < 1.65:
            wt = 1.65
        elif wt_min >= 1.65 and wt_min < 2.11:
            wt = 2.11
        elif wt_min >= 2.11 and wt_min < 2.87:
            wt = 2.87
        elif wt_min >= 2.87 and wt_min < 3.91:
            wt = 3.91
        elif wt_min >= 3.91 and wt_min < 5.56:
            wt = 5.56
        elif wt_min >= 5.56 and wt_min < 7.82:
            wt = 7.82
        else:
            print('spessore tubi non coerente')
    else:
         if wt_min < 1.65:
            wt = 1.65
         elif wt_min >= 1.65 and wt_min < 2.77:
            wt = 2.77
         elif wt_min >= 2.77 and wt_min < 3.38:
            wt = 3.38
         elif wt_min >= 3.38 and wt_min < 4.55:
            wt = 4.55
         elif wt_min >= 4.55 and wt_min < 6.35:
            wt = 6.35
         elif wt_min >= 6.35 and wt_min < 9.09:
            wt = 9.09
         else:
            print('spessore tubi non coerente')
    
    # Conversione dello spessore da millimetri a metri.
    wt = wt/1000 # [m]
    
    return wt

# Calcolo dello spessore minimo del tubo date norme ASME date le condizioni di pressione esterne ed interne al tubo.
def minimum_wall_thickness(P,P_ext,eta_joint,R,S,L,T,E,data):
    # TUBE UNDER INTERNAL PRESSURE: Sezione UG-27
    # Definizione del delta_P critico sul quale effettuare il dimensionamento del tubo quando questo è sotto lo
    # sforzo della pressione interna al tubo.
    P_critical_int = max(P)
    P_critical_ext = min(P_ext)
    dP = data.Selected_Value[17]*(P_critical_int - P_critical_ext)
    if dP > 0:
        # Circumferential stress
        t_min_CS = dP*R/(S*eta_joint - 0.6*dP)
        if dP > 0.385*S*eta_joint or t_min_CS > 0.5*R:
            print('WARNING: il modello a spessore sottile non è valido per il tubo considerato')
        # Longitudinal stress
        t_min_LS = dP*R/(2*S*eta_joint + 0.4*dP)
        if dP > 1.25*S*eta_joint or t_min_CS > 0.5*R:
            print('WARNING: il modello a spessore sottile non è valido per il tubo considerato')        
        t_min_internal = max([t_min_CS,t_min_LS])     
    else:
        t_min_internal = 0
             
    # TUBE UNDER EXTERNAL PRESSURE: Sezione UG-28 ASME Sez.8 Div.1
    # Definizione del delta_P critico sul quale effettuare il dimensionamento del tubo quando questo è sotto lo
    # sforzo della pressione esterna al tubo.
    P_critical_int = min(P)
    P_critical_ext = max(P_ext)
    dP = data.Selected_Value[17]*(P_critical_ext - P_critical_int)
    if dP > 0:
        # 'clf' è l'oggetto che richiama il modello RandomForest che verrà usata per l'interpolazione di tabelle
        # che sono state riportate dal manuale ASME.
        clf = RandomForestRegressor()
        
        # Il metodo iterativo per il calcolo dello spessore parte da un valore di primo tentativo, si calcola
        # la pressione massima consentita per evitare il cedimento dei tubi e si verifica se questa è maggiore o 
        # inferiore alla pressione esercitata sul tubo. Se la pressione consentita è inferiore a quella esercitata
        # si prende uno spessore superiore e si continua la verifica.
        # Step 1: definizione spessore di primo tentativo.
        t_aux = 0.0001 # RIVEDERE VALORE DI PRIMO TENTATIVO
        decision = 'not accepted'
        while decision == 'not accepted':
            # Step 2: definizione dei parametri per l'ingresso nelle charts
            param1 = L/(2*R)
            param2 = 2*R/t_aux
            if param1 > 50:
                param1 = 50
            elif param1 < 0.05:
                param1 = 0.05
            if param2 < 4 or param2 > 1000:
                print('WARNING: D/t è fuori dal range consigliato da ASME. Estrapolazione non è permessa')
            
            if param2 >= 10:
                # Step 3-4-5: Otteniamo i valori dei parametri A e B utilizando le charts ASME
                # L'algoritmo 'RandomForest' è utilizzato per interpolare le tabelle ASME per i parametri A e B
                data = pd.read_excel (r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\Pipe_sizing_A.xlsx')
                clf.fit(data[['D/t','L/D']], data['A'])
                A = clf.predict([[param2,param1]])
                data = pd.read_excel (r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\Pipe_sizing_B.xlsx')
                clf.fit(data[['T','A']], data['B [Mpa]'])
                B = clf.predict([[T,A]])
                
                # Step 6-7: determinazione della massima pressione esterna consentita
                if A > 0.0002:
                    P_a = 4*B*1E6 / 3 / param2
                else:
                    P_a = 2*A*E / 3 / param2
                
                # Step 8: comparazione di P_a con P
                if P_a > dP:
                    decision = 'accepted'
                    print('Allowable pressure is', float(P_a/dP), 'times the exerted pressure')
                else:
                    t_aux = 1.5*t_aux
            
            elif param2 < 10 and param2 > 4:
                # Step 3-4-5: Otteniamo i valori dei parametri A e B utilizando le charts ASME
                # L'algoritmo 'RandomForest' è utilizzato per interpolare le tabelle ASME per i parametri A e B
                data = pd.read_excel (r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\Pipe_sizing_A.xlsx')
                clf.fit(data[['D/t','L/D']], data['A'])
                A = clf.predict([[param2,param1]])
                data = pd.read_excel (r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\Pipe_sizing_B.xlsx')
                clf.fit(data[['T','A']], data['B [Mpa]'])
                B = clf.predict([[T,A]])
                
                # Step 6-7: determinazione della massima pressione esterna consentita
                P_a1 = ((2.167/param2) - 0.0833) * B*1E6 # [Pa]
                P_a2 = 2*S/param2 * (1 - 1/param2) # [Pa]
                P_a = min([P_a1,P_a2]) # [Pa]
                
                # Step 8: comparazione di P_a con P
                if P_a > dP:
                    decision = 'accepted'
                else:
                    t_aux = 1.5*t_aux
                
            elif param2 < 4:
                A = 1.1/(param2)**2
                # L'algoritmo 'RandomForest' è utilizzato per interpolare le tabelle ASME per il parametro B
                data = pd.read_excel (r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\Pipe_sizing_B.xlsx')
                clf.fit(data[['T','A']], data['B [Mpa]'])
                B = clf.predict([[T,A]])
                
                # Step 6-7: determinazione della massima pressione esterna consentita
                P_a1 = ((2.167/param2) - 0.0833) * B * 1E6 # [Pa]
                P_a2 = 2*S/param2 * (1 - 1/param2) # [Pa]
                P_a = min([P_a1,P_a2]) # [Pa]
                
                # Step 8: comparazione di P_a con P
                if P_a > dP:
                    decision = 'accepted'
                else:
                    t_aux = 1.5*t_aux
            
        t_min_external = t_aux
    else:
        t_min_external = 0
    
    # Lo spessore minimo del tubo è il massimo tra gli spessori necessari a reggere la pressione interna ed esterna.
    t_min = max(t_min_external,t_min_internal)
                
    
    return t_min
    