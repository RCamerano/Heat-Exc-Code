# HEAT EXCHANGER MAIN SCRIPT
# CDoICE PER DIMENSIONAMENTO SCAMBIATORI SHELL&TUBE CON METDoO (F)LMTD

import math
import Ft_calc
import TEMA_check
import numpy as np
import pandas as pd
import FlowBC
import Wall_thickness as wt
import shell
import heat_transfer as HT
from CoolProp.CoolProp import PropsSI

# Sez. 1 - DISPOSIZIONE FLUIDI
# All'utente è chiesto di specificare se il fluido di lavoro scorre nel tubo interno (i.e. 'pipe') o
# nel mantello esterno (i.e. 'shell'). Se viene utilizzato un input diverso il cDoice chiede di inserire nuovamente l'input.
# Verrà attribuito l'indice 'wf_position = 0' al fluido di lavoro se questo è all'interno ('pipe') o 'wf_position = 1' se è all'esterno ('shell').
# In tutto il resto del cDoice il primo elemento di una lista farà riferimento al fluido che scorre all'interno, il secondo elemento a 
# quello che scorre all'esterno, indipendentemente dal fatto che il fluido di lavora sia nel pipe o nello shell.
# In tutto il codice le proprietà relative al fluido di lavoro sono caratterizate dall'indice 'wf_position'
# In tutto il codice le proprietà relative al fluido secondario sono caratterizate dall'indice 'sf_position'
decision = 'not accepted'
while decision != 'accepted':
    wf_position = input("Specify working fluid position in the heat exchanger ['pipe' or 'shell']:\n")
    if wf_position == 'pipe':
        wf_position = 0
        sf_position = 1
        decision = 'accepted'
    elif wf_position == 'shell':
        wf_position = 1
        sf_position = 0
        decision = 'accepted'
    else:
        print("input is not available \n")
        
# Sez. 2 - INIZIALIZZAZIONE LISTE
# Dichiarazione variabili per le condizioni al contorno e le costanti del problema
# Ogni variabile è contenuta all'interno di una lista. Il primo elemento delle liste farà riferimento al fluido che scorre nel tubo interno,
# il secondo all'elemento a quello che scorre nel mantello esterno.
# Le condizioni al contorno del problema sono le seguenti:
    # - Temperatura in ingresso ed uscita dello scambiatore per i due fluidi.
    # - Pressione in ingresso ed uscita dello scambiatore per i due fluidi.
    # - Titolo di vapore in ingresso ed uscita dello scambiatore per i due fluidi.
    # - Densità in ingresso ed uscita dello scambiatore per i due fluidi.
    # - Entalpia in ingresso ed uscita dello scambiatore per i due fluidi.
    # - Velocità dei flussi in ingresso e uscita, come conseguenza delle condizioni di flusso e delle scelte ingegneristiche
    # - Numero di Reynolds in ingresso e uscita
# Sono costanti in tutto lo scambiatore i due fluidi e le loro portate massiche.
M = [0,0]
fluid = [0,0]
pressure_inlet = [0,0]
pressure_outlet = [0,0]
temperature_inlet = [0,0]
temperature_outlet = [0,0]
quality_inlet = [0,0]
quality_outlet = [0,0]
density_inlet = [0,0]
density_outlet = [0,0]
enthalpy_inlet = [0,0]
enthalpy_outlet = [0,0]
velocity_inlet = [0,0]
velocity_outlet = [0,0]
viscosity_inlet = [0,0]
viscosity_outlet = [0,0]
Re_inlet = [0,0]
Re_outlet = [0,0]
j_h = [0,0]
alfa = [0,0]
h_fouling = [0,0]

# Sez. 3 - SCELTA MASSIME PERDITE DI CARICO ACCETTABILI
P_loss = [0,0]
P_loss[wf_position] = int(input("Specify maximum allowable pressure loss - WORKING FLUID [Pa]:\n"))
P_loss[sf_position] = int(input("Specify maximum allowable pressure loss - SECONDARY FLUID [Pa]:\n"))

# Sez. 4 - FLUIDO DI LAVORO
#
# Inserimento delle condizioni del FLUIDO DI LAVORO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[wf_position] = 'Water'
[pressure_inlet[wf_position], temperature_inlet[wf_position], quality_inlet[wf_position], density_inlet[wf_position], enthalpy_inlet[wf_position], viscosity_inlet[wf_position], pressure_outlet[wf_position], temperature_outlet[wf_position], quality_outlet[wf_position], density_outlet[wf_position], enthalpy_outlet[wf_position], viscosity_outlet[wf_position]] = FlowBC.InletOutlet(fluid[wf_position], P_loss[wf_position], P_inlet = 200000, T_inlet = 363, T_outlet = 343)

# Sez. 5 - FLUIDO SECONDARIO
#
# Inserimento delle condizioni del FLUIDO SECONDARIO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[sf_position] = 'Water'
[pressure_inlet[sf_position], temperature_inlet[sf_position], quality_inlet[sf_position], density_inlet[sf_position], enthalpy_inlet[sf_position], viscosity_inlet[sf_position], pressure_outlet[sf_position], temperature_outlet[sf_position], quality_outlet[sf_position], density_outlet[sf_position], enthalpy_outlet[sf_position], viscosity_outlet[sf_position]] = FlowBC.InletOutlet(fluid[sf_position], P_loss[sf_position], P_inlet = 100000, T_inlet = 293, T_outlet = 303)

# DEFINIZIONE MASSA FLUSSO SECONDARIO
# La portata massica del fluido secondario viene calcolata sulla base delle condizioni al contorno specificate per i due fluidi
# in mDoo tale che il bilancio energetico sia sempre rispettato.
M[sf_position] = abs(M[wf_position] * (enthalpy_inlet[wf_position] - enthalpy_outlet[wf_position]) / (enthalpy_inlet[sf_position] - enthalpy_outlet[sf_position]))


# Sez. 6 - SCELTE GENERALI DI DESIGN
condition = 'false'
# Selezionare categoria TEMA dello scambiatore (HE_cat)
HE_cat = input("Specify heat exchanger TEMA designation [E,F,X,K]:\n")
# Selezionare tipologia dello scambiatore Shell & Tube (HE_type)
HE_type = input("Specify heat exchanger TEMA designation [1 = fixed-tube plate   2 = U-tube   ]:\n")
while condition == 'false':
    # Specificare numero passaggi lato tubi (np_tubi) e numero passaggi lato mantello (np_shell)
    np_tubi = int(input("Tube Passes number: \n"))
    np_shell = int(input("Shell Passes number: \n"))
    condition = TEMA_check.TEMA_check(HE_cat,np_shell,np_tubi)
    condition = TEMA_check.Ft_check(HE_cat,np_tubi)

# Sez. 7 - DEFINIZIONE PARAMETRI FONDAMENTALI

# Calcolo potenza termica scambiata nello scambiatore
Q = abs(M[wf_position] * (enthalpy_inlet[wf_position] - enthalpy_outlet[wf_position]))

# Calcolo differenza di temperatura medio logaritmica (LMTD) corrispondente a perfetta controcorrente
Delta_T1 = temperature_inlet[wf_position] - temperature_outlet[sf_position]
Delta_T2 = temperature_outlet[wf_position] - temperature_inlet[sf_position]
if Delta_T1 == Delta_T2:
    LMTD = Delta_T1
else:
    LMTD = (Delta_T1 - Delta_T2) / np.log(Delta_T1/Delta_T2)
    
# Calcolo fattore di correzione (Ft) per non perfetta controcorrente
if HE_cat == 'E':
    Ft = Ft_calc.Ft_TEMA_E(temperature_inlet[wf_position], temperature_outlet[wf_position], temperature_inlet[sf_position], temperature_outlet[sf_position],np_shell,np_tubi)
elif HE_cat == 'F':
    Ft = Ft_calc.Ft_TEMA_F()
elif HE_cat == 'K':
    Ft = Ft_calc.Ft_TEMA_K()
elif HE_cat == 'X':
    condition = 'false'
    while condition == 'false':
        flow_mix = [0,0]
        flow_mix[0] = input("Is pipe flow mixed? [Y/N]:\n")
        flow_mix[1] = input("Is shell flow mixed? [Y/N]:\n")
        Ft = Ft_calc.Ft_TEMA_X(flow_mix,temperature_inlet[wf_position], temperature_outlet[wf_position], temperature_inlet[sf_position], temperature_outlet[sf_position]) 
       
# Specificare il valore del coefficiente globale di scambio di primo tentativo (Uo)
# Per stimare il valore di Uo vedi Coulson and Richardson's pag. 637 o https://www.engineersedge.com/thermodynamics/overall_heat_transfer-table.htm
Uo = float(input("Estimate overall heat exchange coefficient Uo [W/m2 K]:\n"))

# Calcolo Area di Scambio
A = Q / (Uo * LMTD * Ft)

# Sez. 8 - SCELTE TUBI

# Calcolo del tube pitch (distanza tra il centro dei tubi) e scelta sulla disposizione dei tubi.
# Comunemente si ha pitch pari a 1.25, 1.33 o 1.5 diametro esterno del tubo.
tube_disposition = input("Select tube disposition between triangular [t] or square [s]:\n")
pitch_ratio = float(input("select the ratio between tube pitch and tube diamater [-]: \n"))

# Scelta lunghezza tubi (L) #mm
L = float(input('Specify tube lenght [mm]:\n'))

# Scelta diametro nominale tubi (Do)
condition = 'false'
while condition == 'false':
    Do = float(input('Specify tubes nominal diameter (0.5 - 0.75 - 1) [inc]:\n'))
    if Do == 0.5 or Do == 0.75 or Do == 1:
        condition = 'true'
    else:
        print('diametro nominale tubi non standard. Scegliere uno dei valori suggeriti')

# Calcolo minimo spessore tubi (wt_min) secondo norma ASME B31.3 --> wt_min = (P.OD) / 2(SEW + PY) + CA
# P = massimo gradiente di pressione che si può verificare nello scambiatore al design point. Eventuali condizioni
# più severe sono considerate scegliendo un opportuno fattore di sovradimensionamento
# Un fattore di sicurezza 'SF' è utilizzato come margine di sicurezza per la progettazione.
P = np.array([pressure_inlet[0],pressure_outlet[0],1])
P = float(input('select oversizing factor for pressure: \n'))*P
# OD = diametro esterno (outer diameter) secondo norme ASME
if Do == 0.5:
    OD = 0.84 #inch
elif Do == 0.75:
    OD = 1.05 #inch
elif Do == 1:
    OD = 1.315 #inch
else:
    print('errore nella scelta del diametro tubi')

pitch = pitch_ratio*OD
    
# S = Maximum allowable stress value for the selected material [MPa]. Riferirsi a tabella 1.A ASME B36.10 per carbon steel e ASME B36.19 per stainless steel
# Alcuni valori sono riportati qui: https://www.cis-inspector.com/asme-code-calculation-allowable-stresses.html
S = 170 # [MPa] 0.9*Yield strength of Carbon steel at 150°C

E = 204E3  #[MPa] Modulo di elasticità

T_wall = np.mean([np.mean([temperature_inlet[0],temperature_outlet[0]]),np.mean([temperature_inlet[1],temperature_outlet[1]])])
k = 8.54 + 0.02*T_wall # Conducibilità del materiale

# E = Efficiency of ligaments
eta_joint = (pitch - OD)/pitch

# CA = corrosion allowance
CA = 0.5 #mm -  Valore esempio. Scegliere opportunamente

# Sf = Safety Factor
Sf = 1.1

# La formula ASME prevede l'utilizzo del raggio interno. Per ragioni di procedura utiliziamo il raggio nominale del tubo
# (sovradimensionamento dello spessore minimo).
P_ext = [pressure_inlet[1],pressure_outlet[1],1]
T = np.mean([np.mean([temperature_inlet[0],temperature_outlet[0]]),np.mean([temperature_inlet[1],temperature_outlet[1]])])
wt_min = wt.minimum_wall_thickness(P,P_ext,eta_joint,Do/2,S,L,T,E) * Sf + CA #mm

# Selezione schedula tubi e rispettivo spessore (wt). Calcolo diametro interno (ID) #mm
wt = wt.wall_thickness(Do,wt_min)
ID = OD - 2*wt

# Sez. 9 - CALCOLI TUBI

# Calcolo n° tubi (nt), arrotondato al successivo numero pari
nt_rough = A / (np.pi * OD * L)
nt = math.ceil(nt_rough / 2.) * 2

# Calcolo velocità lato tubi in ingresso e uscita lato tubi #m/s
velocity_inlet[0] = (4 * M[0] * np) / (np.pi * density_inlet[0] * ID**2 * nt)
velocity_outlet[0] = (4 * M[0] * np) / (np.pi * density_outlet[0] * ID**2 * nt)

# Calcolo numero di Reynolds in ingresso e uscita lato tubi #-
Re_inlet[0] = density_inlet[0] * velocity_inlet[0] * ID * 10**3 / viscosity_inlet[0]
Re_outlet[0] = density_outlet[0] * velocity_outlet[0] * ID * 10**3 / viscosity_outlet[0]
if Re_inlet[0] < 10000:
    print("Attenzione: basso Reynolds all'ingresso dello scambiatore lato tubi!")
if Re_outlet[0] < 10000:
    print("Attenzione: basso Reynolds all'uscita dello scambiatore lato tubi!")

# Sez. 10 - SCELTE MANTELLO

# Calcolo del diametro dello shell. 'K1' ed 'n' sono parametri tabulati.
K1,n = shell.shell_sizing(tube_disposition,np_tubi)
D_bundle = OD*(nt/K1)**(1/n)

# Calcolo del bundle clearance basato su approccio grafico.
# Bisogna inizialmente fissare la tipologia di scambiatore scegliendo tra:
    # - Pull-through Floating head
    # - Split-ring Floating head
    # - Packed head
    # - U-tube
rear_head_type = input("Select between Pull-through floating head [PT], Split-ring floating head [SR],\n Packed head [P] and U-tube [U] :\n")
shell_clearance = shell.clearence(rear_head_type,D_bundle)
D_shell = D_bundle + 2*shell_clearance

# Per diametri dello shell inferiori a 24", lo shell è ottenuto da un tubo commerciale.
# I tubi commerciali non hanno diametri che variano con continuità, ma sono disponibili solo determinate misure.
# Prendiamo queste misure da tabella importata ('tubi.xlsx') e selezioniamo il valore superiore a 
# 'D_shell' più prossimo rispetto a quello calcolato.
if D_shell < 24*25.4: # conversione da inch a mm
    # Dimensioni caratteristiche di tubi commerciali sono nella tabella 'tubi.xlsx' 
    data = pd.read_excel (r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\tubi.xlsx')
    D_shell = max(data.iloc[(data['DN [mm]']-D_shell).abs().argsort()[:2]]['DN [mm]'].tolist())
    idx = max(data.iloc[(data['DN [mm]']-D_shell).abs().argsort()[:2]].index.tolist())
    D_shell_ext = data.D_ext[idx]
    P_ext = [1]
    P = [pressure_inlet[1],pressure_outlet[1],1]
    T = np.mean([temperature_inlet[1],temperature_outlet[1]])
    
    # Dimensionamento dello spessore dei tubi - Procedure ASME
    wt_min = wt.minimum_wall_thickness(P,P_ext,eta_joint,Do/2,S,L,T,E) * Sf + CA #mm
    
    # Utiliziamo come spessore del mantello uno spessore standard tra quelli presenti nella tabella 'tubi.xlsx'.
    wt_serie = data.loc[idx,'XS':'sch.160']
    wt = np.array(wt_serie)
    aux = wt - wt_min
    b = np.argwhere(aux > 0)
    wt = wt_min + min(aux[b])
    
    # Noto il diametro interno dello shell possiamo calcolare il corretto numero di tubi e l'area di scambio.
    D_shell_int = D_shell_ext - 2*wt
    if np_tubi == 1:
        CTP = 0.93
    elif np_tubi == 2:
        CTP = 0.9
    else:
        CTP = 0.85
    if tube_disposition == 't':
        CL = np.sin(np.pi/3)
    else:
        CL = 1
    nt = (CTP/CL)*(D_shell_int**2/pitch**2)*np.pi/4
    A = np.pi * OD * L * nt
    
# Sez. 11 - CALCOLI MANTELLO
    
# Viene inserito dall'utente il distanziamento tra i baffle dello scambiatore. Il distanziamento è un valore
# compreso tra un quinto del diametro dello shell e il diametro dello shell
decision = 'not accepted'
while decision == 'not accepted':
    B = input("Selection of the ratio between baffle spacing and shell diameter:\n")
    if B > 0.2 and B < 1:
        decision = 'accepted'
    else:
        print('the value of the ratio has to be comprised between 0.2 and 1')
B = D_shell/B

# Il numero di baffle è basato sul loro distanziamento e la lunghezza totale dello scambiatore
n_baffle = math.ceil(L/B - 1)

# Calcolo dell'area di attraversamento lato mantello
A_shell = (pitch - OD)*D_shell_int*B / pitch

# Calcolo del numero di Reynolds lato mantello
velocity_inlet[1] = M[1] / (density_inlet[1]*A_shell)
velocity_outlet[1] = M[1] / (density_outlet[1]*A_shell)

if tube_disposition == 't':
    D_eq_shell = 1.1*(pitch**2 - 0.917*OD**2)/OD
else:
    D_eq_shell = 1.27*(pitch**2 - 0.785*OD**2)/OD

Re_inlet[1] = density_inlet[1] * velocity_inlet[1] * D_eq_shell * 10**3 / viscosity_inlet[1]
Re_outlet[1] = density_outlet[1] * velocity_outlet[1] * D_eq_shell * 10**3 / viscosity_outlet[1]

baffle_cut = input('select baffle cut [15-20-25-30-35-40-45]: \n')

# Sez. 12 - VALUTAZIONE HEAT TRANSFER FACTOR

# Lato tubi
j_h[0] = HT.heat_transfer_factor_tube(np.mean([Re_inlet[0],Re_outlet[0]]),L/ID)

# Lato shell
j_h[1] = HT.heat_transfer_factor_shell(np.mean([Re_inlet[1],Re_outlet[1]]),baffle_cut)

# Sez. 13 - VALUTAZIONE HEAT TRANSFER COEFFICIENT

# Calcolo alfa lato tubo
if quality_inlet[0] == quality_outlet[0]:
    # Raffreddamento/Riscaldamento
    alfa[0] = HT.heat_transfer_coefficient_1ph(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,fluid,Re_inlet,Re_outlet,j_h,ID,0)
elif quality_inlet[0] > quality_outlet[0]:
    # Condensazione
    alfa[0] = HT.heat_exchange_coefficient_condensation(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,quality_inlet,quality_outlet,ID,fluid,1,M)
else:
    # Evaporazione
    alfa[0] = HT.heat_exchange_coefficient_evaporation()
  
# Calcolo alfa lato mantello
if quality_inlet[1] == quality_outlet[1]:
    # Raffreddamento/Riscaldamento
    alfa[1] = HT.heat_transfer_coefficient_1ph(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,fluid,Re_inlet,Re_outlet,j_h,ID,0)
elif quality_inlet[1] > quality_outlet[1]:
    # Condensazione
    alfa[1] = HT.heat_exchange_coefficient_condensation(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,quality_inlet,quality_outlet,ID,fluid,1,M)
else:
    # Evaporazione
    alfa[1] = HT.heat_exchange_coefficient_evaporation()
    
# Sez. 14 - SPECIFICA COEFFICIENTI DI SPORCAMENTO

h_fouling[0] = input('Select the fouling coefficient in the pipe side: \n')
h_fouling[1] = input('Select the fouling coefficient in the shell side: \n')

# Sez. 15 - CALCOLO COEFFICIENTE GLOBALE DI SCAMBIO (REALE)

Uo_star = (1/alfa[1] + 1/h_fouling[1] + OD*np.log(OD/ID)/(2*k) + (1/alfa[0])*(OD/ID) + (1/h_fouling[0])*(OD/ID))**-1

# Sez. 16 -  VALUTAZIONE FRICTION COEFFICIENT

