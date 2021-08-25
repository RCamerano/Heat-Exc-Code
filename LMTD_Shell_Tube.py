# HEAT EXCHANGER MAIN SCRIPT
# CDoICE PER DIMENSIONAMENTO SCAMBIATORI SHELL&TUBE CON METDoO (F)LMTD

import math
import Ft_calc
import TEMA_check
import numpy as np
import FlowBC
import Wall_thickness as wt
from CoolProp.CoolProp import PropsSI

# Sez. 1 - DISPOSIZIONE FLUIDI
# All'utente è chiesto di specificare se il fluido di lavoro scorre nel tubo interno (i.e. 'pipe') o
# nel mantello esterno (i.e. 'shell'). Se viene utilizzato un input diverso il cDoice chiede di inserire nuovamente l'input.
# Verrà attribuito l'indice 'wf_position = 0' al fluido di lavoro se questo è all'interno ('pipe') o 'wf_position = 1' se è all'esterno ('shell').
# In tutto il resto del cDoice il primo elemento di una lista farà riferimento al fluido che scorre all'interno, il secondo elemento a 
# quello che scorre all'esterno, indipendentemente dal fatto che il fluido di lavora sia nel pipe o nello shell.
# In tutto il cDoice le proprietà relative al fluido di lavoro sono caratterizate dall'indice 'wf_position'
# In tutto il cDoice le proprietà relative al fluido secondario sono caratterizate dall'indice 'sf_position'
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
[pressure_inlet[wf_position], temperature_inlet[wf_position], quality_inlet[wf_position], density_inlet[wf_position], enthalpy_inlet[wf_position], viscosity_inlet[wf_position], pressure_outlet[wf_position], temperature_outlet[wf_position], quality_outlet[wf_position], density_outlet[wf_position], enthalpy_outlet[wf_position], viscosity_outlet[wf_position]] = FlowBC.InletOutlet(fluid[wf_position], P_loss[wf_position], P_inlet = 100000, T_inlet = 300, T_outlet = 330)

# Sez. 5 - FLUIDO SECONDARIO
#
# Inserimento delle condizioni del FLUIDO SECONDARIO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[sf_position] = 'Water'
[pressure_inlet[sf_position], temperature_inlet[sf_position], quality_inlet[sf_position], density_inlet[sf_position], enthalpy_inlet[sf_position], viscosity_inlet[sf_position], pressure_outlet[sf_position], temperature_outlet[sf_position], quality_outlet[sf_position], density_outlet[sf_position], enthalpy_outlet[sf_position], viscosity_outlet[sf_position]] = FlowBC.InletOutlet(fluid[sf_position], P_loss[sf_position], P_inlet = 100000, T_inlet = 300, T_outlet = 330)

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
Uo = input("Estimate overall heat exchange coefficient Uo [W/m2 K]:\n")

# Calcolo Area di Scambio
A = Q / (Uo * LMTD * Ft)

# Sez. 8 - SCELTE TUBI

# Scelta lunghezza tubi (L) #mm
L = input('Specify tube lenght [mm]:\n')

# Scelta diametro nominale tubi (Do)
condition = 'false'
while condition == 'false':
    Do = input('Specify tubes outer diameter (0.5 - 0.75 - 1) [inc]:\n')
    if Do == 0.5 or Do == 0.75 or Do == 1:
        condition = 'true'
    else:
        print('diametro esterno tubi non standard. Scegliere uno dei valori suggeriti')

# Calcolo minimo spessore tubi (wt_min) secondo norma ASME B31.3 --> wt_min = (P.OD) / 2(SEW + PY) + CA
# P = massimo gradiente di pressione che si può verificare nello scambiatore al design point. Eventuali condizioni
# più severe verranno considerate scegliendo un opportuno fattore di sicurezza
P = max(pressure_inlet[wf_position],pressure_inlet[sf_position],pressure_outlet[wf_position],pressure_outlet[sf_position],1/pressure_inlet[wf_position],1/pressure_inlet[sf_position],1/pressure_outlet[wf_position],1/pressure_outlet[sf_position],abs(pressure_inlet[wf_position]-pressure_inlet[sf_position]),abs(pressure_outlet[wf_position]-pressure_outlet[sf_position]))
# OD = diametro esterno (outer diameter) secondo norme ASME
if Do == 0.5:
    OD = 0.84 #inch
elif Do == 0.75:
    OD = 1.05 #inch
elif Do == 1:
    OD = 1.315 #inch
else:
    print('errore nella scelta del diametro tubi')
# S = Maximum allowable stress value for the selected material [MPa]. Riferirsi a tabella 1.A ASME B36.10 per carbon steel e ASME B36.19 per stainless steel
# Alcuni valori sono riportati qui: https://www.cis-inspector.com/asme-code-calculation-allowable-stresses.html
S = 118 # riferito a carbon steel per tubature di piccole dimensioni, T < 300 °C

# E = Longitudinal Weld Joint Quality factor (1 per seamless pipe, 0.6 per Furnace Butt Welded Pipes, 0.85 per Electric Resistance Welded Pipes)
E = 1

# W = Weld Joint Strength Reduction Factor (incluso solo per completezza. Per noi è 1)
W = 1

# Y = Coefficiente da Tabella 304.1.1. Per tubature con spessore "sottile". Consigliato di applicarlo sempre. Per temperature < 500 °C, Y = 0.4
Y = 0.4

# CA = corrosion allowance
CA = 0.5 #mm -  Valore esempio. Scegliere opportunamente

# Sf = Safety Factor
Sf = 1.1

wt_min = ((P*OD*0.0254) / 2 / (S*1000000*E*W + P*Y) * 0.001) * Sf + CA #mm

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
