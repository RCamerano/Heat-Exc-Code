# HEAT EXCHANGER MAIN SCRIPT
# CDoICE PER DIMENSIONAMENTO SCAMBIATORI SHELL&TUBE CON METDoO (F)LMTD

import math
import Ft_calc
import TEMA_check
import numpy as np
import pandas as pd
import FlowBC
import Wall_thickness
import shell
import heat_transfer as HT
import friction_coefficient
from CoolProp.CoolProp import PropsSI

# Gli input per il dimensionamento sono presi dal file excel 'input.xlsx'. Nel file sono presenti le
# seguenti colonne: 
    # - Category
    # - Parameter: qui è contenuto il nome del parametro di input
    # - Permissible Values: alcuni input non possono variare arbitrariamente, devono essere inseriti dei valori prestabiliti.
    # - Unit: unità di misura del parametro
    # - Selected Value: è il valore del parametro che viene caricato dal codice nel suo funzionamento
data = pd.read_excel(r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\input.xlsx','Input_AC')

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
    wf_position = data.Selected_Value[0]
    if wf_position == 'pipe':
        wf_position = 0
        sf_position = 1
        decision = 'accepted'
    else:
        print("input is not available \n")
    
    
# Sez. 2 - INIZIALIZZAZIONE LISTE
#
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
    # - Numero di Reynolds in ingresso e uscita.
# Sono costanti in tutto lo scambiatore i due fluidi e le loro portate massiche.
M = [0,0] # [Kg/s]
fluid = [0,0] # [-]
pressure_inlet = [0,0] # [Pa]
pressure_outlet = [0,0] # [Pa]
pressure_mean = [0,0] # [Pa]
temperature_inlet = [0,0] # [K]
temperature_outlet = [0,0] # [K]
temperature_mean = [0,0] # [K]
quality_inlet = [0,0] # [-]
quality_outlet = [0,0] # [-]
quality_mean = [0,0] # [-]
density_inlet = [0,0] # [Kg/m^3]
density_outlet = [0,0] # [Kg/m^3]
density_mean = [0,0] # [Kg/m^3]
enthalpy_inlet = [0,0] # [J/Kg]
enthalpy_outlet = [0,0] # [J/Kg]
enthalpy_mean = [0,0] # [J/Kg]
velocity_inlet = [0,0] # [m/s]
velocity_outlet = [0,0] # [m/s]
velocity_mean = [0,0] # [m/s]
viscosity_inlet = [0,0] # [Pa*s]
viscosity_outlet = [0,0] # [Pa*s]
Re_inlet = [0,0] # [-]
Re_outlet = [0,0] # [-]
Re_mean = [0,0] # [-]
j_h = [0,0] # [-]
alfa = [0,0] # [W/(m^2.K)]
h_fouling = [0,0] # [W/(m^2.K)]
f = [0,0] # [-]
delta_P = [0,0] # [Pa]
wt = [0,0] # [m]
wt_min = [0,0] # [m]
A_cross = [0,0] # [m^2]
D_eq = [0,0] # [m]


# Sez. 3 - SCELTA MASSIME PERDITE DI CARICO ACCETTABILI
#
# Le massime perdite di carico accettabili sono parametri di input del problema e sono utilizzate per il
# calcolo delle proprietà all'output dello scambiatore.
P_loss = [0,0]
P_loss[wf_position] = data.Selected_Value[7]
P_loss[sf_position] = data.Selected_Value[13]


# Sez. 4 - FLUIDO DI LAVORO
#
# Inserimento delle condizioni del FLUIDO DI LAVORO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
M[wf_position] = data.Selected_Value[1] # [Kg/s]
fluid[wf_position] = 'Water'
[pressure_inlet[wf_position], temperature_inlet[wf_position], quality_inlet[wf_position], density_inlet[wf_position], enthalpy_inlet[wf_position], viscosity_inlet[wf_position], pressure_outlet[wf_position], temperature_outlet[wf_position], quality_outlet[wf_position], density_outlet[wf_position], enthalpy_outlet[wf_position], viscosity_outlet[wf_position]] = FlowBC.InletOutlet(fluid[wf_position], P_loss[wf_position], P_inlet = data.Selected_Value[2], T_inlet = data.Selected_Value[3], T_outlet = data.Selected_Value[4], Q_inlet = data.Selected_Value[5], Q_outlet = data.Selected_Value[6])
#
# Alcune relazioni presenti nel codice utilizzano il valore medio delle proprietà dei fluidi
temperature_mean[wf_position] = np.mean([temperature_inlet[wf_position],temperature_outlet[wf_position]])
pressure_mean[wf_position] = np.mean([pressure_inlet[wf_position],pressure_outlet[wf_position]])
enthalpy_mean[wf_position] = np.mean([enthalpy_inlet[wf_position],enthalpy_outlet[wf_position]])
density_mean[wf_position] = np.mean([density_inlet[wf_position],density_outlet[wf_position]])
quality_mean[wf_position] = np.mean([quality_inlet[wf_position],quality_outlet[wf_position]])


# Sez. 5 - FLUIDO SECONDARIO
#
# Inserimento delle condizioni del FLUIDO SECONDARIO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[sf_position] = 'Air'
temperature_inlet[sf_position] = data.Selected_Value[9]
pressure_inlet[sf_position] = data.Selected_Value[8]
pressure_outlet[sf_position] = 101325

# Sez. 6 - SCELTE GENERALI DI DESIGN
#
# Selezionare categoria TEMA dello scambiatore (HE_cat)
HE_cat = data.Selected_Value[19]
# Selezionare tipologia dello scambiatore Shell & Tube (HE_type)
HE_type = data.Selected_Value[20]

# Viene specificato il numero passaggi lato tubi (np_tubi) e numero passaggi lato mantello (np_shell)
condition = 'false'
while condition == 'false':
    np_tubi = data.Selected_Value[22]
    np_shell = data.Selected_Value[24]
    # In base alla designazione dello scambiatore si verifica che il numero di passaggi scelto è coerente
    condition = TEMA_check.TEMA_check(HE_cat,np_shell,np_tubi)
    condition = TEMA_check.Ft_check(HE_cat,np_tubi)


# Sez. 7 - DEFINIZIONE PARAMETRI FONDAMENTALI
#
# Calcolo potenza termica scambiata nello scambiatore
Q = abs(M[wf_position] * (enthalpy_inlet[wf_position] - enthalpy_outlet[wf_position]))    
       
# Specificare il valore del coefficiente globale di scambio di primo tentativo (Uo)
# Per stimare il valore di Uo vedi Coulson and Richardson's pag. 637 o https://www.engineersedge.com/thermodynamics/overall_heat_transfer-table.htm
Uo = data.Selected_Value[14] # [W/m^2.K]

# La temperature dell'aria in uscita dallo scambiatore è valutata sulla base del valore di primo tentativo di Uo.
# Sulla base della temperatura di uscita dell'aria è possibile calcolare le altre proprietà
temperature_outlet[sf_position] = FlowBC.air_cooler_outlet(Uo,temperature_inlet,temperature_outlet)
density_outlet[sf_position] = PropsSI('D', 'T', temperature_outlet[sf_position], 'P', pressure_outlet[sf_position], fluid[sf_position])
enthalpy_outlet[sf_position] = PropsSI('H', 'T', temperature_outlet[sf_position], 'P', pressure_outlet[sf_position], fluid[sf_position])
enthalpy_mean[sf_position] = np.mean([enthalpy_inlet[sf_position],enthalpy_outlet[sf_position]])
density_mean[sf_position] = np.mean([density_inlet[sf_position],density_outlet[sf_position]])

# DEFINIZIONE MASSA FLUSSO SECONDARIO
# La portata massica del fluido secondario viene calcolata sulla base delle condizioni al contorno specificate per i due fluidi
# in mDoo tale che il bilancio energetico sia sempre rispettato.
M[sf_position] = abs(M[wf_position] * (enthalpy_inlet[wf_position] - enthalpy_outlet[wf_position]) / (enthalpy_inlet[sf_position] - enthalpy_outlet[sf_position])) # [Kg/s]

# Calcolo differenza di temperatura medio logaritmica (LMTD) corrispondente a perfetta controcorrente
Delta_T1 = temperature_inlet[wf_position] - temperature_outlet[sf_position]
Delta_T2 = temperature_outlet[wf_position] - temperature_inlet[sf_position]
if Delta_T1 == Delta_T2:
    LMTD = Delta_T1
else:
    LMTD = (Delta_T1 - Delta_T2) / np.log(Delta_T1/Delta_T2)

Ft = Ft_calc.Ft_Air_cooler(temperature_inlet, temperature_outlet, wf_position, sf_position, LMTD, np_tubi)

# Calcolo Area di Scambio
A = Q / (Uo * LMTD * Ft) # [m^2]







