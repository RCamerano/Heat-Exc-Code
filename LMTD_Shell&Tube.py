# HEAT EXCHANGER MAIN SCRIPT
# CODICE PER DIMENSIONAMENTO SCAMBIATORI SHELL&TUBE CON METODO (F)LMTD

import Ft_calc
import TEMA_check
import numpy as np
from CoolProp.CoolProp import PropsSI

# Sez. 1 - DISPOSIZIONE FLUIDI
# All'utente è chiesto di specificare se il fluido di lavoro scorre nel tubo interno (i.e. 'pipe') o
# nel mantello esterno (i.e. 'shell'). Se viene utilizzato un input diverso il codice chiede di inserire nuovamente l'input.
# Verrà attribuito l'indice 'wf_position = 0' al fluido di lavoro se questo è all'interno ('pipe') o 'wf_position = 1' se è all'esterno ('shell').
# In tutto il resto del codice il primo elemento di una lista farà riferimento al fluido che scorre all'interno, il secondo elemento a 
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

# Sez. 3 - SCELTA MASSIME PERDITE DI CARICO ACCETTABILI
P_loss = [0,0]
P_loss[wf_position] = input("Specify maximum allowable pressure loss - WORKING FLUID [Pa]:\n")
P_loss[sf_position] = input("Specify maximum allowable pressure loss - SECONDARY FLUID [Pa]:\n")

# Sez. 4 - FLUIDO DI LAVORO
#
# INLET
#
# Inserimento delle condizioni del FLUIDO DI LAVORO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[wf_position] = 'Water'
pressure_inlet[wf_position] = 0 # [Pa]
temperature_inlet[wf_position] = 0 # [K]         # DA VERIFICARE QUESTO APPROCCIO!!!!
quality_inlet[wf_position]  = PropsSI('Q', 'T', temperature_inlet[wf_position], 'P', pressure_inlet[wf_position], fluid[wf_position])

# La portata del fluido di lavoro è un input del problema
M[wf_position] = 0 #[Kg/s]
# Se il flusso è sottoraffreddato o surriscaldato, si calcolano densità e entalpia sulla base di T e P,
# se il flusso è bifase, utilizziamo T e titolo di vapore.
if quality_inlet[wf_position] == -1.0:
    density_inlet[wf_position] = PropsSI('D', 'T', temperature_inlet[wf_position], 'P', pressure_inlet[wf_position], fluid[wf_position])
    enthalpy_inlet[wf_position] = PropsSI('H', 'T', temperature_inlet[wf_position], 'P', pressure_inlet[wf_position], fluid[wf_position])
else:  
    density_inlet[wf_position] = PropsSI('D', 'T', temperature_inlet[wf_position], 'Q', quality_inlet[wf_position], fluid[wf_position]) 
    enthalpy_inlet[wf_position] = PropsSI('H', 'T', temperature_inlet[wf_position], 'Q', quality_inlet[wf_position], fluid[wf_position])
#
# OUTLET
#
# Per il calcolo delle consizioni di uscita, si presuppone che il FLUIDO DI LAVORO subisca perdite di carico pari a quelle massime accettabili
pressure_outlet[wf_position] = pressure_inlet[wf_position] - P_loss[wf_position]   # CAPIRE COME DEFINIRE LO STATO DEL FLUSSO
temperature_outlet[wf_position] = 0 #[k]
quality_outlet[wf_position] = 0
if quality_outlet[wf_position] == -1.0:
    density_outlet[wf_position] = PropsSI('D', 'T', temperature_outlet[wf_position], 'P', pressure_outlet[wf_position], fluid[wf_position]) 
    enthalpy_outlet[wf_position] = PropsSI('H', 'T', temperature_outlet[wf_position], 'P', pressure_outlet[wf_position], fluid[wf_position])
else:
    density_outlet[wf_position] = PropsSI('D', 'T', temperature_outlet[wf_position], 'Q', quality_outlet[wf_position], fluid[wf_position]) 
    enthalpy_outlet[wf_position] = PropsSI('H', 'T', temperature_outlet[wf_position], 'Q', quality_outlet[wf_position], fluid[wf_position])


# Sez. 5 - FLUIDO SECONDARIO
#
# INLET
#
# Inserimento delle condizioni del FLUIDO SECONDARIO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[sf_position] = 'Water'
pressure_inlet[sf_position] = 0 #[Pa]          # DA VERIFICARE QUESTO APPROCCIO!!!!
temperature_inlet[sf_position] = 0 #[K]
quality_inlet[sf_position] = PropsSI('Q', 'T', temperature_inlet[sf_position], 'P', pressure_inlet[sf_position], fluid[sf_position])
# Se il flusso è sottoraffreddato o surriscaldato si calcolano densità ed entlpia sulla base di T e P,
# se il flusso è bifase utilizziamo T e titolo di vapore.
if quality_inlet[sf_position] == -1.0:
    density_inlet[sf_position] = PropsSI('D', 'T', temperature_inlet[sf_position], 'P', pressure_inlet[sf_position], fluid[sf_position]) 
    enthalpy_inlet[sf_position] = PropsSI('H', 'T', temperature_inlet[sf_position], 'P', pressure_inlet[sf_position], fluid[sf_position])
else:  
    density_inlet[sf_position] = PropsSI('D', 'T', temperature_inlet[sf_position], 'Q', quality_inlet[sf_position], fluid[sf_position]) 
    enthalpy_inlet[sf_position] = PropsSI('H', 'T', temperature_inlet[sf_position], 'Q', quality_inlet[sf_position], fluid[sf_position])
#
# OUTLET
#
# Per il calcolo delle consizioni di uscita, si presuppone che il FLUIDO SECONDARIO subisca perdite di carico pari a quelle massime accettabili
pressure_outlet[sf_position] = pressure_inlet[sf_position] - P_loss[sf_position]   # CAPIRE COME DEFINIRE LO STATO DEL FLUSSO
temperature_outlet[sf_position] = 0 #[k]
quality_outlet[sf_position] = 0
if quality_outlet[sf_position] == -1.0:
    density_outlet[sf_position] = PropsSI('D', 'T', temperature_outlet[sf_position], 'P', pressure_outlet[sf_position], fluid[sf_position]) 
    enthalpy_outlet[sf_position] = PropsSI('H', 'T', temperature_outlet[sf_position], 'P', pressure_outlet[sf_position], fluid[sf_position])
else:
    density_outlet[sf_position] = PropsSI('D', 'T', temperature_outlet[sf_position], 'Q', quality_outlet[sf_position], fluid[sf_position]) 
    enthalpy_outlet[sf_position] = PropsSI('H', 'T', temperature_outlet[sf_position], 'Q', quality_outlet[sf_position], fluid[sf_position])

# DEFINIZIONE MASSA FLUSSO SECONDARIO
# La portata massica del fluido secondario viene calcolata sulla base delle condizioni al contorno specificate per i due fluidi
# in modo tale che il bilancio energetico sia sempre rispettato.
M[sf_position] = abs(M[wf_position] * (enthalpy_inlet[wf_position] - enthalpy_outlet[wf_position]) / (enthalpy_inlet[sf_position] - enthalpy_outlet[sf_position]))


# Sez. 6 - SCELTE GENERALI DI DESIGN
condition = 'false'
# Selezionare categoria TEMA dello scambiatore (HE_cat)
HE_cat = input("Specify heat exchanger TEMA designation [E,F,X,K]:\n")
while condition == 'false':
    # Specificare numero passaggi lato tubi (np_tubi) e numero passaggi lato mantello (np_shell)
    np_tubi = input("Tube Passes number: \n")
    np_shell = input("Shell Passes number: \n")
    condition = TEMA_check.TEMA_check(HE_cat,np_shell)
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
    
# Specificare il valore del coefficiente globale di scambio di primo tentativo (Uo)
Uo = input("Estimate overall heat exchange coefficient [Uo]



