# HEAT EXCHANGER MAIN SCRIPT
# CDrICE PER DIMENSIONAMENTO SCAMBIATORI SHELL&TUBE CON METDrO (F)LMTD

import math
import Ft_calc
import numpy as np
import pandas as pd
import FlowBC
import Wall_thickness
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
data = pd.read_excel(r'C:\Users\Utente1\Drcuments\Tifeo\Python\HE\LMTD\input.xlsx','Input_AC')

# Sez. 1 - DISPOSIZIONE FLUIDI
# Verrà attribuito l'indice 'wf_position = 0' al fluiDr di lavoro dato che questo è all'interno ('pipe').
# In tutto il resto del codice il primo elemento di una lista farà riferimento al fluiDr che scorre all'interno, il seconDr elemento a 
# quello che scorre all'esterno, indipendentemente dal fatto che il fluiDr di lavora sia nel pipe o nello shell.
# In tutto il codice le proprietà relative al fluiDr di lavoro sono caratterizate dall'indice 'wf_position'
# In tutto il codice le proprietà relative al fluiDr secondario sono caratterizate dall'indice 'sf_position'
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
# Ogni variabile è contenuta all'interno di una lista. Il primo elemento delle liste farà riferimento al fluiDr che scorre nel tubo interno,
# il seconDr all'elemento a quello che scorre nel mantello esterno.
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


# Sez. 4 - FLUIDr DI LAVORO
#
# Inserimento delle condizioni del FLUIDr DI LAVORO all'ingresso dello scambiatore. Se si tratta di flusso
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
# Inserimento delle condizioni del FLUIDr SECONDARIO all'ingresso dello scambiatore. Se si tratta di flusso
# monofase, specificare Pressione e Temperatura, poi calcolare titolo. Se si tratta di flusso bifase, 
# specificare Pressione O Temperatura (calcolare l'altro) e Titolo.
fluid[sf_position] = 'Air'
temperature_inlet[sf_position] = data.Selected_Value[9]
pressure_inlet[sf_position] = data.Selected_Value[8]
viscosity_inlet[sf_position] = PropsSI('V', 'T', temperature_inlet[sf_position], 'P', pressure_inlet[sf_position], fluid)
pressure_outlet[sf_position] = 101325


# Sez. 6 - SCELTE GENERALI DI DESIGN
#
# Viene specificato il numero passaggi lato tubi (np_tubi) e numero passaggi lato mantello (np_shell)
np_tubi = data.Selected_Value[19]
#
draft = data.Selected_Value[31]


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
viscosity_outlet[sf_position] = PropsSI('V', 'T', temperature_outlet[sf_position], 'P', pressure_outlet[sf_position], fluid)
enthalpy_mean[sf_position] = np.mean([enthalpy_inlet[sf_position],enthalpy_outlet[sf_position]])
density_mean[sf_position] = np.mean([density_inlet[sf_position],density_outlet[sf_position]])
viscosity_mean = np.mean([viscosity_inlet[sf_position],viscosity_outlet[sf_position]])

# DEFINIZIONE MASSA FLUSSO SECONDARIO
# La portata massica del fluiDr secondario viene calcolata sulla base delle condizioni al contorno specificate per i due fluidi
# in mDro tale che il bilancio energetico sia sempre rispettato.
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


# Sez. 8 - SCELTE TUBI
#
# Scelta lunghezza tubi (L) 
L = data.Selected_Value[23] # [m]
# Altezza fins
h_fins = data.Selected_Value[27] # [m]
# Numero di fins per metro
Nf = data.Selected_Value[26]
# Numero di file di tubi
Nr = data.Selected_Value[24]
# Tube arrangement
tube_disposition = data.Selected_Value[21]
# Fins thickness
t_fins = data.Selected_Value[28]
# Conduttività fins
k_fins = data.Selected_Value[32]

# Scelta diametro nominale tubi (Dr, root diameter)
condition = 'false'
while condition == 'false':
    Dr = data.Selected_Value[20] # [inch]
    if Dr == 0.25 or Dr == 0.5 or Dr == 0.75 or Dr == 1:
        condition = 'true'
    else:
        print('diametro nominale tubi non standard. Scegliere uno dei valori suggeriti')

# Calcolo minimo spessore tubi (wt_min) seconDr norma ASME B31.3 --> wt_min = (P.OD) / 2(SEW + PY) + CA
# P = massimo gradiente di pressione che si può verificare nello scambiatore al design point. Eventuali condizioni
# più severe sono considerate sceglienDr un opportuno fattore di sovradimensionamento.
# Un fattore di sicurezza 'SF' è utilizzato come margine di sicurezza per la progettazione.
P = np.array([pressure_inlet[0],pressure_outlet[0],100000])
Dr = Dr*0.0254 # [m]
#
# Calcolo del tube pitch (distanza tra il centro dei tubi) e scelta sulla disposizione dei tubi.
# Comunemente si ha pitch pari a 1.25, 1.33 o 1.5 diametro esterno del tubo. Questo è un valore preso come
# input da foglio Excel.
tube_disposition = data.Selected_Value[21]
pitch = data.Selected_Value[22]*Dr # [m]
    
# S = Maximum allowable stress Selected_Value for the selected material [Pa]. Riferirsi a tabella 1.A ASME B36.10 per carbon steel e ASME B36.19 per stainless steel
# Alcuni valori sono riportati qui: https://www.cis-inspector.com/asme-code-calculation-allowable-stresses.html
S = 170E6 # [Pa] 0.9*Yield strength of Carbon steel at 150°C

E = 204E9  # [Pa] Modulo di elasticità

# La temperatura della parete del tubo è approssimata come la media tra le temperature medie dei fluidi nei
# tubi e nel mantello.
T_wall = np.mean([np.mean([temperature_inlet[0],temperature_outlet[0]]),np.mean([temperature_inlet[1],temperature_outlet[1]])]) # [K]
k = 8.54 + 0.02*T_wall # Conducibilità del materiale # [W/m.K]

# E = Efficiency of ligaments
eta_joint = (pitch - Dr)/pitch # [-]

# CA = corrosion allowance
CA = data.Selected_Value[29] # [m]

# Sf = Safety Factor
Sf = data.Selected_Value[18] # [-]

# Dimensionamento dello spessore dei tubi - Procedure ASME
# 'P_ext' contiene le 3 pressioni esterne che il tubo potrebbe Drver sostenere durante le operazioni:
    # - Pressione del fluiDr in ingresso al mantello (critica nel caso in cui P_int < P_ext)
    # - Pressione del fluiDr in uscita al mantello (critica nel caso in cui P_int > P_ext)
    # - Pressione atmosferica, nel caso in cui i tubi si scoprissero (mantello "vuoto")
P_ext = [101325]
# Nel calcolo dello spessore minimo la formula ASME prevede l'utilizzo del raggio interno. 
# Per ragioni di procedura utiliziamo il raggio nominale del tubo# (sovradimensionamento dello spessore minimo).
wt_min[0] = Wall_thickness.minimum_wall_thickness(P,P_ext,eta_joint,Dr/2,S,L,T_wall,E,data) * Sf + CA # [m]

# Selezione schedula tubi e rispettivo spessore (wt). Calcolo diametro interno (ID) e diametro equivalente.
wt[0] = Wall_thickness.wall_thickness(Dr,wt_min[0]) #[m]
ID = Dr - 2*wt[0] # [m]
D_eq[0] = ID/np_tubi
OD = Dr + 2*h_fins

# Superficie di tubo liscio per metro di tubatura [m2/m]
A_root = np.pi*Dr
# Approssimazione della superficie si 1 (una) aletta [m2]
A_f1 = np.pi*(OD**2 - Dr**2)/2
# Area superficiale delle alette per metro di tubo [m2/m]
A_f2 = Nf*A_f1
# Superficie degli apici delle alette per metro di tubo [m2/m]
A_f3 = np.pi*Nf*h_fins*OD/2
# Area della tubatura non occupata da alette x m di tubo [m2/m]
A_r1 = (1 - Nf*t_fins)*A_root
# Area per meter [m2/m]
APM = A_f3 + A_f2 + A_r1
# Area per square meter per tube [m2/m2]
APMFPR = APM/pitch
# Area ratio
AR = APM/A_root
# Mean outside area
A_mean = np.mean([np.pi*Dr,np.pi*ID])
# Outside area
A_out = APM + A_root*(1 - t_fins*Nf)
# Inside area
A_in = np.pi*ID


# Sez. 9 - CALCOLI TUBI
#
# Calcolo n° tubi (nt), arrotondato al successivo numero pari
nt = A / Nr / APMFPR
nt = math.ceil(nt / 2.) * 2

# Nota le geometria dei tubi si calcola l'area di attraversamento del fluiDr.
A_cross[0] = (np.pi*ID**2/4)*nt/np_tubi

# Calcolo velocità lato tubi in ingresso e uscita lato tubi #m/s
velocity_inlet[0] = M[0] / (density_inlet[0] * A_cross[0])
velocity_outlet[0] = M[0] / (density_outlet[0] * A_cross[0])
velocity_mean[0] = np.mean([velocity_inlet[0],velocity_outlet[0]])

# Calcolo numero di Reynolds in ingresso e uscita lato tubi #-
Re_inlet[0] = density_inlet[0] * velocity_inlet[0] * ID  / viscosity_inlet[0]
Re_outlet[0] = density_outlet[0] * velocity_outlet[0] * ID / viscosity_outlet[0]
Re_mean[0] = np.mean([Re_inlet[0],Re_outlet[0]])
if Re_inlet[0] < 10000:
    print("Attenzione: basso Reynolds all'ingresso dello scambiatore lato tubi!")
if Re_outlet[0] < 10000:
    print("Attenzione: basso Reynolds all'uscita dello scambiatore lato tubi!")

# Sez. 10 - CALCOLI LATO ARIA

# Area frontale della prima bancata libera
Ar_fa1 = (pitch - Dr)*L*(nt/Nr - 1)

# Area frontale della prima bancata inclusa tra le alette
S = (1 - Nf*t_fins)/Nf
if S < 0:
    print('The number of fins per unit meter and fins thickness are not compatible')
Ar_fa2 = 2*h_fins*S*Nf*L*nt/Nr
#
Air_far = Ar_fa1 + Ar_fa2
# Air velocity
velocity_mean[1] = M[1]/density_mean[1]*Air_far


# Sez. 11 - SPECIFICA COEFFICIENTI DI SPORCAMENTO
#
h_fouling[0] = data.Selected_Value[15]
h_fouling[1] = data.Selected_Value[16]


# Sez. 12 - CALCOLO COEFFICIENTE DI SCAMBIO LATO ARIA
#
# Calcolo Reynolds
Re_mean[1] = Dr*density_mean[1]* velocity_mean[1]/viscosity_mean[1]

# Calcolo alfa
if draft == 'forced':
    cp = PropsSI('Cpmass', 'P', pressure_mean[1], 'T', temperature_mean[1], fluid[1])
    k_air = PropsSI('L', 'P', pressure_mean[1], 'T', temperature_mean[1], fluid[1])
    Pr_air = cp*k_air/viscosity_mean[1]
    
    Nu = 0.134*Re_mean[1]**0.681*Pr_air**0.33*(S/h_fins)**0.2*(S/t_fins)**0.1143*(1 + 8.467E-5*velocity_mean[1]/11811*Nr**2)**-0.14
    
    alfa[1] = Nu*Dr/k_air
    
elif draft == 'induced':
    cp = PropsSI('Cpmass', 'P', pressure_mean[1], 'T', temperature_mean[1], fluid[1])
    k_air = PropsSI('L', 'P', pressure_mean[1], 'T', temperature_mean[1], fluid[1])
    Pr_air = cp*k_air/viscosity_mean[1]   
    
    Nu = 0.287*Re_mean**0.685*Pr_air**0.33*AR**-0.31*(Nr/6)**0.138
    
    alfa[1] = Nu*Dr/k_air

# Calcolo resistenza termica delle alette
m = h_fins*np.sqrt(2)/((1/alfa[1 + h_fouling[1]])*k_fins*t_fins)
n_a = 1/(1 + m**2/3*np.sqrt(OD/Dr))
h_fin = (1 - n_a)*(1/alfa[1] + h_fouling[1])/(A_root/APM + n_a)

    
# Sez. 13 - CALCOLO COEFFICIENTE DI SCAMBIO LATO TUBI

# Calcolo alfa lato tubo. Abbiamo 3 procedure differenti in base al fenomeno di scambio termica che si presenta.
# - Raffreddamento
if quality_inlet[0] == quality_outlet[0]:
    alfa[0] = HT.heat_transfer_coefficient_1ph(temperature_mean,pressure_mean,fluid,Re_mean,j_h,D_eq[0],0,L)
# - Condensazione
elif quality_inlet[0] > quality_outlet[0]:
    alfa[0] = HT.heat_exchange_coefficient_condensation(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,quality_inlet,quality_outlet,D_eq,fluid,0,M)


# Sez. 14 - CALCOLO COEFFICIENTE GLOBALE DI SCAMBIO

Uo_star = 1/( 1/alfa[1] + h_fouling[1] + h_fin + wt*A_out/k/A_mean + (h_fouling[0] + 1/alfa[0])*A_out/A_mean )
Q_star = A_out * nt * Uo_star * LMTD * Ft


# Sez. 15 - CALCOLO PERDITE DI CARICO LATO FAN

# Friction factor lato aria
if tube_disposition == 't':
    f[1] = 18.93*Re_mean[1]**-0.316*(pitch/Dr)**-0.927
elif tube_disposition == 's':
    f[1] = 18.93*Re_mean[1]**-0.316*(pitch/Dr)**-0.927*0.7071**-0.52

# Le perdite di carico sono calcolate approssimando inizialmente il dimensionamento del fan usando la sola perdita
# di carica per attrito. Dopo un primo dimensionamento di massima del fan si calcolano le perdite di carico
# dovute all'effetto di ingresso e un secondo dimensionamento del fan è effettuato.

# Perdita di carico lato aria - First guess
delta_P[1] = f[1]*Nr*density_mean[1]*velocity_mean[1]**2/9.81
eff = data.Selected_value[0]
# Potenza richiesta dal fan per le solite perdite di attrito
BHP = delta_P[1]*M[1]/eff/density_mean[1]/201 # [hp]

# Dimensionamento fan.
n_fans = BHP/25
F_a = A_out / Nr / APMFPR
FAPF = 0.4 * F_a / n_fans
D_fan = np.sqrt(FAPF*np.pi/4)
AVFF = M[1] / n_fans / density_mean[1] / (np.pi/4) / (D_fan**2)

# Calcolo delle perdite per effetti di ingresso
delta_P[1] = delta_P[1] + density_mean[1]*AVFF**2/9.81

# Potenza complessiva richiesta dal fan
BHP = delta_P[1]*M[1]/eff/density_mean[1]/201 # [hp]


# Sez. 15 - CALCOLO PERDITE DI CARICO LATO TUBI

f[0] = friction_coefficient.friction_coeff_tube(Re_mean[0])
delta_P[0] = ( 4*f[0]*L*velocity_mean[0]**2 / D_eq[0]*np.mean([density_inlet[0],density_outlet[0]]) ) + velocity_mean[0]**2*density_mean[0]


