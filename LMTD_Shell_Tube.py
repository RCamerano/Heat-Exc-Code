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

# Gli input per il dimensionamento sono presi dal file excel 'input.xlsx'. Nel file sono presenti le
# seguenti colonne: 
    # - Category
    # - Parameter: qui è contenuto il nome del parametro di input
    # - Permissible Values: alcuni input non possono variare arbitrariamente, devono essere inseriti dei valori prestabiliti.
    # - Unit: unità di misura del parametro
    # - Selected Value: è il valore del parametro che viene caricato dal codice nel suo funzionamento
data = pd.read_excel(r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\input.xlsx')

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
    elif wf_position == 'shell':
        wf_position = 1
        sf_position = 0
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
fluid[sf_position] = 'Water'
[pressure_inlet[sf_position], temperature_inlet[sf_position], quality_inlet[sf_position], density_inlet[sf_position], enthalpy_inlet[sf_position], viscosity_inlet[sf_position], pressure_outlet[sf_position], temperature_outlet[sf_position], quality_outlet[sf_position], density_outlet[sf_position], enthalpy_outlet[sf_position], viscosity_outlet[sf_position]] = FlowBC.InletOutlet(fluid[sf_position], P_loss[sf_position], P_inlet = data.Selected_Value[8], T_inlet = data.Selected_Value[9], T_outlet = data.Selected_Value[10], Q_inlet = data.Selected_Value[11], Q_outlet = data.Selected_Value[12])
#
# Alcune relazioni presenti nel codice utilizzano il valore medio delle proprietà dei fluidi
temperature_mean[sf_position] = np.mean([temperature_inlet[sf_position],temperature_outlet[sf_position]])
pressure_mean[sf_position] = np.mean([pressure_inlet[sf_position],pressure_outlet[sf_position]])
enthalpy_mean[sf_position] = np.mean([enthalpy_inlet[sf_position],enthalpy_outlet[sf_position]])
density_mean[sf_position] = np.mean([density_inlet[sf_position],density_outlet[sf_position]])
quality_mean[sf_position] = np.mean([quality_inlet[sf_position],quality_outlet[sf_position]])

# DEFINIZIONE MASSA FLUSSO SECONDARIO
# La portata massica del fluido secondario viene calcolata sulla base delle condizioni al contorno specificate per i due fluidi
# in mDoo tale che il bilancio energetico sia sempre rispettato.
M[sf_position] = abs(M[wf_position] * (enthalpy_inlet[wf_position] - enthalpy_outlet[wf_position]) / (enthalpy_inlet[sf_position] - enthalpy_outlet[sf_position])) # [Kg/s]

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
        flow_mix[0] = data.Selected_Value[7]
        flow_mix[1] = data.Selected_Value[8]
        Ft = Ft_calc.Ft_TEMA_X(flow_mix,temperature_inlet[wf_position], temperature_outlet[wf_position], temperature_inlet[sf_position], temperature_outlet[sf_position]) 
       
# Specificare il valore del coefficiente globale di scambio di primo tentativo (Uo)
# Per stimare il valore di Uo vedi Coulson and Richardson's pag. 637 o https://www.engineersedge.com/thermodynamics/overall_heat_transfer-table.htm
Uo = data.Selected_Value[14] # [W/m^2.K]

# Calcolo Area di Scambio
A = Q / (Uo * LMTD * Ft) # [m^2]


# Sez. 8 - SCELTE TUBI
#
# Scelta lunghezza tubi (L) 
L = data.Selected_Value[27] # [m]

# Scelta diametro nominale tubi (Do)
condition = 'false'
while condition == 'false':
    Do = data.Selected_Value[23] # [inch]
    if Do == 0.25 or Do == 0.5 or Do == 0.75 or Do == 1:
        condition = 'true'
    else:
        print('diametro nominale tubi non standard. Scegliere uno dei valori suggeriti')

# Calcolo minimo spessore tubi (wt_min) secondo norma ASME B31.3 --> wt_min = (P.OD) / 2(SEW + PY) + CA
# P = massimo gradiente di pressione che si può verificare nello scambiatore al design point. Eventuali condizioni
# più severe sono considerate scegliendo un opportuno fattore di sovradimensionamento.
# Un fattore di sicurezza 'SF' è utilizzato come margine di sicurezza per la progettazione.
P = np.array([pressure_inlet[0],pressure_outlet[0],100000])
# OD = diametro esterno (outer diameter) secondo norme ASME
if Do == 0.25:
    OD = 0.54 # [inch]
elif Do == 0.5:
    OD = 0.84 # [inch]
elif Do == 0.75:
    OD = 1.05 # [inch]
elif Do == 1:
    OD = 1.315 # [inch]
else:
    print('errore nella scelta del diametro tubi')
#
# Conversione in metri del diametro esterno dei tubi.
OD = OD*0.0254 # [m]
Do = Do*0.0254 # [m]
#
# Calcolo del tube pitch (distanza tra il centro dei tubi) e scelta sulla disposizione dei tubi.
# Comunemente si ha pitch pari a 1.25, 1.33 o 1.5 diametro esterno del tubo. Questo è un valore preso come
# input da foglio Excel.
tube_disposition = data.Selected_Value[25]
pitch = data.Selected_Value[26]*OD # [m]
    
# S = Maximum allowable stress Selected_Value for the selected material [Pa]. Riferirsi a tabella 1.A ASME B36.10 per carbon steel e ASME B36.19 per stainless steel
# Alcuni valori sono riportati qui: https://www.cis-inspector.com/asme-code-calculation-allowable-stresses.html
S = 170E6 # [Pa] 0.9*Yield strength of Carbon steel at 150°C

E = 204E9  # [Pa] Modulo di elasticità

# La temperatura della parete del tubo è approssimata come la media tra le temperature medie dei fluidi nei
# tubi e nel mantello.
T_wall = np.mean([np.mean([temperature_inlet[0],temperature_outlet[0]]),np.mean([temperature_inlet[1],temperature_outlet[1]])]) # [K]
k = 8.54 + 0.02*T_wall # Conducibilità del materiale # [W/m.K]

# E = Efficiency of ligaments
eta_joint = (pitch - OD)/pitch # [-]

# CA = corrosion allowance
CA = data.Selected_Value[30] # [m]

# Sf = Safety Factor
Sf = data.Selected_Value[18] # [-]

# Dimensionamento dello spessore dei tubi - Procedure ASME
# 'P_ext' contiene le 3 pressioni esterne che il tubo potrebbe dover sostenere durante le operazioni:
    # - Pressione del fluido in ingresso al mantello (critica nel caso in cui P_int < P_ext)
    # - Pressione del fluido in uscita al mantello (critica nel caso in cui P_int > P_ext)
    # - Pressione atmosferica, nel caso in cui i tubi si scoprissero (mantello "vuoto")
P_ext = [pressure_inlet[1],pressure_outlet[1],100000]
# Nel calcolo dello spessore minimo la formula ASME prevede l'utilizzo del raggio interno. 
# Per ragioni di procedura utiliziamo il raggio nominale del tubo# (sovradimensionamento dello spessore minimo).
wt_min[0] = Wall_thickness.minimum_wall_thickness(P,P_ext,eta_joint,OD/2,S,L,T_wall,E,data) * Sf + CA # [m]

# Selezione schedula tubi e rispettivo spessore (wt). Calcolo diametro interno (ID) e diametro equivalente.
wt[0] = Wall_thickness.wall_thickness(Do,wt_min[0]) #[m]
ID = OD - 2*wt[0] # [m]
D_eq[0] = ID/np_tubi


# Sez. 9 - CALCOLI TUBI
#
# Calcolo n° tubi (nt), arrotondato al successivo numero pari
nt = A / (np.pi * OD * L)
nt = math.ceil(nt / 2.) * 2

# Nota le geometria dei tubi si calcola l'area di attraversamento del fluido.
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

# Sez. 10 - SCELTE MANTELLO
#
# Calcolo del diametro dello shell. 'K1' ed 'n' sono parametri tabulati.
K1,n = shell.shell_sizing(tube_disposition,np_tubi)
D_bundle = OD*(nt/K1)**(1/n)

# Calcolo del bundle clearance basato su approccio grafico.
# Bisogna inizialmente fissare la tipologia di scambiatore ('rear_head_type') scegliendo tra:
    # - Pull-through Floating head
    # - Split-ring Floating head
    # - Packed head
    # - U-tube
rear_head_type = data.Selected_Value[21]
shell_clearance = shell.clearance(rear_head_type,D_bundle)
D_shell_int = D_bundle + 2*shell_clearance # [m]

# Per diametri dello shell inferiori a 24", lo shell è ottenuto da un tubo commerciale.
# I tubi commerciali non hanno diametri che variano con continuità, ma sono disponibili solo determinate misure.
# Prendiamo queste misure da tabella importata ('tubi.xlsx') e selezioniamo il valore superiore a 
# 'D_shell' più prossimo rispetto a quello calcolato.
if D_shell_int < 24*0.0254: # conversione da inch a m
    # Dimensioni caratteristiche di tubi commerciali sono nella tabella 'tubi.xlsx', caricate poi nella
    # variabile 'data_tube' 
    data_tube = pd.read_excel(r'C:\Users\Utente1\Documents\Tifeo\Python\HE\LMTD\tubi.xlsx')

    # In questo passaggio si cerca nella tabella il diametro del tubo più prossimo al diametro interno
    # dello shell che è stato ottenuto. 
    # N.B.: la relazione che segue serve a decidere se selezioniamo il tubo con diametro inferiore o superiore
    # rispetto al diametro interno dello shell. Se all'interno della relazione è presente il segno '>' andremo a
    # selezionare un tubo con diametro INFERIORE rispetto al diametro calcolato. Viceversa se è presente il segno
    # '<' andremo a selezionare un tubo con diametro SUPERIORE rispetto al diametro calcolato.
    # Nel caso in cui viene selezionato un tubo con diametro inferiore rispetto al diametro interno dello shell
    # calcolato, la clearence precedentemente calcolata non verrà rispetta.
    data_tube['DN [mm]'][(data_tube['DN [mm]'] - D_shell_int) > 0] = 0
    D_shell_0 = data_tube.iloc[(data_tube['DN [mm]'] - D_shell_int).abs().argsort()[:1]]['DN [mm]'].tolist()[0]
    
    # Otteniamo la riga della tabela tubi cui corrisponde il diametro selezionato, tramite questo è possibile
    # risalire all altre caratteristiche del tubo selezionato.
    idx = max(data_tube.iloc[(data_tube['DN [mm]'] - D_shell_int).abs().argsort()[:1]].index.tolist())
    D_shell_ext = data_tube.D_ext[idx]
    
    # Dimensionamento dello spessore dei tubi - Procedure ASME
    # 'P' contiene le 3 pressioni interne del mantello che potrebbero dover essere sostenute durante le operazioni:
    # - Pressione del fluido in ingresso al mantello (critica nel caso in cui P_int > P_ext)
    # - Pressione del fluido in uscita al mantello (critica nel caso in cui P_int < P_ext)
    # - Pressione atmosferica, nel caso in cui il mantello fosse vuoto (mantello "vuoto")
    P_ext = [100000]
    P = [pressure_inlet[1],pressure_outlet[1],100000]
    wt_min[1] = Wall_thickness.minimum_wall_thickness(P,P_ext,eta_joint,D_shell_ext/2,S,L,temperature_mean,E,data) * Sf + CA # [m]
    
    # Una volta che il diametro dello shell è stato definito ed è stato calcolata lo spessore minimo che deve
    # avere lo shell per resistere all'effetto della pressione, entriamo nella tabella e selezioniamo tra gli spessori
    # disponibili quello più prossimo al minimo calcolato.
    wt_serie = data_tube.loc[idx,'XS':'sch.160'] # [mm]
    wt[1] = np.array(wt_serie)
    aux = wt[1] - wt_min[1] # [m]
    b = np.argwhere(aux > 0)
    wt[1] = wt_min[1] + float(min(aux[b])) # [m]
    # La schedula del tubo è ottenuta tramite l'indice dello spessore selezionato.
    Schedula_shell = max(wt_serie.iloc[(wt_serie - wt[1]).abs().argsort()[:1]].index.tolist())
    
    # Noto il diametro esterno dello shell e lo spessore del tubo da cui questo è ricavato possiamo calcolare
    # il diametro interno reale dello shell. Quando questo è stato definito si può calcolare il numero di tubi 
    # che possiamo inserire all'interno dello shell selezionato.
    # N.B.: il numero di tubi già calcolato viene sovrascritto con il nuovo numero di tubi
    D_shell_int = D_shell_ext - 2*wt[1]
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
    nt = math.ceil(nt / 2.) * 2
    
# Quando la geometria dello shell e dei tubi è fissata si calcola l'area di scambio risultante.
A = np.pi * OD * L * nt   
    

# Sez. 11 - CALCOLI MANTELLO
#    
# Viene caricato da tabella di input il rapporto tra la distanza fra i baffle e il diametro interno dello shell.
# Questo deve essere un valore compresso tra 0.2 ed 1.
decision = 'not accepted'
while decision == 'not accepted':
    B = data.Selected_Value[28]
    if B >= 0.2 and B <= 1:
        decision = 'accepted'
    else:
        print('the Selected_Value of the ratio has to be comprised between 0.2 and 1')
B = D_shell_int*B

# Il numero di baffle è basato sul loro distanziamento e la lunghezza totale dello scambiatore
n_baffle = math.ceil(L/B - 1)

# Calcolo dell'area di attraversamento lato mantello
A_cross[1] = (pitch - OD)*D_shell_int*B / pitch

# Calcolo della velocità di attraversamento lato mantello
velocity_inlet[1] = M[1] / (density_inlet[1]*A_cross[1])
velocity_outlet[1] = M[1] / (density_outlet[1]*A_cross[1])
velocity_mean[1] = np.mean([velocity_inlet[1],velocity_outlet[1]])

# Calcolo del diametro equivalente lato mantello.
if tube_disposition == 't':
    D_eq[1] = 1.1*(pitch**2 - 0.917*OD**2)/OD
else:
    D_eq[1] = 1.27*(pitch**2 - 0.785*OD**2)/OD

# Calcolo del numero di Reynolds lato mantello
Re_inlet[1] = density_inlet[1] * velocity_inlet[1] * D_eq[1]  / viscosity_inlet[1]
Re_outlet[1] = density_outlet[1] * velocity_outlet[1] * D_eq[1]  / viscosity_outlet[1]
Re_mean[1] = np.mean([Re_inlet[1],Re_outlet[1]])

# Il valore del baffle cut è caricato da tabella degli input.
baffle_cut = data.Selected_Value[29]


# Sez. 12 - VALUTAZIONE HEAT TRANSFER FACTOR
#
# Lato tubi
j_h[0] = HT.heat_transfer_factor_tube(np.mean([Re_inlet[0],Re_outlet[0]]),L/ID)
# Lato shell
j_h[1] = HT.heat_transfer_factor_shell(np.mean([Re_inlet[1],Re_outlet[1]]),baffle_cut)


# Sez. 13 - SPECIFICA COEFFICIENTI DI SPORCAMENTO
#
h_fouling[0] = data.Selected_Value[15]
h_fouling[1] = data.Selected_Value[16]


# Sez. 14 - VALUTAZIONE HEAT TRANSFER COEFFICIENT
#
# Calcolo alfa lato tubo. Abbiamo 3 procedure differenti in base al fenomeno di scambio termica che si presenta.
# - Raffreddamento/Riscaldamento
if quality_inlet[0] == quality_outlet[0]:
    alfa[0] = HT.heat_transfer_coefficient_1ph(temperature_mean,pressure_mean,fluid,Re_mean,j_h,D_eq[0],0,L)
# - Condensazione
elif quality_inlet[0] > quality_outlet[0]:
    alfa[0] = HT.heat_exchange_coefficient_condensation(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,quality_inlet,quality_outlet,D_eq,fluid,0,M)
# - Evaporazione
else:
    # Per il calcolo di alfa nel caso evaporativo è necessario conoscere il coefficiente di scambio del fluido
    # nell'altro lato della parete, in questo caso dobbiamo quindi calcolare alfa lato mantello nel caso di
    # raffreddamento/condensazione (non può esserci evaporazione da entrambi i lati della parete)
    if quality_inlet[1] == quality_outlet[1]:
        alfa[1] = HT.heat_transfer_coefficient_1ph(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,fluid,Re_inlet,Re_outlet,j_h,D_eq[1],1,L)
    elif quality_inlet[1] > quality_outlet[1]:
        alfa[1] = HT.heat_exchange_coefficient_condensation(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,quality_inlet,quality_outlet,D_eq[1],fluid,1,M)
    # Noto alfa dall'altra lato della parete possiamo calcolare alfa di evaporazione.
    alfa[0] = HT.heat_exchange_coefficient_evaporation(temperature_mean,pressure_mean,quality_mean,enthalpy_mean,density_mean,D_eq,velocity_mean,fluid,0,alfa[1],OD,ID,k,h_fouling,LMTD)
    alfa[0] = HT.heat_exchange_coefficient_evaporation_Liu(temperature_mean,pressure_mean,quality_mean,enthalpy_mean,density_mean,D_eq,velocity_mean,fluid,0,alfa[1],OD,ID,k,h_fouling,LMTD,Q/A,M,data,Ft)

    
# Calcolo alfa lato tubo. Abbiamo 3 procedure differenti in base al fenomeno di scambio termica che si presenta.
# - Raffreddamento/Riscaldamento
if quality_inlet[1] == quality_outlet[1] and alfa[1] == 0:
    alfa[1] = HT.heat_transfer_coefficient_1ph(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,fluid,Re_inlet,Re_outlet,j_h,D_eq[1],1,L)
# - Condensazione
elif quality_inlet[1] > quality_outlet[1] and alfa[1] == 0:
    alfa[1] = HT.heat_exchange_coefficient_condensation(temperature_inlet,temperature_outlet,pressure_inlet,pressure_outlet,quality_inlet,quality_outlet,D_eq[1],fluid,1,M)
# - Evaporazione
elif quality_inlet[1] < quality_outlet[1]:
    # In questo caso alfa lato tubi è già noto e può essere inserito direttamente tra gli input della funzione.
    alfa[1] = HT.heat_exchange_coefficient_evaporation(temperature_mean,pressure_mean,quality_mean,enthalpy_mean,density_mean,D_eq,velocity_mean,fluid,1,alfa[0],OD,ID,k,h_fouling,LMTD)
    alfa[1] = HT.heat_exchange_coefficient_evaporation_Liu(temperature_mean,pressure_mean,quality_mean,enthalpy_mean,density_mean,D_eq,velocity_mean,fluid,1,alfa[0],OD,ID,k,h_fouling,LMTD,Q/A,M,data,Ft)


# Sez. 15 - CALCOLO COEFFICIENTE GLOBALE DI SCAMBIO (REALE)
#
# Il valore del coefficiente di scambio globale viene calcolato per poi essere confrontao con il valore Uo di 
# primo tentativo.
Uo_star = (1/alfa[1] + 1/h_fouling[1] + OD*np.log(OD/ID)/(2*k) + (1/alfa[0])*(OD/ID) + (1/h_fouling[0])*(OD/ID))**-1
# Il calore realmente scambiato è calcolato tramite il valore del coefficiente di scambio reale appena calcolato.
Q_star = A * Uo_star * LMTD * Ft


# Sez. 16 -  VALUTAZIONE PERDITE DI CARICO
#
# Calcolo del friction coefficient tabellato
f[0] = friction_coefficient.friction_coeff_tube(np.mean([Re_inlet[0],Re_outlet[0]]))
f[1] = friction_coefficient.friction_coeff_shell(np.mean([Re_inlet[1],Re_outlet[1]]),baffle_cut)

# Calcolo delle perdite di carico.
delta_P[0] = ( 4*f[0]*L*velocity_mean[0]**2 / D_eq[0]*np.mean([density_inlet[0],density_outlet[0]]) ) + velocity_mean[0]**2*np.mean([density_inlet[0],density_outlet[0]])
delta_P[1] = ( 4*f[1]*L*velocity_mean[1]**2 / D_eq[1]*np.mean([density_inlet[1],density_outlet[1]]) )

# SEZ. 17 - DESIGN OUTPUT

print('Uo_star is:', float("{0:.2f}".format(Uo_star)))
print('Uo_star is', 100*float("{0:.2f}".format((Uo_star-Uo)/Uo)), '% with respect to Uo')
print('pressure loss in the tube side is', float("{0:.2f}".format(delta_P[0])), '[Pa]. Selected threshold Selected_Value is', P_loss[0])
print('pressure loss in the shell side is', float("{0:.2f}".format(delta_P[1])), '[Pa]. Selected threshold Selected_Value is', P_loss[1])
print('heat transfer coefficinet in the tube side is', float("{0:.2f}".format(alfa[0])), 'W/m^2 K')
print('heat transfer coefficinet in the shell side is', float("{0:.2f}".format(alfa[1])), 'W/m^2 K')