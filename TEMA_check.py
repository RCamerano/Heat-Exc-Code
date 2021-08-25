# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:13:57 2021

@author: Riccardo
"""

def TEMA_check(HE_cat,np_shell,np_tubi,):
    # TEMA E --> si accetta solo numero dispari di passaggi lato mantello. No limiti lato tubi
    if HE_cat == 'E' and np_shell%2 != 0:
        condition = 'true'
   
    # TEMA F --> caratterizzato da un lungo baffer longitudinale. In teoria è accettabile qualunque numero di passaggi
    #            lato mantello e tubi, ma solo il 2:2 è rilevante per i nostri fini. Solo questo sarà gestito dal codice 
    elif HE_cat == 'F' and np_shell == 2 and np_tubi == 2:
        condition = 'true'
    elif HE_cat == 'F' and (np_shell != 2 or np_tubi != 2):
        print("l'unico scambiatore TEMA F gestibile per ora è il 2:2. Correggere numero passaggi lato mantello e/o tubi")
    
    # TEMA K --> il kettle può essere utilizzato solo come evaporatore. Inoltre non è preciso parlare di numero di passaggi
    #            Per uniformare il codice, il kettle deve essere associato alla configurazione 1:2
    elif HE_cat == 'K'and (np_shell == 1 and np_tubi == 2):
        condition = 'true'
    elif HE_cat == 'K'and (np_shell != 1 or np_tubi != 2):
        print("l'evaporatore Kettle è associato a una configurazione 1:2. Correggere numero passaggi lato mantello e/o tubi")    
    
    
    elif HE_cat == 'X':
        condition = 'true'
    else:
        condition = 'false'
        print('numero di passaggi lato mantello non coerente con scambiatore scelto')
    return condition

def Ft_check(HE_cat,np_tubi):
    if HE_cat == 'E' and (np_tubi%2 != 0 and np_tubi != 1):
        condition = 'false'
        print('non esiste una formula per calcolare il valore di Ft. Scambiatore non usuale')
    return condition
        