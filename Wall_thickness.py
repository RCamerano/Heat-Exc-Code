# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 09:25:17 2021

@author: Riccardo
"""

def wall_thickness(Do,wt_min):
    if Do == 0.5:
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
    return wt