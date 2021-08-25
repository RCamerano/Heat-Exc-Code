# -*- coding: utf-8 -*-
"""
Created on Wed May 26 10:40:55 2021

@author: Paolo
"""

import numpy as np
import matplotlib as plt

def HE_plot(CVs):
    length_HE = np.linspace(0,len(CVs)-1,len(CVs))
    plot_variable = input('choose variable to plot:\n')
    a = []
    for control_volume in CVs:
        a.append(getattr(control_volume,plot_variable))
    plt.pyplot.plot(length_HE,a)