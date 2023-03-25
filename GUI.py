# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 08:51:39 2023

@author: Robert SM4391
"""

import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Ellipse
import numpy as np

window = tk.Tk()
window.title("Confidence Ellipses Calculator - Prototype")
window.resizable(0,0)

def importMatrixButton():
    """prints back the buttons text"""
    print("Import Matrix")
    
def saveButton():
    """prints back the buttons text"""
    print("Save")
    
def importVariablesButton():
    """prints back the buttons text"""
    print("Import Variables")

def plotButton():
    """prints back the buttons text"""
    print("Plot")

def exitButton():
    """prints back the buttons text"""
    print("Exit")

def plot():
    mean = (0, 0)
    a = 2
    b = 1
    confidenceOne = 1
    confidenceTwo = 2
    theta = 45

    ellipseOne = Ellipse(xy=mean, width=a * confidenceOne, height=b * confidenceOne, angle=theta, fill=False, color='blue')
    ellipseTwo = Ellipse(xy=mean, width=a * confidenceTwo, height=b * confidenceTwo, angle=theta, fill=False, color='blue', linestyle = "dotted")
    
    cov = np.array([[0.2,0.7],[0.7,0.3]])
    points = np.random.multivariate_normal(mean, cov , size=500)
    
    fig, ax = plt.subplots()
    ax.add_patch(ellipseOne)
    ax.add_artist(ellipseTwo)
    plt.plot(points[:,0], points[:,1], ".", color='black', alpha=0.3)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
       
    sigmaNameOne = "σ-1"
    sigmaOneBox = plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="blue")
    sigmaNameTwo = "σ-2"
    sigmaTwoBox = plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="blue", linestyle="--")
    plt.legend(handles=[sigmaOneBox, sigmaTwoBox], labels=[sigmaNameOne, sigmaNameTwo], framealpha=0.8, loc="upper right")
    
    ax.set_xlabel("omeghach2")
    ax.set_ylabel("omeghabh2")
    
    ax.plot(0, 0, "ro")
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, padx= 10, pady=10)


customLabelX = tk.Label(window, text = "Variable X").grid(row=5, column = 0, pady=10)
customEntryX = tk.Entry(window).grid(row = 5, column = 1, columnspan = 2, pady=10)

customLabelY = tk.Label(window, text = "Variable Y").grid(row=6, column = 0, pady=10)
customEntryY = tk.Entry(window).grid(row = 6, column = 1, columnspan = 2, pady=10)

sigmaOneLabel = tk.Label(window, text = "Center X").grid(row=7, column = 0, pady=10)
sigmaOneEntry = tk.Entry(window).grid(row = 7, column = 1, columnspan = 2, pady=10)

sigmaTwoLabel = tk.Label(window, text = "Center Y").grid(row=5, column = 4, pady=10)
sigmaTwoEntry = tk.Entry(window).grid(row = 5, column = 5, pady=10)

centreCValue = tk.Label(window, text = "σ-1").grid(row=6, column = 4, pady=10)
centreXEntry = tk.Entry(window).grid(row = 6, column = 5, pady=10)

centreYValue = tk.Label(window, text = "σ-2").grid(row=7, column = 4, pady=10)
centreYEntry = tk.Entry(window).grid(row = 7, column = 5, pady=10)

importMatrix = tk.Button(window, text = "Import Matrix", command = importMatrixButton, width = 20).grid(row = 0, column=6, padx=10)
importVariables = tk.Button(window, text = "Import Variables", command = importVariablesButton, width = 20).grid(row = 1, column=6, padx=10)
save = tk.Button(window, text = "Save", command = saveButton, width = 20).grid(row = 2, column=6, padx=10)
plottingButton = tk.Button(window, text = "Plot", command = plotButton, width = 20).grid(row = 3, column=6, padx=10)
exitProgramButton = tk.Button(window, text = "Exit", command = exitButton, width = 20).grid(row = 4, column=6, padx=10)

plot()

window.mainloop()