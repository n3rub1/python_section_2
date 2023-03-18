# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 08:51:39 2023

@author: User
"""

import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Ellipse

window = tk.Tk()
window.title("Confidence Ellipses Calculator - Prototype")
window.resizable(0,0)
"""
# , "calPlanck", "acib217", "xi", "asz143", "aps100", "aps143"
#                      , "aps143217", "aps217", "aksz", "kgal100", "kgal143", "kgal143217", "kgal217", "galfTE100", "galfTE100143", "galfTE100217"
#                      , "galfTE143", "galfTE143", "galfTE143217", "galfTE217", "cal0", "cal2"


# def print_selected_value():
#     # Get the selected value from the radio button
#     selected_value = radioButtonsVariableX.get()
#     print(f"Selected value: {selected_value}")
# submit_button = tk.Button(window, text="Print selected value", command=print_selected_value)
# submit_button.grid(row=3, column=0)

# button= tk.Button(window, text="New")
# # button.pack()

importButton = tk.Button(window, text="Import", command=partial(buttonCallBack, "Import"))

"""

def importButtonCallBack():
    """prints back the buttons text"""
    print("Import")
    
def clearButtonCallBack():
    """prints back the buttons text"""
    print("Clear")

def calculateButtonCallBack():
    """prints back the buttons text"""
    print("Calculate")

def exitButtonCallBack():
    """prints back the buttons text"""
    print("Exit")

def plot():
    mean = (0, 0)
    a = 2
    b = 1
    theta = 45

    ellipseOne = Ellipse(xy=mean, width=a, height=b, angle=180-theta, fill=False, color='blue')
    ellipseTwo = Ellipse(xy=mean, width=a+1, height=b+1, angle=180-theta, fill=False, color='blue', linestyle = "dotted")
    fig, ax = plt.subplots()
    ax.add_patch(ellipseOne)
    ax.add_artist(ellipseTwo)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, rowspan=4, columnspan=len(radioButtonValuesX))

radioButtonValuesX = ["omegabh2", "omegach2", "theta", "tau", "w", "logA", "ns"]
radioButtonValuesY = ["omegabh2", "omegach2", "theta", "tau", "w", "logA", "ns"]
radioButtonsVariableX = tk.StringVar()
radioButtonsVariableX.set("omegabh2")
radioButtonsVariableY = tk.StringVar()
radioButtonsVariableY.set("omegabh2")
index = 0

#for the X
for button in radioButtonValuesX:
    b = tk.Radiobutton(window, text = button, variable = radioButtonsVariableX, value = button).grid(row=4, column=index)
    index = index + 1

#make sure versus is always in the middle
tk.Label(window, text = "VS").grid(row=5, column=int(len(radioButtonValuesX)/2))

index = 0
#for the Y
for button in radioButtonValuesY:
    b = tk.Radiobutton(window, text = button, variable = radioButtonsVariableY, value = button).grid(row=6, column=index)
    index = index + 1

customLabelX = tk.Label(window, text = "Custom X value").grid(row=4, column=len(radioButtonValuesX) + 1)
customLabelY = tk.Label(window, text = "Custom Y value").grid(row=6, column=len(radioButtonValuesX) + 1)
customEntryX = tk.Entry(window).grid(row = 4, column=len(radioButtonValuesX)+ 2)
customEntryY = tk.Entry(window).grid(row = 6, column=len(radioButtonValuesX)+ 2)

importButton = tk.Button(window, text = "Import", command = importButtonCallBack).grid(row = 0, column=len(radioButtonValuesX) + 1, columnspan=2)
clearButton = tk.Button(window, text = "Clear", command = clearButtonCallBack).grid(row = 1, column=len(radioButtonValuesX) + 1, columnspan=2)
calculateButton = tk.Button(window, text = "Calculate", command = calculateButtonCallBack).grid(row = 2, column=len(radioButtonValuesX) + 1, columnspan=2)
exitButton = tk.Button(window, text = "Exit", command = exitButtonCallBack).grid(row = 3, column=len(radioButtonValuesX) + 1, columnspan=2)

plot()

window.mainloop()