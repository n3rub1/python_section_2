# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 09:43:53 2023

@author: User
"""

import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Ellipse
import numpy as np

window = tk.Tk()
window.title("Error Ellipses - Prototype")
window.resizable(0,0)

def getHeaderVariables():
    """get the list of all the variable names from the file"""
    with open("base_w_plikHM_TTTEEE_lowl_lowE_BAO_Riess18_Pantheon.covmat", "r") as file:
        headerLine = file.readline().strip()
        variableList = headerLine.split("\t") 
    return variableList


def getFullDataMatrix():
    """This method is used to bring in the original matrix"""
    #get the data using the genfromtxt method, do not skip the header, and get the data based on the tab spaces    
    return np.genfromtxt("base_w_plikHM_TTTEEE_lowl_lowE_BAO_Riess18_Pantheon.covmat",skip_header=0, names=True, delimiter="\t")

    
def printCovarianceAndFisherValues(covarianceMatrix, fisherMatrix, variableNameOne, variableNameTwo):
    """Prints the covariance and fisher matrix values"""
    horizontalLines = "-" * 50
    print("{0} VS {1}".format(variableNameOne, variableNameTwo))
    print(horizontalLines)
    print("Convariance Matrix")
    print(covarianceMatrix)
    print(horizontalLines)
    print("Fisher Matrix")
    print(fisherMatrix, end="\n\n")
    

def getFisherMatrix(variableNameOne, variableNameTwo, covarianceMatrix):
    """transforms the two arrays into a fisher matrix and prints out the values"""
    fisherMatrix = np.linalg.inv(covarianceMatrix)
    return(fisherMatrix)


def getCovarianceMatrix(variableNameOne, variableNameTwo, fullDataMatrix):
    """Gets the covariance matrix values depeneding on the variable names passed"""
    indexOfVariableOne = variableNames.index(variableNameOne)
    indexOfVariableTwo = variableNames.index(variableNameTwo)
    
    x = fullDataMatrix[indexOfVariableOne][indexOfVariableOne]
    xy = fullDataMatrix[indexOfVariableTwo][indexOfVariableOne]
    y = fullDataMatrix[indexOfVariableTwo][indexOfVariableTwo]
    
    firstRow = np.array([x, xy])
    secondRow = np.array([xy, y])
    matrix = np.stack((firstRow, secondRow), axis = 0)
    return matrix


def _calculateInnerFormula(sigmaX, sigmaY, sigmaXY):
    """calculate the inner parts of the ellipse parameters"""
    firstPart = ((sigmaX + sigmaY) / 2)
    innerPart = ((((sigmaX - sigmaY)**2)/4)+(sigmaXY**2))
    root = np.sqrt(innerPart)
    return [firstPart, root]

    
def _calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Get the confidence ellipse parameters for a^2"""
    aParts = _calculateInnerFormula(sigmaX, sigmaY, sigmaXY)
    aSquared = aParts[0] + aParts[1]
    a = np.sqrt(aSquared)
    a = a * levelOfConfidence
    return a


def _calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Get the confidence ellipse parameters for b^2"""
    bParts = _calculateInnerFormula(sigmaX, sigmaY, sigmaXY)
    bSquared = bParts[0] - bParts[1]
    b = np.sqrt(bSquared)
    b = b * levelOfConfidence
    return b

    
def _calculateTanTwoTheta(sigmaX, sigmaY, sigmaXY):
    """Get the confidence ellipse parameters for tan(20)"""
    tanTwoTheta = ((2*sigmaXY)/(sigmaX-sigmaY))
    twoTheta = np.arctan(tanTwoTheta)
    theta = twoTheta / 2
    thetaInDegrees = (theta * 360) / (2 * np.pi)
    return thetaInDegrees


def printEllipse_params(params, name):
    """prints the width, height and inclination for all ellipses with α = 1.52"""
    print("{0}  at α = 1.52 \n height (a): {1}, width(b): {2}, inclination in degrees(θ): {3}\n".format(name, params[0], params[1], params[2]))
    
def ellipse_params(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Calculate the width, height and inclination depending of the level of confidence"""
    a = _calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence)
    b = _calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence)
    thetaInDegrees = _calculateTanTwoTheta(sigmaX, sigmaY, sigmaXY)
    result = np.array([b, a, thetaInDegrees])
    return result
    
def importVariablesButton():
    """prints back the buttons text"""
    print("Import Variables")
    
    
def saveButton():
    """prints back the buttons text"""
    print("Save")
    

def plotButton():
    """prints back the buttons text"""
    plot(ellipseResultFromOmegabh2VsOmegach2_A[0], ellipseResultFromOmegabh2VsOmegach2_B[0],
         ellipseResultFromOmegabh2VsOmegach2_A[1], ellipseResultFromOmegabh2VsOmegach2_B[1],
         ellipseResultFromOmegabh2VsOmegach2_A[2], ellipseResultFromOmegabh2VsOmegach2_B[2],
         0.022, 0.12, omegabh2VsOmegach2CovarianceMatrix, "Ωbh2", "Ωch2")

def exitButton():
    """prints back the buttons text"""
    print("Exit")
    
def makeBlankCanvas():
    """makes a blank canvas at startup"""
    
    fig, ax = plt.subplots()
    
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, padx= 10, pady=10)

def plot(widthA, widthB, heightA, heightB, inclinationA, inclinationB, x, y, covarianceMatrix, variableNameX, variableNameY):
    
    mean = np.array([x, y])
    ellipseOneSigma = Ellipse(xy=(x, y), width=widthA, height=heightA, angle=inclinationA, fill=False, color='blue')
    ellipseTwoSigma = Ellipse(xy=(x, y), width=widthB, height=heightB, angle=inclinationB, fill=False, color='blue', linestyle = "dotted")

    # Generate random points using the provided covariance matrix
    points = np.random.multivariate_normal(mean, covarianceMatrix, size=500)

    # Create the plot and add the ellipse to it
    fig, ax = plt.subplots()
    ax.add_patch(ellipseOneSigma)
    ax.add_artist(ellipseTwoSigma)

    # Plot the random points
    ax.scatter(points[:, 0], points[:, 1], alpha=1, c = "black", s=1)

    ax.plot(x, y, "ro", markersize = 2)
    ax.set_xlabel(variableNameX)
    ax.set_ylabel(variableNameY)
    
    sigmaNameOne = "σ-1"
    sigmaOneBox = plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="blue")
    sigmaNameTwo = "σ-2"
    sigmaTwoBox = plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="blue", linestyle="--")
    plt.legend(handles=[sigmaOneBox, sigmaTwoBox], labels=[sigmaNameOne, sigmaNameTwo], framealpha=0.8, loc="upper right")
    
    ax.set_xlim(xmin=min(x-widthA, x-widthB), xmax=max(x+widthA, x+widthB))
    ax.set_ylim(ymin=min(y-heightA, y-heightB), ymax=max(y+heightA, y+heightB))
    
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, padx= 10, pady=10)
    

#initilalize the program and fill up the matrix values in the dictionary
variableNames = getHeaderVariables()
fullDataMatrix = getFullDataMatrix()

#get the converiance matrices of the requested values
omegabh2VsOmegach2CovarianceMatrix = getCovarianceMatrix("omegabh2", "omegach2", fullDataMatrix)
omegach2VsWCovarianceMatrix = getCovarianceMatrix("omegach2", "w", fullDataMatrix)
logAVsNsCovarianceMatrix = getCovarianceMatrix("logA", "ns", fullDataMatrix)
tauVsWCovarianceMatrix = getCovarianceMatrix("tau", "w", fullDataMatrix)

#get the fisher matrices of the requested values
omegabh2VsOmegach2FisherMatrix = getFisherMatrix("omegabh2", "omegach2", omegabh2VsOmegach2CovarianceMatrix)
omegach2VsWFisherMatrix = getFisherMatrix("omegach2", "w", omegach2VsWCovarianceMatrix)
logAVsNsFisherMatrix = getFisherMatrix("logA", "ns", logAVsNsCovarianceMatrix)
tauVsWFisherMatrix = getFisherMatrix("tau", "w", tauVsWCovarianceMatrix)

#print converiance matrix and fisher matrix
printCovarianceAndFisherValues(omegabh2VsOmegach2CovarianceMatrix, omegabh2VsOmegach2FisherMatrix, "omegabh2", "omegach2")
printCovarianceAndFisherValues(omegach2VsWCovarianceMatrix, omegach2VsWFisherMatrix, "omegach2", "w")
printCovarianceAndFisherValues(logAVsNsCovarianceMatrix, logAVsNsFisherMatrix, "logA", "ns")
printCovarianceAndFisherValues(tauVsWCovarianceMatrix, tauVsWFisherMatrix, "tau", "w")

CONFIDENCE_A = 1.52
CONFIDENCE_B = 2.48

ellipseResultFromOmegabh2VsOmegach2_A = ellipse_params(omegabh2VsOmegach2CovarianceMatrix[0][0],
                                                       omegabh2VsOmegach2CovarianceMatrix[1][1],
                                                       omegabh2VsOmegach2CovarianceMatrix[0][1], CONFIDENCE_A)

ellipseResultFromOmegabh2VsOmegach2_B = ellipse_params(omegabh2VsOmegach2CovarianceMatrix[0][0],
                                                       omegabh2VsOmegach2CovarianceMatrix[1][1],
                                                       omegabh2VsOmegach2CovarianceMatrix[0][1], CONFIDENCE_B)

ellipseResultFromOmegach2VsW_A = ellipse_params(omegach2VsWCovarianceMatrix[0][0],
                                                omegach2VsWCovarianceMatrix[1][1],
                                                omegach2VsWCovarianceMatrix[0][1], CONFIDENCE_A)

ellipseResultFromOmegach2VsW_B = ellipse_params(omegach2VsWCovarianceMatrix[0][0],
                                              omegach2VsWCovarianceMatrix[1][1],
                                              omegach2VsWCovarianceMatrix[0][1], CONFIDENCE_B)

ellipseResultFromlogAVsNs_A = ellipse_params(logAVsNsCovarianceMatrix[0][0],
                                             logAVsNsCovarianceMatrix[1][1],
                                             logAVsNsCovarianceMatrix[0][1], CONFIDENCE_A)

ellipseResultFromlogAVsNs_B = ellipse_params(logAVsNsCovarianceMatrix[0][0],
                                             logAVsNsCovarianceMatrix[1][1],
                                             logAVsNsCovarianceMatrix[0][1], CONFIDENCE_B)

ellipseResultFromtauVsW_A = ellipse_params(tauVsWCovarianceMatrix[0][0],
                                             tauVsWCovarianceMatrix[1][1],
                                             tauVsWCovarianceMatrix[0][1], CONFIDENCE_A)

ellipseResultFromtauVsW_B = ellipse_params(tauVsWCovarianceMatrix[0][0],
                                             tauVsWCovarianceMatrix[1][1],
                                             tauVsWCovarianceMatrix[0][1], CONFIDENCE_B)

printEllipse_params(ellipseResultFromOmegabh2VsOmegach2_A, "Ωbh2 vs Ωch2")
printEllipse_params(ellipseResultFromOmegach2VsW_A, "Ωch2 vs w")
printEllipse_params(ellipseResultFromlogAVsNs_A, "ln(A) vs ns")
printEllipse_params(ellipseResultFromtauVsW_A, "τ vs w")
   

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

importMatrix = tk.Button(window, text = "Import Matrix", command = getFullDataMatrix, width = 20).grid(row = 0, column=6, padx=10)
importVariables = tk.Button(window, text = "Import Variables", command = getHeaderVariables, width = 20).grid(row = 1, column=6, padx=10)
save = tk.Button(window, text = "Save", command = saveButton, width = 20).grid(row = 2, column=6, padx=10)
plottingButton = tk.Button(window, text = "Plot", command = plotButton, width = 20).grid(row = 3, column=6, padx=10)
exitProgramButton = tk.Button(window, text = "Exit", command = exitButton, width = 20).grid(row = 4, column=6, padx=10)

makeBlankCanvas()

window.mainloop()