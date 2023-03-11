# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 07:55:58 2023

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# variableNameDictionary = {}
variableNames = ["omegabh2", "omegach2", "theta", "tau", "w", "logA", "ns"]

def getFullDataMatrix():
    """This method is used to bring in the original matrix"""
    #get the data using the genfromtxt method, do not skip the header, and get the data based on the tab spaces
    return np.genfromtxt("base_w_plikHM_TTTEEE_lowl_lowE_BAO_Riess18_Pantheon.covmat",skip_header=0, names=True, delimiter="\t")
    
def _invertMatrix(matrix):
    """Inverts the data of the matrix"""
    return np.linalg.inv(matrix)

def _fisherMatrixAndPrint(matrix, variableNameOne, variableNameTwo):
    """transforms the two arrays into a fisher matrix and prints out the values"""
    horizontalLines = "-" * 50
    fisherMatrix = _invertMatrix(matrix)
    
    print("{0} VS {1}".format(variableNameOne, variableNameTwo))
    print(horizontalLines)
    print("Convariance Matrix")
    print(matrix)
    print(horizontalLines)
    print("Fisher Matrix")
    print(fisherMatrix, end="\n\n")
    return(fisherMatrix)

def getMatrixValues(variableNameOne, variableNameTwo, fullDataMatrix):
    """Gets the matrix values depeneding on the variable names passed"""
    indexOfVariableOne = variableNames.index(variableNameOne)
    indexOfVariableTwo = variableNames.index(variableNameTwo)
    
    x = fullDataMatrix[indexOfVariableOne][indexOfVariableOne]
    xy = fullDataMatrix[indexOfVariableTwo][indexOfVariableOne]
    y = fullDataMatrix[indexOfVariableTwo][indexOfVariableTwo]
    
    firstRow = np.array([x, xy])
    secondRow = np.array([xy, y])
    matrix = np.stack((firstRow, secondRow), axis = 0)
    result = _fisherMatrixAndPrint(matrix, variableNameOne, variableNameTwo)
    
    return result
    
def _calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Get the confidence ellipse parameters for a^2"""
    firstPart = ((sigmaX**2) + (sigmaY**2)) / 2
    innerPart = (((sigmaX**2 - sigmaY**2)/4)+(sigmaX**2 * sigmaY))
    root = np.sqrt(innerPart)
    return (firstPart + root) * levelOfConfidence

def _calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Get the confidence ellipse parameters for b^2"""
    firstPart = ((sigmaX**2) + (sigmaY**2)) / 2
    innerPart = (((sigmaX**2 - sigmaY**2)/4)+(sigmaX**2 * sigmaY))
    root = np.sqrt(innerPart)
    return (firstPart - root) * levelOfConfidence
    

def _calculateTanTwotheta(sigmaX, sigmaY, sigmaXY):
    """Get the confidence ellipse parameters for tan(20)"""
    return ((2*sigmaX*sigmaY)/(sigmaX**2-sigmaY**2))
    
def ellipse_params(fisherMatrix, levelOfConfidence):
    """Calculate the width, height and inclination depending of the level of confidence"""
    #DO A PRINT!!
    return np.array([_calculateASquared(fisherMatrix[0][0], fisherMatrix[1][1], fisherMatrix[0][1], levelOfConfidence)
                     , _calculateBSqaured(fisherMatrix[0][0], fisherMatrix[1][1], fisherMatrix[0][1], levelOfConfidence)
                     , _calculateTanTwotheta(fisherMatrix[0][0], fisherMatrix[1][1], fisherMatrix[0][1])])
    
    
#initilalize the program and fill up the matrix values in the dictionary
fullDataMatrix = getFullDataMatrix()

#calculate and print the covariance and fisher matrices
omegabh2VsOmegach2 = getMatrixValues("omegabh2", "omegach2", fullDataMatrix)
omegach2VsW = getMatrixValues("omegach2", "w", fullDataMatrix)
logAVsNs = getMatrixValues("logA", "ns", fullDataMatrix)
tauVsW = getMatrixValues("tau", "w", fullDataMatrix)

ellipseOmegabh2VsOmegach2_1 = ellipse_params(omegabh2VsOmegach2, 1.52)
ellipseOmegach2VsW_1 = ellipse_params(omegach2VsW, 1.52)
ellipseLogAVsNs_1 = ellipse_params(logAVsNs, 1.52)
ellipseTauVsW_1 = ellipse_params(tauVsW, 1.52)

ellipseOmegabh2VsOmegach2_2 = ellipse_params(omegabh2VsOmegach2, 2.48)
ellipseOmegach2VsW_2 = ellipse_params(omegach2VsW, 2.48)
ellipseLogAVsNs_2 = ellipse_params(logAVsNs, 2.48)
ellipseTauVsW_2 = ellipse_params(tauVsW, 2.48)

fig, ax = plt.subplots()
ellipse = Ellipse(xy=(0, 0), width=ellipseOmegabh2VsOmegach2_1[0], height=ellipseOmegabh2VsOmegach2_1[1], angle=ellipseOmegabh2VsOmegach2_1[2])
ax.add_artist(ellipse)
plt.show()

















