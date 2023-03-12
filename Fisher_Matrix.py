# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 07:55:58 2023

@author: Robert Gatt SM4391
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
    # result = _fisherMatrixAndPrint(matrix, variableNameOne, variableNameTwo)
    return matrix

def _calculateInnerFormula(sigmaX, sigmaY, sigmaXY):
    """calculate the inner parts of the ellipse parameters"""
    firstPart = ((sigmaX**2 + sigmaY**2) / 2)
    innerPart = ((((sigmaX**2 - sigmaY**2)**2)/4)+(sigmaXY**2))
    root = np.sqrt(innerPart)
    return [firstPart, root]
    
def _calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Get the confidence ellipse parameters for a^2"""
    a = _calculateInnerFormula(sigmaX, sigmaY, sigmaXY)
    return (a[0] + a[1]) * levelOfConfidence

def _calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Get the confidence ellipse parameters for b^2"""
    b = _calculateInnerFormula(sigmaX, sigmaY, sigmaXY)
    return (b[0] - b[1]) * levelOfConfidence
    
def _calculateTanTwoTheta(sigmaX, sigmaY, sigmaXY):
    """Get the confidence ellipse parameters for tan(20)"""
    return ((2*sigmaXY)/(sigmaX**2-sigmaY**2))
    
def ellipse_params(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Calculate the width, height and inclination depending of the level of confidence"""
    a = np.sqrt(np.abs(_calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence)))
    b = np.sqrt(np.abs(_calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence)))
    tan2Theta = 0.5 * np.arctan(_calculateTanTwoTheta(sigmaX, sigmaY, sigmaXY))
    result = np.array([a, b, tan2Theta])
    print("Width={0}, Height={1}, Inclination={2}".format(a, b, tan2Theta))
    return result
    
#initilalize the program and fill up the matrix values in the dictionary
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

#[0][0] sigmaX, [0][1] sigma xy [1][0] sigma xy [1][1] sigma y

# sigmaX = float(input("Enter sigmaX: "))
# sigmaY = float(input("Enter sigmaY: "))
# sigmaXY = float(input("Enter sigmaXY: "))
confidenceA = 1.52
# # confidenceB = 2.48

print("SigmaX {0}".format(omegabh2VsOmegach2FisherMatrix[0][0]))
print("SigmaY {0}".format(omegabh2VsOmegach2FisherMatrix[1][1]))
print("SigmaXY {0}".format(omegabh2VsOmegach2FisherMatrix[0][1]))


# ellipseResultFromUser = ellipse_params(sigmaX, sigmaY, sigmaXY, confidenceA)
ellipseResultFromOmegabh2VsOmegach2 = ellipse_params(omegabh2VsOmegach2CovarianceMatrix[0][0], omegabh2VsOmegach2CovarianceMatrix[1][1], omegabh2VsOmegach2CovarianceMatrix[0][1], confidenceA)
# ellipseResultFromOmegabh2VsOmegach2 = ellipse_params(omegabh2VsOmegach2FisherMatrix[0][0], omegabh2VsOmegach2FisherMatrix[1][1], omegabh2VsOmegach2FisherMatrix[0][1], confidenceA)

# Define the ellipse parameters
width = ellipseResultFromOmegabh2VsOmegach2[0]
height = ellipseResultFromOmegabh2VsOmegach2[1]
inclination = ellipseResultFromOmegabh2VsOmegach2[2]

# Calculate the angle of the semi-major axis
theta = inclination * np.pi / 180

# Calculate the x and y coordinates of the center of the ellipse
x0, y0 = 0.022, 0.12

# Create the ellipse object
ellipse = Ellipse(xy=(x0, y0), width=width, height=height, angle=inclination, fill=False)

# Create the plot and add the ellipse to it
fig, ax = plt.subplots(figsize=(8,8))
ax.add_artist(ellipse)
ax.scatter(x0, y0, color='red')

# Set the plot limits to ensure that the entire ellipse is visible
x_pad = width / 2
y_pad = height / 2
ax.set_xlim(x0 - width - x_pad, x0 + width + x_pad)
ax.set_ylim(y0 - height - y_pad, y0 + height + y_pad)

# Show the plot
plt.show()







