# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 07:55:58 2023

@author: Robert Gatt SM4391
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import matplotlib.transforms as transforms

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
    
def ellipse_params(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Calculate the width, height and inclination depending of the level of confidence"""
    a = _calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence)
    b = _calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence)
    thetaInDegrees = _calculateTanTwoTheta(sigmaX, sigmaY, sigmaXY)
    result = np.array([a, b, thetaInDegrees])
    print("(a) Width={0}, (b) Height={1}, (theta) Inclination={2}".format(a, b, thetaInDegrees))
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


# sigmaX = float(input("Enter sigmaX: "))
# sigmaY = float(input("Enter sigmaY: "))
# sigmaXY = float(input("Enter sigmaXY: "))
confidenceA = 1.52
confidenceB = 2.48

#[0][0] sigmaX, [0][1] sigma xy [1][0] sigma xy [1][1] sigma y
# ellipseResultFromUser = ellipse_params(sigmaX, sigmaY, sigmaXY, confidenceA)
ellipseResultFromOmegabh2VsOmegach2_A = ellipse_params(omegabh2VsOmegach2CovarianceMatrix[0][0],
                                                       omegabh2VsOmegach2CovarianceMatrix[1][1],
                                                       omegabh2VsOmegach2CovarianceMatrix[0][1], confidenceA)

ellipseResultFromOmegabh2VsOmegach2_B = ellipse_params(omegabh2VsOmegach2CovarianceMatrix[0][0],
                                                       omegabh2VsOmegach2CovarianceMatrix[1][1],
                                                       omegabh2VsOmegach2CovarianceMatrix[0][1], confidenceB)


#############################################################################################################

# Define the ellipse parameters
width_A = ellipseResultFromOmegabh2VsOmegach2_A[0]
height_A = ellipseResultFromOmegabh2VsOmegach2_A[1]
inclination_A = ellipseResultFromOmegabh2VsOmegach2_A[2]

width_B = ellipseResultFromOmegabh2VsOmegach2_B[0]
height_B = ellipseResultFromOmegabh2VsOmegach2_B[1]
inclination_B = ellipseResultFromOmegabh2VsOmegach2_B[2]

x, y = 0.022, 0.12
mean = np.array([x, y])
points = np.random.multivariate_normal(mean, omegabh2VsOmegach2CovarianceMatrix, size=500)
ellipseOneSigma = Ellipse(xy=(x, y), width=width_A, height=height_A, angle=180-inclination_A, fill=False, color='blue')
ellipseTwoSigma = Ellipse(xy=(x, y), width=width_B, height=height_B, angle=180-inclination_B, fill=False, color='blue', linestyle = "dotted")

print(points[:,0])

# Create the plot and add the ellipse to it
fig, ax = plt.subplots()
ax.add_patch(ellipseOneSigma)
ax.add_artist(ellipseTwoSigma)

# Plot the random points
plt.scatter(points[:,0], points[:,1], color='black', s=1)

ax.plot(x, y, "ro")

ax.set_xlim(xmin=x-width_A*2, xmax=x+width_A*2)
ax.set_ylim(ymin=y-height_A*2, ymax=y+height_A*2)

# ax.set_xlim(xmin=min(x-width_A, x-width_B), xmax=max(x+width_A, x+width_B))
# ax.set_ylim(ymin=min(y-height_A, y-height_B), ymax=max(y+height_A, y+height_B))
# Show the plot
plt.show()






