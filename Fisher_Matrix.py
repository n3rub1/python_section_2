# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 07:55:58 2023

@author: Robert Gatt SM4391
"""

"""
The two required packages are imported first by the program; they are all 
included in the Python’s standard library.
•	Matplotlib: used to plot the error ellipses and the Gaussian 500 random points.
•	NumPy: used throughout the program to work various mathematical formulas 
    such as the inverse of the covariance matrix, square roots, inverse tangent 
    and converting arrays to matrices. 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

"""
This program consists of various function, each of which has a specific purpose 
and is designed to follow the principle of single-responsibility programming.
•	getHeaderVariables
•	getFullDataMatrix
•	printCovarianceAndFisherValues
•	getFisherMatrix
•	getCovarianceMatrix
•	_calculateInnerFormula
•	_calculateASquared
•	_calculateBSqaured
•	_calculateTanTwoTheta
•	ellipse_params
•	printEllipse_params
•	plot
"""

def getHeaderVariables():
    """get the list of all the variable names from the file"""
    with open("base_w_plikHM_TTTEEE_lowl_lowE_BAO_Riess18_Pantheon.covmat", "r") as file:
        headerLine = file.readline().strip()
        variableList = headerLine.split("\t") 
    return variableList


def getFullDataMatrix():
    """This method is used to bring in the original matrix"""  
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
    

def getFisherMatrix(covarianceMatrix):
    """transforms the covariance matrix into a fisher matrix and returns"""
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


def ellipse_params(sigmaX, sigmaY, sigmaXY, levelOfConfidence):
    """Calculate the width, height and inclination depending of the level of confidence"""
    a = _calculateASquared(sigmaX, sigmaY, sigmaXY, levelOfConfidence)
    b = _calculateBSqaured(sigmaX, sigmaY, sigmaXY, levelOfConfidence)
    thetaInDegrees = _calculateTanTwoTheta(sigmaX, sigmaY, sigmaXY)
    result = np.array([b, a, thetaInDegrees])
    return result


def printEllipse_params(params, name):
    """prints the width, height and inclination for all ellipses with α = 1.52"""
    print("{0}  at α = 1.52 \n height (a): {1}, width(b): {2}, inclination in degrees(θ): {3}\n".format(name, params[0], params[1], params[2]))
    
"""
First, the function creates two error ellipses using the Ellipse() class from
matplotlib library.  Afterwards, it generates 500 random points using the NumPy 
random.multivariate_normal method based on the covariance matrix.
Next, the function adds the random points and ellipses to the plot and adds a 
red dot to mark the center of the ellipses. It creates a legend in the top 
right-hand corner of the screen to show the two standard deviations.
Finally, the function plots the graphs using the plot.show() method.
"""

def plot(widthA, widthB, heightA, heightB, inclinationA, inclinationB, x, y, covarianceMatrix, variableNameX, variableNameY):
    """Plot the ellipses and the generated points"""
    
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

    plt.show()
 
    
#initilalize the program and fill up the matrix values
variableNames = getHeaderVariables()
fullDataMatrix = getFullDataMatrix()

#get the covariance matrices of the requested values
omegabh2VsOmegach2CovarianceMatrix = getCovarianceMatrix("omegabh2", "omegach2", fullDataMatrix)
omegach2VsWCovarianceMatrix = getCovarianceMatrix("omegach2", "w", fullDataMatrix)
logAVsNsCovarianceMatrix = getCovarianceMatrix("logA", "ns", fullDataMatrix)
tauVsWCovarianceMatrix = getCovarianceMatrix("tau", "w", fullDataMatrix)

#get the fisher matrices of the requested values
omegabh2VsOmegach2FisherMatrix = getFisherMatrix(omegabh2VsOmegach2CovarianceMatrix)
omegach2VsWFisherMatrix = getFisherMatrix(omegach2VsWCovarianceMatrix)
logAVsNsFisherMatrix = getFisherMatrix(logAVsNsCovarianceMatrix)
tauVsWFisherMatrix = getFisherMatrix(tauVsWCovarianceMatrix)

#print converiance matrix and fisher matrix
printCovarianceAndFisherValues(omegabh2VsOmegach2CovarianceMatrix, omegabh2VsOmegach2FisherMatrix, "omegabh2", "omegach2")
printCovarianceAndFisherValues(omegach2VsWCovarianceMatrix, omegach2VsWFisherMatrix, "omegach2", "w")
printCovarianceAndFisherValues(logAVsNsCovarianceMatrix, logAVsNsFisherMatrix, "logA", "ns")
printCovarianceAndFisherValues(tauVsWCovarianceMatrix, tauVsWFisherMatrix, "tau", "w")

#set the confidence levels as final variables
CONFIDENCE_A = 1.52
CONFIDENCE_B = 2.48

#get the ellipse parameters
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

#print the ellipse parameters
printEllipse_params(ellipseResultFromOmegabh2VsOmegach2_A, "Ωbh2 vs Ωch2")
printEllipse_params(ellipseResultFromOmegach2VsW_A, "Ωch2 vs w")
printEllipse_params(ellipseResultFromlogAVsNs_A, "ln(A) vs ns")
printEllipse_params(ellipseResultFromtauVsW_A, "τ vs w")

#plot the ellipses
plot(ellipseResultFromOmegabh2VsOmegach2_A[0], ellipseResultFromOmegabh2VsOmegach2_B[0],
     ellipseResultFromOmegabh2VsOmegach2_A[1], ellipseResultFromOmegabh2VsOmegach2_B[1],
     ellipseResultFromOmegabh2VsOmegach2_A[2], ellipseResultFromOmegabh2VsOmegach2_B[2],
     0.022, 0.12, omegabh2VsOmegach2CovarianceMatrix, "Ωbh2", "Ωch2")

plot(ellipseResultFromOmegach2VsW_A[0], ellipseResultFromOmegach2VsW_B[0],
     ellipseResultFromOmegach2VsW_A[1], ellipseResultFromOmegach2VsW_B[1],
     ellipseResultFromOmegach2VsW_A[2], ellipseResultFromOmegach2VsW_B[2],
     0.12, -1, omegach2VsWCovarianceMatrix, "Ωch2", "w")

plot(ellipseResultFromlogAVsNs_A[0], ellipseResultFromlogAVsNs_B[0],
     ellipseResultFromlogAVsNs_A[1], ellipseResultFromlogAVsNs_B[1],
     ellipseResultFromlogAVsNs_A[2], ellipseResultFromlogAVsNs_B[2],
     -19.94, 0.96, logAVsNsCovarianceMatrix, "ln(A)", "ns")

plot(ellipseResultFromtauVsW_A[0], ellipseResultFromtauVsW_B[0],
     ellipseResultFromtauVsW_A[1], ellipseResultFromtauVsW_B[1],
     ellipseResultFromtauVsW_A[2], ellipseResultFromtauVsW_B[2],
     0.09, -1, tauVsWCovarianceMatrix, "τ", "w")

#check if user would like to enter their own matrix to get the values
isMenu = True

while(isMenu):
    userInput = int(input("Would you like to enter your own parameters?\n 1: Yes \n 2: No \n"))
    if(userInput == 1):
        userSigmaX = float(input("Enter σx"))
        userSigmaY = float(input("Enter σy"))
        userSigmaXY = float(input("Enter σxy"))
        userLevelOfConfidence = float(input("Enter level of confidence"))       
        userEllipseResult = ellipse_params(userSigmaX, userSigmaY, userSigmaXY, userLevelOfConfidence)
        printEllipse_params(userEllipseResult, "User Inputs")       
    elif(userInput == 2):
        print("Thank you and goodbye")
        isMenu = False
    else:
        print("invalid input")
