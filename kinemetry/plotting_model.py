#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 23:40:23 2019

@author: Duc Hoang, Rhodes College
"""

#Import necessary library
from astropy.io import fits
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def main(name):

    
    """
    The models and the residual has been saved in velModels.fits in the following order:
        - Index 0: Velocity model constructed by KINEMETRY
        - Index 1: Velocity residual (Subtract the model from the original velocity)
        - Index 2: Original Velocity map
        - Index 3: The Center's coordinates
        - Index 4: Radii of best fitted ellipses
        - Index 5: 1 - Ellipticity of the fitted ellipses
        - Index 6: Position Angle of the ellipses
        - Index 7: Signal to Noise analysis
        - Index 8: Name of the galaxy
    """
    #First open the fits file
    hdul = fits.open("velModel.fits")
    
    #Take the models out of the fits file as mentioned
    velModel = hdul[0].data
    residual = hdul[1].data
    velOriginal= hdul[2].data
    center = hdul[3].data
    radii = hdul[4].data
    q = hdul[5].data
    posAngle = hdul[6].data
    signalNoise = hdul[7].data
    
    
    #Plot all of graphs
    
    #Create the axes to be displayed (Spaxel to kpc)
    modelshapeX = velModel.shape[1]
    modelshapeY = velModel.shape[0]
    #Get the coordinates in physical units
    xtickPos, xtickVal = Coordinate_Conversion(center[1]+0.5,1.3,modelshapeX)
    ytickPos, ytickVal = Coordinate_Conversion(center[0]+0.5,1.3,modelshapeY)
    
    
    #1) Plot the Original Velocity map
    h = plt.figure(1)
    plt.figure(figsize = (12.4,8))
    sns.heatmap(velOriginal,cmap="RdBu_r",cbar_kws={'label': r'Velocity ($km \cdot s^{-1}$)'},vmin = -300, vmax = 300)
    title1 = "Original Velocity Map of " + name
    plt.title(title1)
    #X axis
    plt.xlabel("Distance From Nucleus (kpc)")
    plt.xticks(xtickPos,xtickVal)
    #Y axis
    plt.ylabel("Distance From Nucleus (kpc)")
    plt.yticks(ytickPos,ytickVal)
    h.show()
    
    #2) For the model, it also has to have the ellipses inscribed on the 2D map
    #Calculate the height and width for the ellipses (heights = radii, width = rad - (1-q)*rad)
    
    radii = (radii/0.8)*1.3
    width = np.zeros(radii.shape) #initialize the width array
    #Looping though to populate it
    for i in range(radii.shape[0]):
        width[i] = radii[i] - (1-q[i])*radii[i]
    
    #Plot the Model
    f, ax = plt.subplots(1)
    plt.figure(figsize = (12.4,8))
    ax = sns.heatmap(velModel,cmap="RdBu_r", cbar_kws={'label': r'Velocity ($km \cdot s^{-1}$)'},vmin = -300, vmax = 300)
    title2 = "Velocity map of " + name
    plt.title(title2)
    
    #X axis
    plt.xlabel("Distance From Nucleus (kpc)")
    plt.xticks(xtickPos,xtickVal)
    #Y axis
    plt.ylabel("Distance From Nucleus (kpc)")
    plt.yticks(ytickPos,ytickVal)
    
    #Plot the best fitted rings 
    for i in range(radii.shape[0]):
        ax.add_patch(patches.Ellipse((center[1]+0.5, center[0]+0.5), # Center Coordintate, the additional 0.5 account for the fact that the ticks are at the middle of each unit
                                     width[i], #width
                                     radii[i], #height
                                     angle = posAngle[i],
                                     alpha = 1, fill = False, 
                                     edgecolor="black", linewidth=1, linestyle='solid'))
    
    plt.show()
    
    #3) Plot the residual
    g = plt.figure(3)
    plt.figure(figsize = (12.4,8))
    sns.heatmap(residual,cmap="RdBu_r", cbar_kws={'label':  r'Velocity ($km \cdot s^{-1}$)'},vmin = -100, vmax = 100)
    plt.title("Residual between Original Velocity and generated Model")
    #X axis
    plt.xlabel("Distance From Nucleus (kpc)")
    plt.xticks(xtickPos,xtickVal)
    #Y axis
    plt.ylabel("Distance From Nucleus (kpc)")
    plt.yticks(ytickPos,ytickVal)
    g.show()
    
    #4) Plot the signal to noise ratio
    g = plt.figure(4)
    signalNoise = np.absolute(signalNoise)
    plt.figure(figsize = (12.4,8))
    sns.heatmap(signalNoise,cmap="Reds", cbar_kws={'label': 'Ratio'},vmin = 0, vmax = 2, annot = True)
    plt.title("Singal to Noise on each pixel")
    #X axis
    plt.xlabel("Distance From Nucleus (kpc)")
    plt.xticks(xtickPos,xtickVal)
    #Y axis
    plt.ylabel("Distance From Nucleus (kpc)")
    plt.yticks(ytickPos,ytickVal)
    g.show()


def Coordinate_Conversion(center, ratio, length):
    """
    Generate the coverted coordinates for the graphs (Originally it's pixel but now represent it in real physical units)
    
    Parameters:
        - Center: Contains x or y values of the center (Depending on the axis we're trying to convert)
        - Ratio: Ratio between the pixel and the real physical unit. 1 pixel = ... physical unit (I assume that it's kpc but it can be any)
        - Length: The length of the coordinate in pixel
    Returns: 
        - TickPos: The positions at which the ticks are placed (in pixel)
        - TickVal: The values at each tick.
    """
    
    unit = 1/ratio #Interval in pixel
    
    #Do a loop to create higher end of the axis
    higher_track = center #Initialize a tracker
    TickPos = np.array([center]) #Initialize the axis
    TickVal = np.array([0]) #Initialize labels for the ticks 
    i = 1 #Iterator
    #While loop to create higher end of the axis
    while higher_track < (length - unit):
        TickPos = np.append(TickPos,TickPos[-1]+unit)
        TickVal = np.append(TickVal,i)
        i += 1 #Increment couter
        higher_track += unit #Increment the tracker
    
    #Do a loop to create lower end of the axis
    lower_track = center #Initialize a tracker
    j = 1 #Iterator
    #While loop to create higher end of the axis
    while lower_track > unit:
        TickPos = np.insert(TickPos,0,TickPos[0]-unit)
        TickVal = np.insert(TickVal,0,j)
        j += 1 #Increment couter
        lower_track -= unit #Increment the tracker
        
    return TickPos, TickVal

main("PG1411")

