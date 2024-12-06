# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:16:07 2022

@author: ritzkowf
"""
import sys
import os
import shutil

import pandas as pd
import time
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import seabreeze
from seabreeze.spectrometers import Spectrometer, list_devices
import shelve as shv
import os
import elliptec
import avaSpectrometer as ava 

from scipy import signal

from scipy.signal import find_peaks
# import gaussian filter 1d
from scipy.ndimage import gaussian_filter1d

sys.path.append(os.path.abspath(r"C:\PythonCoding\refractive-index-database"))



def alignment(intTime):
    spec = ava.spectrometer()
    spec.set_integration_time(intTime)

    # Get wavelength axis
    wavelength = spec.get_wavelengths()

    idxLow = np.abs(wavelength-500).argmin()
    idxHigh = np.abs(wavelength-800).argmin()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    
    line, = ax.plot(wavelength[idxLow:idxHigh],spec.intensities()[idxLow:idxHigh])
    ax.set_ylim([0,65000])
    ax.set_xlabel("Wavelength (nm)")

 

    try:
        while True:
            line.set_ydata(spec.intensities()[idxLow:idxHigh])
            fig.canvas.draw()        
            fig.canvas.flush_events()
     
    except KeyboardInterrupt:
        pass
    line.set_ydata(spec.intensities()[idxLow:idxHigh])
    fig.canvas.draw()        
    fig.canvas.flush_events()
     
    
    spec.close()
     
     
     
    return


def main():
    dirName = r"C:\Users\ritzk\OneDrive\Desktop\Data\dscan\25102024\Measurement_"
    i=0

    while os.path.exists(dirName+f"{i}"):
        i += 1
    dirName = dirName + f"{i}"
    os.makedirs(dirName, exist_ok=True)
    
    dirPath = dirName + "/"
    





    intTime =3
    avgFactor = 20



  




    
    
    posList = np.arange(5,40,0.2)
    

    spec = ava.spectrometer()
    spec.set_integration_time(intTime)

    # Get wavelength axis
    wavelength = spec.get_wavelengths()

    # Get wavelength axis

    idxLow = np.abs(wavelength-500).argmin()
    idxHigh = np.abs(wavelength-850).argmin()
    

    # Initialize the controller with the correct device name
    controller = elliptec.Controller('COM4')  # Adjust device name if needed
    ls = elliptec.Linear(controller)

    # Home the linear stage before usage
    ls.home()
    
    print("Please blockthe beam!")
    print("When done, press enter:")
    input("..")
    
    darkspectrum = spec.get_spectrum()
    
    print("Please unblock the beam!")
    print("When done, press enter:")
    input("..")
    
    # Init data array
    data = np.zeros([len(posList),len(wavelength)])
    
    # open live plot
    fig = plt.figure()
    ax =  plt.subplot(111)
    im = ax.imshow(data[:,idxLow:idxHigh],aspect = 'auto',origin='lower',extent=[500,800,posList[0],posList[-1]])
    im.set_clim(vmin=0, vmax=65000)
    timeAxis = np.zeros(len(posList))
    startT = time.time()    
    for i in range(len(posList)):
        ls.set_distance(posList[i])
        time.sleep(0.05)
        # averaging over 10 spectra
        for j in range(avgFactor):
            data[i,:] += spec.get_spectrum()
        data[i,:] /= avgFactor
        #data[i,:] = spec.get_spectrum()
        timeAxis[i]=time.time()
        print(f" Progress {100*(i/len(posList))}% " )
        im.set_data(data[:,idxLow:idxHigh])
        fig.canvas.draw()        
        fig.canvas.flush_events()
    
    IntData = np.sum(data[:,idxLow:idxHigh],axis=1)
    IntDataSmoothed = gaussian_filter1d(IntData, sigma=6)
    peak =np.abs(IntDataSmoothed).argmax()
    ls.set_distance(posList[peak])


    

    del(ls)
    del(spec)
    plt.savefig(dirPath + "data.png")
    pd.DataFrame(data).to_csv(dirPath +"data.csv")
    pd.DataFrame(wavelength).to_csv(dirPath +"wavelength.csv")
    pd.DataFrame(posList).to_csv(dirPath +"posList.csv")
    pd.DataFrame(timeAxis).to_csv(dirPath +"timeAxis.csv")
    pd.DataFrame(darkspectrum).to_csv(dirPath +"darkspectrum.csv")

    shutil.copy(__file__, dirPath+'measurementScript.py') 
    
        # fig1, ax1 = plt.subplots(1,1)
        
        # ax1.plot(wavelength[idxLow:idxHigh],np.sum(data[:,idxLow:idxHigh]),axis=0)
        # ax1.set_ylim([-200,200])
        # secax1 = ax1.secondary_xaxis('top', functions =(NIR2MIR,MIR2NIR))
        # secax.set_xlabel('Wavelength MIR (nm)')
        
    return




if __name__ == '__main__':
    main()

