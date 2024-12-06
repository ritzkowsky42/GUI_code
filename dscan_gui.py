"""
This script should generate a graphical user interface that allows a user to acquire dscan measurements
Author: Adina Bechhofer
Contact for support: adinabec [at] mit.edu
"""

import tkinter
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,  NavigationToolbar2Tk 
import numpy as np
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

def confirm_beam_blocked():
    darkspectrum = spec.get_spectrum()
    print("unblock")
    Text_block = tkinter.Text(root, height = 5, width = 52)
    l = Label(root, text = "Attention!").grid(row=6)
    Text_block.grid(row=7, column=1)
    Text_block .insert(tk.END, "Please unblock the beam! click confirm when your are done.")
    unblocked =tkinter.Button(root,text = 'Confirm', command = confirm_beam_unblocked)
    unblocked.grid(row=8, column=1)

def confirm_beam_unblocked():
    # Init data array
    data = np.zeros([len(sys_vars["posList"]),len(sys_vars["wavelength"])])
    
    # open live plot
    fig = plt.figure()
    ax =  plt.subplot(111)
    im = ax.imshow(data[:,idxLow:idxHigh],aspect = 'auto',origin='lower',extent=[500,800,sys_vars["posList"][0],sys_vars["posList"][-1]])
    im.set_clim(vmin=0, vmax=65000)
    timeAxis = np.zeros(len(sys_vars["posList"]))
    canvas = FigureCanvasTkAgg(fig, 
                               master = root)
    canvas.get_tk_widget().grid(row=9,column=2,columnspan=10,rowspan=10)   
    toolbarFrame = tkinter.Frame(master=root)
    toolbarFrame.grid(row=18,column=2)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    startT = time.time()    
    for i in range(len(sys_vars["posList"])):
        ls.set_distance(sys_vars["posList"][i])
        time.sleep(0.05)
        # averaging over 10 spectra
        for j in range(sys_vars["avgFactor"]):
            data[i,:] += spec.get_spectrum()
        data[i,:] /= sys_vars["avgFactor"]
        #data[i,:] = spec.get_spectrum()
        timeAxis[i]=time.time()
        print(f" Progress {100*(i/len(sys_vars["posList"]))}% " )
        im.set_data(data[:,idxLow:idxHigh])
        canvas.draw()        
        fig.canvas.flush_events()
    
    IntData = np.sum(data[:,idxLow:idxHigh],axis=1)
    IntDataSmoothed = gaussian_filter1d(IntData, sigma=6)
    peak =np.abs(IntDataSmoothed).argmax()
    ls.set_distance(sys_vars["posList"][peak])

    del(ls)
    del(spec)
    plt.savefig(sys_vars["dirPath"] + "data.png")
    pd.DataFrame(data).to_csv(sys_vars["dirPath"] +"data.csv")
    pd.DataFrame(wavelength).to_csv(sys_vars["dirPath"] +"wavelength.csv")
    pd.DataFrame(posList).to_csv(sys_vars["dirPath"] +"posList.csv")
    pd.DataFrame(timeAxis).to_csv(sys_vars["dirPath"] +"timeAxis.csv")
    pd.DataFrame(darkspectrum).to_csv(sys_vars["dirPath"] +"darkspectrum.csv")

    shutil.copy(__file__, sys_vars["dirPath"]+'measurementScript.py') 
    
        # fig1, ax1 = plt.subplots(1,1)
        
        # ax1.plot(wavelength[idxLow:idxHigh],np.sum(data[:,idxLow:idxHigh]),axis=0)
        # ax1.set_ylim([-200,200])
        # secax1 = ax1.secondary_xaxis('top', functions =(NIR2MIR,MIR2NIR))
        # secax.set_xlabel('Wavelength MIR (nm)')
        
    return


def get_inputs():
    ##  get directory name
    dirName = dirName_var.get()

    while os.path.exists(dirName+f"{i}"):
        i += 1
    dirName = dirName + f"{i}"
    os.makedirs(dirName, exist_ok=True)
    
    # calculate new path
    dirPath = dirName + "/"
    intTime = intTime_var.get()
    avgFactor = avgFactor_var.get()

    sys_vars["dirPath"] = dirPath
    sys_vars["intTime"] = intTime
    avgFactor["avgFactor"] = avgFactor
    sys_vars["dirName"] = dirName

    print("entered variables:") 
    print("dir:", dirPath)
    print("intTime:", intTime, ", avgFactor:", avgFactor)
    
    spec.set_integration_time(sys_vars["intTime"])
    
    dirName_var.set("")


def main():

    posList = np.arange(5,40,0.2)
    sys_vars["posList"] = posList
    
    spec = ava.spectrometer()
    get_inputs() # this also sets the spectrometer integration time
    
    # Get wavelength axis
    wavelength = spec.get_wavelengths()
    sys_vars["wavelength"] = wavelength

    # Get wavelength axis

    idxLow = np.abs(wavelength-500).argmin() ## consider changing 
    idxHigh = np.abs(wavelength-850).argmin()
    

    # Initialize the controller with the correct device name
    controller = elliptec.Controller('COM4')  # Adjust device name if needed
    ls = elliptec.Linear(controller)

    # Home the linear stage before usage
    ls.home()
    
    Text_block = tkinter.Text(root, height = 5, width = 52)
    l = Label(root, text = "Attention!").grid(row=6)
    Text_block.grid(row=7, column=1)
    Text_block .insert(tk.END, "Please block the beam! click confirm when your are done")
    blocked =tkinter.Button(root,text = 'Confirm', command = confirm_beam_blocked)
    blocked.grid(row=8,column=1)
    

if __name__ == '__main__':
    root = tkinter.Tk()
    root.geometry("900x700")
    ## Don't know if we need menur yet, but worth keeping asround
    menu = tkinter.Menu(root)
    root.config(menu=menu)
    filemenu = tkinter.Menu(menu)
    menu.add_cascade(label='File', menu=filemenu)
    filemenu.add_command(label='New')
    filemenu.add_command(label='Open...')
    filemenu.add_separator()
    filemenu.add_command(label='Exit', command=root.quit)
    helpmenu = tkinter.Menu(menu)
    menu.add_cascade(label='Help', menu=helpmenu)
    helpmenu.add_command(label='About')

    dirName_var = tkinter.StringVar()
    intTime_var = tkinter.DoubleVar()
    avgFactor_var = tkinter.IntVar()

    sys_vars = {}

    tkinter.Label(root, text='Directory name').grid(row=0)
    tkinter.Label(root, text='Integration time (s)').grid(row=1)
    tkinter.Label(root, text='Averaging factor').grid(row=2)
    e1 = tkinter.Entry(root, textvariable=dirName_var)
    e2 = tkinter.Entry(root, textvariable=intTime_var)
    e3 = tkinter.Entry(root, textvariable=avgFactor_var)
    e1.grid(row=0, column=1)
    e2.grid(row=1, column=1)
    e3.grid(row=2, column=1)

    sub_btn=tkinter.Button(root,text = 'Start acquisition', command = main)
    sub_btn.grid(row=5,column=1)

    tkinter.mainloop()


  