"""
This script should generate a graphical user interface that allows a user to acquire dscan measurements
Author: Adina Bechhofer
Contact for support: adinabec [at] mit.edu
"""

import tkinter
from tkinter import ttk
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,  NavigationToolbar2Tk 
import numpy as np
import sys
import os

import pandas as pd
import time
import numpy as np
import matplotlib.pyplot as plt


from scipy import signal

from scipy.signal import find_peaks
# import gaussian filter 1d
from scipy.ndimage import gaussian_filter1d

sys.path.append(os.path.abspath(r"C:\PythonCoding\refractive-index-database"))

def browse_location():
    filename = tkinter.filedialog.askdirectory()
    dirName_var.set(filename)
    print("set filename")


def confirm_beam_blocked():
    print("Get dark spectrum")
    blocked.grid_forget()
    sys_vars["darkspectrum"] = np.random.rand(3, 8)
    print("unblock")
    Text_block.delete('1.0', tkinter.END)
    Text_block.insert(tkinter.END, "Please unblock the beam! \nClick confirm when your are done.")
    unblocked.grid(row=8, column=1)

def confirm_beam_unblocked():
    # hide buttons and features from before 
    attn_label.grid_forget()
    Text_block.grid_forget()
    unblocked.grid_forget()
    # Init data array
    data = np.zeros([len(sys_vars["posList"]),len(sys_vars["wavelength"])])
    
    # make progress bar appear
    progress_label.grid(row=7, column=0)
    progress.grid(row=7, column=1)
    # open live plot
    fig = plt.figure()
    ax =  plt.subplot(111)
    canvas = FigureCanvasTkAgg(fig, 
                               master = root)
    im = ax.imshow(data[:, sys_vars["idxLow"]:sys_vars["idxHigh"]],aspect = 'auto',origin='lower',extent=[500,800,sys_vars["posList"][0],sys_vars["posList"][-1]])
    im.set_clim(vmin=0, vmax=65000)
    
    timeAxis = np.zeros(len(sys_vars["posList"]))
    canvas.get_tk_widget().grid(row=8,column=1,columnspan=10,rowspan=10)   
    toolbarFrame = tkinter.Frame(master=root)
    toolbarFrame.grid(row=17,column=1)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    startT = time.time()    
    for i in range(len(sys_vars["posList"])):
        print("ls.set_distance(poslist[i])")
        time.sleep(0.05)
        # averaging over 10 spectra
        for j in range(sys_vars["avgFactor"]):
            data[i,:] += np.random.rand((data.shape[-1]))*65000
        data[i,:] /= sys_vars["avgFactor"]
        #data[i,:] = spec.get_spectrum()
        timeAxis[i]=time.time()
        print(f" Progress {100*(i/len(sys_vars['posList']))}% " )
        progress["value"] = 100*(i/len(sys_vars['posList']))
        im.set_data(data[:,sys_vars["idxLow"]:sys_vars["idxHigh"]])
        canvas.draw()        
        fig.canvas.flush_events()

    sys_vars["canvas"] = canvas
    sys_vars["toolbarFrame"] =toolbarFrame
    progress["value"] = 100

    IntData = np.sum(data[:,sys_vars["idxLow"]:sys_vars["idxHigh"]],axis=1)
    IntDataSmoothed = gaussian_filter1d(IntData, sigma=6)
    peak =np.abs(IntDataSmoothed).argmax()
    print("ls.set_distance(sys_vars[posList][peak])")

    print("del(ls)")
    print("del(spec)")
    plt.savefig(sys_vars["dirPath"] + "data.png")
    pd.DataFrame(data).to_csv(sys_vars["dirPath"] +"data.csv")
    pd.DataFrame(sys_vars["wavelength"]).to_csv(sys_vars["dirPath"] +"wavelength.csv")
    pd.DataFrame(sys_vars["posList"]).to_csv(sys_vars["dirPath"] +"posList.csv")
    pd.DataFrame(timeAxis).to_csv(sys_vars["dirPath"] +"timeAxis.csv")
    pd.DataFrame(sys_vars["darkspectrum"]).to_csv(sys_vars["dirPath"] +"darkspectrum.csv")

    print("Shuthil")
    repeat.grid(row=18, column=10)
    
        # fig1, ax1 = plt.subplots(1,1)
        
        # ax1.plot(wavelength[idxLow:idxHigh],np.sum(data[:,idxLow:idxHigh]),axis=0)
        # ax1.set_ylim([-200,200])
        # secax1 = ax1.secondary_xaxis('top', functions =(NIR2MIR,MIR2NIR))
        # secax.set_xlabel('Wavelength MIR (nm)')
        
    return

def reset_all():
    dirName_var.set("")
    intTime_var.set(0)
    avgFactor_var.set(0)
    start_btn.grid(row=5,column=1)
    progress.grid_forget()
    sys_vars["canvas"].get_tk_widget().grid_forget()
    sys_vars["toolbarFrame"] .grid_forget()
    repeat.grid_forget()
    progress_label.grid_forget()


def get_inputs():
    ##  get directory name
    dirName = dirName_var.get()
    i = 0
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
    sys_vars["avgFactor"] = avgFactor
    sys_vars["dirName"] = dirName

    print("entered variables:") 
    print("dir:", dirPath)
    print("intTime:", intTime, ", avgFactor:", avgFactor)
    
    print("Set spectrometer int time")
    start_btn.grid_forget()
    #dirName_var.set("")


def main():

    posList = np.arange(5,40,0.2)
    sys_vars["posList"] = posList
    
    print("initialize spectrometer ")
    get_inputs() # this also sets the spectrometer integration time
    
    # Get wavelength axis
    print("get wavelength range from spectrometers")
    sys_vars["wavelength"] =  np.linspace(600, 2000, 100)

    # Get wavelength axis
    idxLow = np.abs(sys_vars["wavelength"] -500).argmin() ## consider changing 
    idxHigh = np.abs(sys_vars["wavelength"] -850).argmin()
    sys_vars["idxLow"] = idxLow 
    sys_vars["idxHigh"] = idxHigh 
    

    # Initialize the controller with the correct device name
    print("set controller")
    print("assign controller")

    # Home the linear stage before usage
    print("ls.home()")
    
    attn_label.grid(row=6)
    Text_block.grid(row=7, column=1)
    Text_block.insert(tkinter.END, "Please block the beam! \nClick confirm when your are done.")
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
    browse_btn = tkinter.Button(root, text="Browse", command=browse_location).grid(row=0, column=2) 
    tkinter.Label(root, text='Integration time (s)').grid(row=1)
    tkinter.Label(root, text='Averaging factor').grid(row=2)
    e1 = tkinter.Entry(root, textvariable=dirName_var)
    e2 = tkinter.Entry(root, textvariable=intTime_var)
    e3 = tkinter.Entry(root, textvariable=avgFactor_var)
    e1.grid(row=0, column=1)
    e2.grid(row=1, column=1)
    e3.grid(row=2, column=1)

    start_btn=tkinter.Button(root,text = 'Start acquisition', command = main)
    start_btn.grid(row=5,column=1)

    blocked = tkinter.Button(root, text = 'Confirm', command = confirm_beam_blocked)
    attn_label = tkinter.Label(root, text = "Attention!")
    Text_block = tkinter.Text(root, height = 5, width = 60)
    unblocked =tkinter.Button(root,text = 'Confirm', command = confirm_beam_unblocked)
    repeat =tkinter.Button(root,text = 'Done', command = reset_all)
    progress = ttk.Progressbar(root, length=100, mode='determinate')
    progress_label = tkinter.Label(root, text='Progress')

    tkinter.mainloop()


  