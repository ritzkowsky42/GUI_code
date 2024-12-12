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
import shutil
import pandas as pd
import time
import matplotlib.pyplot as plt
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

def browse_location():
    filename = tkinter.filedialog.askdirectory()
    dirName_var.set(filename)
    print("set filename")

def confirm_beam_blocked():
    blocked.grid_forget()
    sys_vars["darkspectrum"] = sys_vars["spec"].get_spectrum()
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
    im = ax.imshow(data[:,sys_vars["idxLow"]:sys_vars["idxHigh"]],aspect = 'auto',origin='lower',extent=[500,800,sys_vars["posList"][0],sys_vars["posList"][-1]])
    im.set_clim(vmin=0, vmax=65000)
    
    timeAxis = np.zeros(len(sys_vars["posList"]))
    canvas.get_tk_widget().grid(row=8,column=1,columnspan=10,rowspan=10)   
    toolbarFrame = tkinter.Frame(master=root)
    toolbarFrame.grid(row=17,column=1)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    startT = time.time()    
    for i in range(len(sys_vars["posList"])):
        sys_vars["ls"].set_distance(sys_vars["posList"][i])
        time.sleep(0.05)
        # averaging over 10 spectra
        for j in range(sys_vars["avgFactor"]):
            data[i,:] += sys_vars["spec"].get_spectrum()
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
    sys_vars["ls"].set_distance(sys_vars["posList"][peak])

    plt.savefig(sys_vars["dirPath"] + "data.png")
    pd.DataFrame(data).to_csv(sys_vars["dirPath"] +"data.csv")
    pd.DataFrame(sys_vars["wavelength"]).to_csv(sys_vars["dirPath"] +"wavelength.csv")
    pd.DataFrame(sys_vars["posList"]).to_csv(sys_vars["dirPath"] +"posList.csv")
    pd.DataFrame(timeAxis).to_csv(sys_vars["dirPath"] +"timeAxis.csv")
    pd.DataFrame(sys_vars["darkspectrum"]).to_csv(sys_vars["dirPath"] +"darkspectrum.csv")

    shutil.copy(__file__, sys_vars["dirPath"]+'measurementScript.py') 
    repeat.grid(row=18, column=0)
    dscan.grid(row=18, column=1)
    
        # fig1, ax1 = plt.subplots(1,1)
        
        # ax1.plot(wavelength[idxLow:idxHigh],np.sum(data[:,idxLow:idxHigh]),axis=0)
        # ax1.set_ylim([-200,200])
        # secax1 = ax1.secondary_xaxis('top', functions =(NIR2MIR,MIR2NIR))
        # secax.set_xlabel('Wavelength MIR (nm)')
        
    return

def reset_all():
    #dirName_var.set("")
    #intTime_var.set(0)
    #avgFactor_var.set(0)
    start_btn.grid(row=5,column=1)
    progress.grid_forget()
    sys_vars["canvas"].get_tk_widget().grid_forget()
    sys_vars["toolbarFrame"] .grid_forget()
    repeat.grid_forget()
    dscan.grid_forget()
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
    min_insert = min_insert_var.get()
    max_insert = max_insert_var.get()

    sys_vars["dirPath"] = dirPath
    sys_vars["intTime"] = intTime
    sys_vars["avgFactor"] = avgFactor
    sys_vars["dirName"] = dirName
    sys_vars["min_insert"] = min_insert
    sys_vars["max_insert"] = max_insert

    print("entered variables:") 
    print("dir:", dirPath)
    print("intTime:", intTime, ", avgFactor:", avgFactor)
    
    spec.set_integration_time(sys_vars["intTime"])
    start_btn.grid_forget()
    #dirName_var.set("")

def on_gui_close():
    if "ls" in sys_vars:
        del sys_vars["ls"]
    if "spec" in sys_vars:
        del sys_vars["spec"]

def main():

    global spec 
    spec = ava.spectrometer()
    sys_vars["spec"] = spec
    
    get_inputs() # this also sets the spectrometer integration time
    posList = np.arange(sys_vars["min_insert"], sys_vars["max_insert"], 0.2)
    sys_vars["posList"] = posList
    
    
    # Get wavelength axis
    wavelength = spec.get_wavelengths()
    sys_vars["wavelength"] = wavelength

    # Get wavelength axis
    idxLow = np.abs(wavelength-500).argmin() ## consider changing 
    idxHigh = np.abs(wavelength-850).argmin()
    sys_vars["idxLow"] = idxLow 
    sys_vars["idxHigh"] = idxHigh 
    

    # Initialize the controller with the correct device name
    controller = elliptec.Controller('COM4')  # Adjust device name if needed
    global ls

    ls = elliptec.Linear(controller)
    sys_vars["ls"] = ls

    # Home the linear stage before usage
    ls.home()
    
    attn_label.grid(row=6)
    Text_block.grid(row=7, column=1)
    Text_block.insert(tkinter.END, "Please block the beam! \nClick confirm when your are done.")
    blocked.grid(row=8,column=1)
    
def do_dscan_retrieval():
    print("Doing D-scan retreival")
    pass

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
    min_insert_var = tkinter.IntVar()
    max_insert_var = tkinter.IntVar()

    sys_vars = {}

    tkinter.Label(root, text='Directory name').grid(row=0)
    browse_btn = tkinter.Button(root, text="Browse", command=browse_location).grid(row=0, column=2)
    tkinter.Label(root, text='Integration time (s)').grid(row=1)
    tkinter.Label(root, text='Averaging factor').grid(row=2)
    tkinter.Label(root, text='Wedge insertion (mm)').grid(row=3)
    tkinter.Label(root, text='min', justify="right").grid(row=3, column=1)
    tkinter.Label(root, text='max').grid(row=3, column=3)

    e1 = tkinter.Entry(root, textvariable=dirName_var)
    e2 = tkinter.Entry(root, textvariable=intTime_var)
    e3 = tkinter.Entry(root, textvariable=avgFactor_var)
    e4 = tkinter.Entry(root, textvariable=min_insert_var, width=5)
    e5 = tkinter.Entry(root, textvariable=max_insert_var, width=5)
    e1.grid(row=0, column=1)
    e2.grid(row=1, column=1)
    e3.grid(row=2, column=1)
    e4.grid(row=3, column=2)
    e5.grid(row=3, column=4)

    start_btn=tkinter.Button(root,text = 'Start acquisition', command = main)
    start_btn.grid(row=5,column=1)

    blocked = tkinter.Button(root, text = 'Confirm', command = confirm_beam_blocked)
    attn_label = tkinter.Label(root, text = "Attention!", font=('Helvetica',12, 'bold'))
    Text_block = tkinter.Text(root, height = 5, width = 45)
    unblocked =tkinter.Button(root,text = 'Confirm', command = confirm_beam_unblocked)
    repeat = tkinter.Button(root,text = 'New acquisition', command = reset_all)
    dscan = tkinter.Button(root, text = 'Phase retrival', command = do_dscan_retrieval)
    progress = ttk.Progressbar(root, length=100, mode='determinate')
    progress_label = tkinter.Label(root, text='Progress')

    ## initialize default insertion width
    min_insert_var.set(5)
    max_insert_var.set(40)
    root.protocol("WM_DELETE_WINDOW", on_gui_close)
    
    tkinter.mainloop()


# #FXIME make closing 
#     del(sys_vars["ls"])
#     del(sys_vars["spec"])
  