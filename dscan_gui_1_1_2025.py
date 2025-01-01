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
import pypret
from scipy.ndimage import gaussian_filter1d
# importing 2d filter
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
# import gaussian filter 1d
## Import pypret stuff
from scipy.fftpack import next_fast_len
from pypret.fourier import FourierTransform
from scipy.interpolate import interp1d

class HardWare:
    def __init__(self, spec = None, lin_stage = None):
        self.spec = spec
        self.lin_stage = lin_stage
        self.pos_list = None
        self.wavelength = None
        self.idxLow = None
        self.idxHigh = None
        self.darkspectrum = None
    
    def initialize_spec_and_lin_stage(self, spec, lin_stage):
        self.spec = spec
        self.lin_stage = lin_stage

    def set_spec_int_time(self, intTime):
        self.spec.set_integration_time(intTime)

    def set_position_list(self, minpos, maxpos):
        self.pos_list = np.arange(minpos, maxpos, 0.2)

    def get_spec_wavelength(self):
        """
        TODO: This method name might be ambigious. People will expect the method to return the wavelength. 
        """
        self.wavelength = self.spec.get_wavelengths()
        self.idxLow = np.abs(self.wavelength-500).argmin() ## consider changing 
        self.idxHigh = np.abs(self.wavelength-850).argmin()

    def home_lin_stage(self):
        self.lin_stage.home()

    def set_dark_spectrum(self):
        self.darkspectrum = self.spec.get_spectrum()

    def dataset_dim(self):
        return [len(self.pos_list),len(self.wavelength)]   
    
    def get_low_and_high_wavelenghth_inds(self):
        return self.idxLow, self.idxHigh
    
    def measure_spectrum(self):
        return self.spec.get_spectrum()
    
    def set_linstage_dist(self, dist):
        self.lin_stage.set_distance(dist)
    
    def __del__(self):
        del(self.spec)
        del(self.lin_stage)
        del(self.pos_list)

    


sys.path.append(os.path.abspath(r"C:\PythonCoding\refractive-index-database"))

def browse_location():
    filename = tkinter.filedialog.askdirectory()
    dirName_var.set(filename)
    print("set filename")

def confirm_beam_blocked(hardware, sys_vars):
    blocked.grid_forget()
    hardware.set_dark_spectrum()
    print("unblock")
    Text_block.delete('1.0', tkinter.END)
    Text_block.insert(tkinter.END, "Please unblock the beam! \nClick confirm when your are done.")
    unblocked.grid(row=8, column=1)

def confirm_beam_unblocked(hardware, sys_vars):
    # hide buttons and features from before 
    attn_label.grid_forget()
    Text_block.grid_forget()
    unblocked.grid_forget()
    # Init data array
    data = np.zeros(hardware.dataset_dim())
    
    # make progress bar appear
    progress_label.grid(row=7, column=0)
    progress.grid(row=7, column=1)
    # open live plot
    fig = plt.figure()
    ax =  plt.subplot(111)
    canvas = FigureCanvasTkAgg(fig, 
                               master = root)
    idxLow, idxHigh = hardware.get_low_and_high_wavelenghth_inds()
    im = ax.imshow(data[:,idxLow:idxHigh],aspect = 'auto',origin='lower',extent=[500,800,hardware.pos_list[0],hardware.pos_list[-1]])
    im.set_clim(vmin=0, vmax=65000)
    
    timeAxis = np.zeros(len(hardware.pos_list))
    canvas.get_tk_widget().grid(row=8,column=1,columnspan=10,rowspan=10)   
    toolbarFrame = tkinter.Frame(master=root)
    toolbarFrame.grid(row=17,column=1)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    startT = time.time()    
    for i in range(len(hardware.pos_list)):
        hardware.set_linstage_dist(hardware.pos_list[i])
        time.sleep(0.05)
        # averaging over 10 spectra
        for j in range(sys_vars["avgFactor"]):
            data[i,:] += hardware.measure_spectrum()
        data[i,:] /= sys_vars["avgFactor"]
        #data[i,:] = spec.get_spectrum()
        timeAxis[i]=time.time()
        print(f" Progress {100*(i/len(hardware.pos_list))}% " )
        progress["value"] = 100*(i/len(hardware.pos_list))
        im.set_data(data[:,idxLow:idxHigh])
        canvas.draw()        
        fig.canvas.flush_events()

    sys_vars["canvas"] = canvas
    sys_vars["toolbarFrame"] =toolbarFrame
    progress["value"] = 100

    IntData = np.sum(data[:,idxLow:idxHigh],axis=1)
    IntDataSmoothed = gaussian_filter1d(IntData, sigma=6)
    peak =np.abs(IntDataSmoothed).argmax()
    hardware.set_linstage_dist(hardware.pos_list[peak])

    plt.savefig(sys_vars["dirPath"] + "data.png")
    pd.DataFrame(data).to_csv(sys_vars["dirPath"] +"data.csv")
    pd.DataFrame(hardware.wavelength).to_csv(sys_vars["dirPath"] +"wavelength.csv")
    pd.DataFrame(hardware.pos_list).to_csv(sys_vars["dirPath"] +"posList.csv")
    pd.DataFrame(timeAxis).to_csv(sys_vars["dirPath"] +"timeAxis.csv")
    pd.DataFrame(hardware.darkspectrum).to_csv(sys_vars["dirPath"] +"darkspectrum.csv")

    shutil.copy(__file__, sys_vars["dirPath"]+'measurementScript.py') 
    repeat.grid(row=18, column=0)
    dscan.grid(row=18, column=1)
    
        # fig1, ax1 = plt.subplots(1,1)
        
        # ax1.plot(wavelength[idxLow:idxHigh],np.sum(data[:,idxLow:idxHigh]),axis=0)
        # ax1.set_ylim([-200,200])
        # secax1 = ax1.secondary_xaxis('top', functions =(NIR2MIR,MIR2NIR))
        # secax.set_xlabel('Wavelength MIR (nm)')
        
    return

def reset_all(hardware, sys_vars):
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
    #del()
    ## TODO: add del


def get_inputs(hardware, sys_vars):
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

    print("entered variables:") 
    print("dir:", dirPath)
    print("intTime:", intTime, ", avgFactor:", avgFactor)
    
    hardware.set_spec_int_time(intTime)
    start_btn.grid_forget()
    return min_insert, max_insert

    

def on_gui_close(hardware):
    del(hardware)

def main(hardware, sys_vars):

    # Initialize the controller with the correct device name
    controller = elliptec.Controller('COM4')  # Adjust device name if needed
    hardware.initialize_spec_and_lin_stage(ava.spectrometer(), elliptec.Linear(controller))
    
    min_in, max_in = get_inputs(hardware, sys_vars) # this also sets the spectrometer integration time
    hardware.set_position_list(min_in, max_in)

    # Get wavelength axis
    hardware.get_spec_wavelength()

    # Home the linear stage before usage
    hardware.home_lin_stage()
    
    attn_label.grid(row=6)
    Text_block.grid(row=7, column=1)
    Text_block.insert(tkinter.END, "Please block the beam! \nClick confirm when your are done.")
    blocked.grid(row=8,column=1)

def load_fundamentalSpectraLight(fname):
    """Load a fundamental spectrum from the csv files."""
    # skiprows simply ignores the header
    wavelength, intensity = np.loadtxt(fname, skiprows=4, delimiter="\t", unpack=True)
    return pypret.Spectrum(wavelength * 1e-9, intensity)

def load_dscan_data(fname,wlim,sigma=3):
    """Load a single frog file in .frg format into a MeasurementTrace object.
    
    Assumes that the trace is equidistantly spaced in terms of the wavelength.
    # """
    # # load the whole file and store all values in a 1d-array 
    # with open(fname, "rt") as f:
    #     file = f.read().replace("\r\n", " ")
    # raw = np.fromstring(file, sep=" ")
    
    # M = int(raw[0])  # number of delay steps
    # N = int(raw[1])  # number of wavelength/frequency points
    # # scale to basic units
    # dtau = 1e-15 * raw[2]  # delay discretization 
    # dwl = 1e-9 * raw[3]    # wavelength discretization
    # wl0 = 1e-9 * raw[4]    # central wavelength
    # data = raw[5:].reshape(M, N)  # load in the data

    # # create the specific wavelength/delay grid
    # wl = (np.arange(N) - (N // 2)) * dwl + wl0    
    # tau = (np.arange(M) - (M // 2)) * dtau
    data = pd.read_csv(fname+"/data.csv", sep=",", header=None)
    dataOut = data.to_numpy()
    dataOut = dataOut[1:,1:]
    print(dataOut.shape)
    wl = pd.read_csv(fname+"/wavelength.csv", sep=",", header=None)
    wl = wl.to_numpy()
    wl = wl[1:,1]*1e-9

    IntData = np.sum(dataOut, axis=1)
    # smooth the intData
    IntData_smooth = gaussian_filter1d(IntData, sigma=5)
    Data_smooth = gaussian_filter(dataOut, sigma=sigma)
    
    # find max intensity
    maxInt = np.abs(IntData_smooth).argmax()

    print(wl.size)
    iStart = np.argmin(np.abs(wl-wlim[0]))
    print(iStart)
    iStop = np.argmin(np.abs(wl-wlim[1]))
    print(iStop)

    tau = pd.read_csv(fname+"/posList.csv", sep=",", header=None)
    tau = tau.to_numpy()
    tau = tau[1:,1]*1e-3
    #tau = tau - tau[len(tau)//2]
    # scaling to wedge insertion
    tau = tau * np.arctan(4.5/50)
    tau -= tau[maxInt] +0.0005e-3
    # plt.plot(tau,IntData_smooth)
    # plt.show()
    # scalingFactorArray = IntData_smooth/IntData
    # dataOut = dataOut * scalingFactorArray[:,None]
    
    print(tau.size)

    # convert data (M, N) array into list of M Spectrum instances
    spectra = []
    for i in range(len(tau)):
        # intensity = dataOut[i,:]
        intensity= Data_smooth[i,:]
        tempSpec = pypret.Spectrum(wl[iStart:iStop], intensity[iStart:iStop])
        #tempSpec.remove_dark_signal((450e-9,500e-9))
        spectra.append(tempSpec)
    # create MeasurementTrace with correct method specifier
    trace = pypret.MeasurementTrace(tau, spectra, method="dscan", process="thg")
    return trace
    
def do_dscan_retrieval(hardware, sys_vars):
    progress.grid_forget()
    sys_vars["canvas"].get_tk_widget().grid_forget()
    sys_vars["toolbarFrame"] .grid_forget()
    repeat.grid_forget()
    dscan.grid_forget()
    progress_label.grid_forget()

    measurement = load_dscan_data(sys_vars["dirPath"], (550e-9,850e-9), sigma=3) ## TODO: replace hardcoded value with input
    fundamental = load_fundamentalSpectraLight("sp_2130_04.dat") ## TODO: why hard coded 

    specfile =r"sp_2130_04.dat"
    wavelength, intensity = np.loadtxt(specfile, skiprows=4, delimiter="\t", unpack=True)
    ### Not sure what this is for and I'm a little worried about the rescaling of wavelength 
    # wavelength= wavelength*1e-9
    # y = intensity * wavelength * wavelength
    # y[y < 0.0] = 0.0
    # y = np.sqrt(y + 0.0j)
    # y /= y.max()
    # # calculate angular frequencies
    # w = pypret.spectrum.wl2om(wavelength)

    # # create pulse parameters from the measured spectrum
    # # choose center wavelength as the mean of the intensity
    # w0 = pypret.lib.mean(w, pypret.lib.abs2(y))
    # wl0 = pypret.spectrum.wl2om(w0)
    # # choose simulation grid that encompasses the measured spectrum
    # dw = abs(np.mean(np.diff(w)))
    # N= 768
    # ft = FourierTransform(N, dw=dw, w0=w.min() - w0)
    # pulse2 = pypret.Pulse(ft, w0, unit="om")
    # # interpolate
    # smooth=4
    # s = interp1d(wavelength, y, kind="linear", bounds_error=False, fill_value=(y[0], y[-1]))(pulse2.wl)

    # m = (pulse2.wl <= wavelength[0])
    # s[m] = s[m] * np.exp(-(pulse2.wl[m] - wavelength[0])**2 / (2 * smooth**2))
    # m = (pulse2.wl >= wavelength[-1])
    # s[m] = s[m] * np.exp(-(pulse2.wl[m] - wavelength[-1]) ** 2 / (2 * smooth** 2))
    # pulse2.spectrum = s
    pypret.save((fundamental, measurement), "dscan_ps.h5")
    ## TODO: we need to incorporate this into the GUI framework. Not sure how to.
    fig = plt.figure()
    ax =  plt.subplot(111)
    canvas1 = FigureCanvasTkAgg(fig, 
                               master = root)
    canvas1.get_tk_widget().grid(row=8,column=1,columnspan=10,rowspan=10)   
    toolbarFrame = tkinter.Frame(master=root)
    toolbarFrame.grid(row=17,column=1)
    toolbar = NavigationToolbar2Tk(canvas1, toolbarFrame)
    pypret.graphics.TracePlot(measurement, marginals=True, cmap="nipy_spectral",
                          title="raw trace")
    canvas1.draw()        
    fig.canvas.flush_events()

    measurement.limit_parameter((None, 2.2e-3))  
    measurement.remove_dark_signal((550e-9,570e-9))
    pypret.graphics.TracePlot(measurement, marginals=True, cmap="nipy_spectral",
                            title="removed < Nothing")
    canvas1.draw()        
    fig.canvas.flush_events()
    pulseMeas = measurement.grid(process="thg") 
    pypret.graphics.TracePlot(measurement, marginals=True, cmap="nipy_spectral", title="interpolated")
    canvas1.draw()        
    fig.canvas.flush_events()

    fundamental.normalize()
    #fundamental.show(title="raw fundamental")

    # remove remaining CW
    fundamental.limit((None, None), (None, None))  # limits estimated from the plot
    #fundamental.show(title="CW removed")

    spectrum = fundamental.copy()  # make copy, keep the original Spectrum instance
    spectrum.smooth(50e-9)  # perform a smoothing with 0.5nm bandwidth
    spectrum.bandpass((1.2e-6, 2.52e-6), dx=10e-9)  # apodize around the edges of the spectrum
    #spectrum.show(title="smoothed and apodized")

    # here we use the Pulse instance created when gridding the FROG trace
    tl_pulse = spectrum.create_pulse(pulseMeas.copy())

    # plot the TL pulse with 4 x interpolation
    pypret.graphics.PulsePlot(tl_pulse, oversampling=4, suptitle="TL pulse")
    print("transform-limited pulse duration (FWHM): %.1f fs" % (tl_pulse.fwhm(dt=tl_pulse.dt / 100) * 1e15))

    pnps = pypret.PNPS(pulseMeas, "dscan", "thg", material= pypret.material.ZNSE)
    objective = pypret.LeastSquaresObjective(pnps, measurement)

    trace = objective.trace(tl_pulse.spectrum)
    pypret.graphics.CompareTraceMeshPlot(
        measurement, trace, suptitle="comparison with TL pulse", marginals=True
    )
    plt.show()

    # perform the retrieval, uses COPRA by default. Alternative is method="l-bfgs-b"
    # this may take quite a while to execute (~1min), as the FROG trace is very large and retrieval is thus slow
    fig2 = plt.figure()
    axs =  plt.subplots(3, 1)
    canvas2 = FigureCanvasTkAgg(fig2, 
                               master = root)
    canvas2.get_tk_widget().grid(row=8,column=11,columnspan=10,rowspan=10)   
    toolbarFrame2 = tkinter.Frame(master=root)
    toolbarFrame2.grid(row=17,column=11)
    toolbar2 = NavigationToolbar2Tk(canvas2, toolbarFrame2)
    res = pypret.retrieve(
        pnps,
        measurement,
        repetitions= 16,  # retrieve 16 times from random initial guesses
        verbose=True,  # print out status updates
        maxiter=300,  # perform 300 COPRA iterations
        measured_spectrum=fundamental,  # compare everything to the measured fundamental spectrum
        method = "copra",
        guess = tl_pulse.spectrum,  # use the TL pulse as initial guess
    )

    # save the full results for later usage (file is quite large)
    res.save("retrieved_dscan_ps.hdf5")

    # %% obtain and plot the mean result
    mean = res.mean()
    rpulse = mean.pulse
    mean.show()
    canvas2.draw()        
    fig2.canvas.flush_events()
    print(
        "Retrieved pulse duration (FWHM): %.1f fs"
        % (rpulse.fwhm(dt=tl_pulse.dt / 100) * 1e15)
    )
    repeat.grid(row=18, column=0)



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
    tkinter.Label(root, text='Integration time (ms)').grid(row=1)
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

    ## TODO: add a spectrometer window to see current spectrum 

    # initialize hardware
    
    hardware = HardWare()

    start_btn=tkinter.Button(root,text = 'Start acquisition', command = lambda: main(hardware, sys_vars))
    start_btn.grid(row=5,column=1)
    
    blocked = tkinter.Button(root, text = 'Confirm', command =  lambda: confirm_beam_blocked(hardware, sys_vars))
    attn_label = tkinter.Label(root, text = "Attention!", font=('Helvetica',12, 'bold'))
    Text_block = tkinter.Text(root, height = 5, width = 45)
    unblocked =tkinter.Button(root,text = 'Confirm', command =  lambda: confirm_beam_unblocked(hardware, sys_vars))
    repeat = tkinter.Button(root,text = 'New acquisition', command =  lambda:  reset_all(hardware, sys_vars))
    dscan = tkinter.Button(root, text = 'Phase retrival', command =  lambda:  do_dscan_retrieval(hardware, sys_vars))
    progress = ttk.Progressbar(root, length=100, mode='determinate')
    progress_label = tkinter.Label(root, text='Progress')

    ## initialize default insertion width
    min_insert_var.set(5)
    max_insert_var.set(40)
    intTime_var.set(10)
    avgFactor_var.set(2)
    root.protocol("WM_DELETE_WINDOW", on_gui_close)
    
    tkinter.mainloop()


# #FXIME make closing 
#     del(sys_vars["ls"])
#     del(sys_vars["spec"])
  