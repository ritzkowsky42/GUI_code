from avaspec import *
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import time as time


class spectrometer:

    def __init__(self):
                
        #(-1 Eth and USB, 0 USB, 256 Eth)
        self.ret = AVS_Init(-1)
        self.ret = AVS_UpdateUSBDevices()
        self.ret = AVS_GetList()
        self.handle = AVS_Activate(self.ret[0])

        config = DeviceConfigType
        ret = AVS_GetParameter(self.handle)
        
        pixels = AVS_GetNumPixels(self.handle)

        lamb = AVS_GetLambda(self.handle)
        self.wavelengths = []
        for pix in range(pixels):
            self.wavelengths.append(lamb[pix])
            pixels = AVS_GetNumPixels(self.handle)

        ret = AVS_UseHighResAdc(self.handle, True)

        self.measconfig = MeasConfigType()
        self.measconfig.m_StartPixel = 0
        self.measconfig.m_StopPixel = pixels - 1
        self.measconfig.m_IntegrationTime = 50 #in milliseconds
        self.measconfig.m_IntegrationDelay = 0 #in FPGA clock cycles
        self.measconfig.m_NrAverages = 1
        self.measconfig.m_CorDynDark_m_Enable = 0  # nesting of types does NOT work!!
        self.measconfig.m_CorDynDark_m_ForgetPercentage = 100
        self.measconfig.m_Smoothing_m_SmoothPix = 0
        self.measconfig.m_Smoothing_m_SmoothModel = 0
        self.measconfig.m_SaturationDetection = 0
        self.measconfig.m_Trigger_m_Mode = 0
        self.measconfig.m_Trigger_m_Source = 0
        self.measconfig.m_Trigger_m_SourceType = 0
        self.measconfig.m_Control_m_StrobeControl = 0
        self.measconfig.m_Control_m_LaserDelay = 0
        self.measconfig.m_Control_m_LaserWidth = 0
        self.measconfig.m_Control_m_LaserWaveLength = 0.0
        self.measconfig.m_Control_m_StoreToRam = 0

        self.ret = AVS_PrepareMeasure(self.handle, self.measconfig)

    def get_wavelengths(self):
       
        return np.array(self.wavelengths)

    def get_spectrum(self):
        self.scans = 1
        self.ret = AVS_Measure(self.handle, 0, self.scans)
        self.dataready = False
        while not self.dataready:
            self.dataready = AVS_PollScan(self.handle)
            time.sleep(self.measconfig.m_IntegrationTime/1000)
                
        self.spectrum = []        
        self.timestamp, self.scopedata = AVS_GetScopeData(self.handle)
        for i,pix in enumerate(self.wavelengths):
            self.spectrum.append(self.scopedata[i])

        return self.spectrum
    
    def set_integration_time(self, time):
        self.measconfig.m_IntegrationTime = time
        self.ret = AVS_PrepareMeasure(self.handle, self.measconfig)

    def __del__(self):
        AVS_StopMeasure(self.handle)
        AVS_Deactivate(self.handle)
        AVS_Done()

def main():
    # Testing the spectrometer
    spec = spectrometer()
    spec.set_integration_time(50)

    wavelength = spec.get_wavelengths()
    idxLow = np.abs(wavelength-200).argmin()
    idxHigh = np.abs(wavelength-1100).argmin()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    
    line, = ax.plot(wavelength[idxLow:idxHigh],spec.get_spectrum()[idxLow:idxHigh])
    ax.set_ylim([0,65000])
    ax.set_xlabel("Wavelength (nm)")



    try:
        while True:
            line.set_ydata(spec.get_spectrum()[idxLow:idxHigh])
            fig.canvas.draw()        
            fig.canvas.flush_events()
    
    except KeyboardInterrupt:
        pass
    # line.set_ydata(spec.get_spectrum()[idxLow:idxHigh])
    # fig.canvas.draw()        
    # fig.canvas.flush_events()
    



    
if __name__ == '__main__':
    main()

