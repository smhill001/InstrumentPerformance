# -*- coding: utf-8 -*-
"""
Updated on Tue Mar 07 07:22:23 2017

NAME:       SystemResponse.py

PURPOSE:    Analyzes response of 1260mm-CLR-ST2000 system.

@author: steven.hill
Update: 4/8/2018

This script plots the responses of the 1260mm-CLR-ST2000 combinations. 
It also computes the equivalent widths and mean sensitivities of filter 
windows. Input data are text files of spectra. Program 
control is by a series of configuration files.

****INPUT PARAMETERS (hard coded as of 6/8/16):
    Filter - This string identifies the filter by a alphanumeric index
             that includes the approximate central wavelength, e.g., 380NUV
    Target - This string identifies the target, e.g., "Jupiter", "Vega", etc.
             This parameter allows the lookup of the target type, e.g.,
             "Planet", "Star", etc. that permits construction of directory
             paths.
    DateUT - The UT date of the observation. Combined with Target, this forms
             a unique key to the observation, assuming that most unique parameters
             are invariant over a single observing night.
            
INPUT FILES:
    FluxCalPlotConfig.txt - Master list of target plots, parameters and data.
        It provides pointers to other tables controlling inputs data. See
        the GitHub wiki for the Techniques project for more details on the
        metadata.
    Spectral Files - The actual data plotted are spectral response files. 
        These files can have been produced by one of several spectra-from-FITS
        programs.
             
OUTPUTS:
    Graphics Files - PNG plots are sent to the 
        ../Projects/Techniques/Flux Calibration directory
    Equivalent widths - text files are sent to the 
        ../Projects/Techniques/Flux Calibration directory
"""
import sys
drive='f:'
sys.path.append(drive+'\\Astronomy\Python Play')
sys.path.append(drive+'\\Astronomy\Python Play\Techniques Library')
sys.path.append(drive+'\\Astronomy\Python Play\Galaxies')

import matplotlib.pyplot as pl
import pylab
import numpy as np
import SysRespLIB as SRL

#RETRIEVE ST2000 RESPONSES#####################################################
path=drive+"/Astronomy/Projects/Techniques/Flux Calibration/"
# Read and reshape spectral data files    

ManCamData=SRL.manufacturer_camera_data(path+"CameraResponse-ST2000Data.txt")
ManCamData.load_all_data()
ManCamData.uniform_wave_grid()

C8_StarBright_Data=SRL.manufacturer_Celestrom_data(path+"CameraResponse - StarBright.txt")
C8_StarBright_Data.load_all_data()
C8_StarBright_Data.uniform_wave_grid()

C8_StarBrightXLT_Data=SRL.manufacturer_Celestrom_data(path+"CameraResponse - StarBrightXLT.txt")
C8_StarBrightXLT_Data.load_all_data()
C8_StarBrightXLT_Data.uniform_wave_grid()

Atmosphere_Data=SRL.atmosphere_data(path+"CameraResponse - Atmosphere.txt")
Atmosphere_Data.load_all_data()
Atmosphere_Data.uniform_wave_grid()

CLRPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
CLRPlotParams.loadplotparams(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(CLRPlotParams.DataFile)
CLR550_ObsList.load_select_data("1260mm200lpm")
Mean200linespermm1260mm=SRL.SpectrumAggregation("f:",CLR550_ObsList)
Mean200linespermm1260mm.ComputeAverageandStats()
tmp=SRL.Compute_EWs(path,"1260mm200lpm-550CLR-EW",Mean200linespermm1260mm.MeanSpec,1.0189)
np.savetxt(path+'SystemResponseCLR-1260mm200lpm.txt',Mean200linespermm1260mm.MeanSpec,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

CLR550_ObsList.load_select_data("1260mm100lpm")
Mean100linespermm1260mm=SRL.SpectrumAggregation("f:",CLR550_ObsList)
Mean100linespermm1260mm.ComputeAverageandStats()

CLR550_ObsList.load_select_data("135mm100lpm")
Mean100linespermm135mm=SRL.SpectrumAggregation("f:",CLR550_ObsList)
Mean100linespermm135mm.ComputeAverageandStats()
tmp=SRL.Compute_EWs(path,"135mm100lpm-550CLR-EW",Mean100linespermm135mm.MeanSpec,1.0224)
np.savetxt(path+'SystemResponseCLR-135mm100lpm.txt',Mean100linespermm135mm.MeanSpec,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

CLR550_ObsList.load_select_data("135mm200lpm")
Mean200linespermm135mm=SRL.SpectrumAggregation("f:",CLR550_ObsList)
Mean200linespermm135mm.ComputeAverageandStats()

#RETRIEVE FILTER RESPONSES#####################################################
NIRPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
NIRPlotParams.loadplotparams(drive,"685NIR","TBD")
NIR685_ObsList=SRL.measurement_list(NIRPlotParams.DataFile)
NIR685_ObsList.load_all_data()
MeanNIRtmp=SRL.SpectrumAggregation("f:",NIR685_ObsList)
MeanNIRtmp.ComputeAverageandStats()
MeanNIR135mm=MeanNIRtmp.MeanSpec[430:2325,:]

#MAKE ST2000 MEASURED RESPONSE PLOT#####################################################
CLRPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(Mean200linespermm1260mm.MeanSpec,1.0189,'C0','1260mm200lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm1260mm.MeanSpec,1.0,'b','1260mm100lpm')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'Response1260mmCLR-Measured.png',dpi=300)


CLRPlotParams.Setup_Plot()

#Plot Mfg. QE -> response to photon counts
#pl.scatter(ManCamData.Wavelength,ManCamData.ST2000_Norm,
#           label='ST2000 SBIG QE',s=20,color='b')
#Plot Mfg. response to energy (~QE/wavelength)
EnergyResponse=ManCamData.griddata/ManCamData.wavegrid*450.

#pl.scatter(ManCamData.Wavelength,np.array(ManCamData.ST2000_Norm)/
#           np.array(ManCamData.Wavelength)*450.,
#           label='ST2000 SBIG Ergs',s=20,color='g')

pl.plot(ManCamData.wavegrid,ManCamData.griddata,
           label='ST2000 SBIG Response',linewidth=1.0,color='C7')

pl.plot(Atmosphere_Data.wavegrid,Atmosphere_Data.griddata,
           label='Atmospheric Transmission',linewidth=1.0,color='b')


pl.plot(C8_StarBrightXLT_Data.wavegrid,C8_StarBright_Data.griddata*EnergyResponse,
           label='StarBright Response',linewidth=1.0,color='C1')
pl.plot(C8_StarBrightXLT_Data.wavegrid,C8_StarBrightXLT_Data.griddata*EnergyResponse,
           label='StarBrightXLT Response',linewidth=1.0,color='C1')

#PLOT TOTAL SYSTEM RESPONSE MODELS
pl.plot(C8_StarBright_Data.wavegrid,
        (C8_StarBright_Data.griddata*EnergyResponse*Atmosphere_Data.griddata)/0.52,
        label='StarBright Energy Response',linewidth=1.0,color='g')
pl.plot(C8_StarBright_Data.wavegrid,
        (C8_StarBright_Data.griddata*ManCamData.griddata*Atmosphere_Data.griddata)/0.62,
        label='Total System Response',linewidth=1.0,color='k')
pl.plot(C8_StarBrightXLT_Data.wavegrid,
        (C8_StarBrightXLT_Data.griddata*EnergyResponse*Atmosphere_Data.griddata)/0.6,
        label='XLT Model Response',linewidth=1.0,color='r')

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'Response1260mmCLR-Modeled.png',dpi=300)

