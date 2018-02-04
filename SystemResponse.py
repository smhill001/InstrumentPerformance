# -*- coding: utf-8 -*-
"""
Updated on Tue Mar 07 07:22:23 2017

NAME:       FluxCalibration_V001.py

PURPOSE:    Produces the relative response of a given filter or line normalized
            to Green (GRN550?)

PURPOSE: To create data and plots about the response of filters relative to
         the peak response in the 550GRN filter.

@author: steven.hill
FUNCTION: N/A

This module reads the output text files created by FilterTransmissionfromFITS
for given filters. It then averages, normalizes and plots the resulting 
response curves. Additionally, a text file is produced characterizing each
filter.

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
    PlotParameters - List of plotting configuration info related to the 
             filter-target combination not the specific observation
    ProcessingConfigFile.txt - Two parameters that control processing, plus
             two others that are currently unused. This permits code reuse.
    FileList - List of individual FITS files for processing
    ObsBands - List of spectral bands for which equiv. widths will be calculated
    FITS Images - 2D image files with necessary metadata for flux and wavelength
             calibration            
             
OUTPUTS:
    1D Spectrum - txt file
    1D Spectrum - png file
    Equivalent widths - txt file
****
 
"""
import sys
drive='f:'
sys.path.append(drive+'\\Astronomy\Python Play')
sys.path.append(drive+'\\Astronomy\Python Play\Techniques Library')
sys.path.append(drive+'\\Astronomy\Python Play\Galaxies')

import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate
from copy import deepcopy
import GalaxyLIB as GL
import SysRespLIB as SRL

###############################################################################
path=drive+"/Astronomy/Projects/Stars/Vega/Spectral Data/1D Spectra/"
# Read and reshape spectral data files    
Vega20130921UT = scipy.fromfile(file=path+"VegaResponse20130921UT.txt", dtype=float, count=-1, sep='\t')    
Vega20130921UT=scipy.reshape(Vega20130921UT,[Vega20130921UT.size/2,2])

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("1260mm200lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean200linespermm1260mm=MeanCLR

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("135mm100lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean100linespermm135mm=MeanCLR

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("1260mm100lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean100linespermm1260mm=MeanCLR
###############################################################################
PlotParams=SRL.SysResp_plot_params(drive,"685NIR","TBD")
NIR685_ObsList=SRL.measurement_list(PlotParams.DataFile)
NIR685_ObsList.load_all_data()
MeanNIR=SRL.Average_Spectrum("f:",NIR685_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"650RED","TBD")
RED650_ObsList=SRL.measurement_list(PlotParams.DataFile)
RED650_ObsList.load_all_data()
MeanRED=SRL.Average_Spectrum("f:",RED650_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"550GRN","TBD")
GRN550_ObsList=SRL.measurement_list(PlotParams.DataFile)
GRN550_ObsList.load_all_data()
MeanGRN=SRL.Average_Spectrum("f:",GRN550_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"450BLU","TBD")
BLU450_ObsList=SRL.measurement_list(PlotParams.DataFile)
BLU450_ObsList.load_all_data()
MeanBLU=SRL.Average_Spectrum("f:",BLU450_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"380NUV","TBD")
NUV380_ObsList=SRL.measurement_list(PlotParams.DataFile)
NUV380_ObsList.load_all_data()
print "*********NUV380_ObsList.FileList",NUV380_ObsList.FileList
MeanNUV=SRL.Average_Spectrum("f:",NUV380_ObsList)

###############################################################################
pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0,x1,xtks=350,1050,15
y0,y1=1e-3,1.2e0

pl.xlim(x0,x1)
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))

pl.ylim(y0,y1)
pl.yscale('log')

pl.grid()
pl.tick_params(axis='both', which='major', labelsize=7)

pl.ylabel(r"$Normalized$ $Response$",fontsize=7)
pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
pl.title("Normalized Response",fontsize=9)

#pl.plot(Vega20130921UT[:,0]/10.,Vega20130921UT[:,1],label='20130921UT',linewidth=1)

tmp=SRL.Draw_with_Conf_Level(Mean200linespermm1260mm,1.0189,'k','1260mm200lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm135mm,1.0244,'0.5','135mm100lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm1260mm,1.0,'C5','1260mm100lpm')
tmp=SRL.Draw_with_Conf_Level(MeanNIR,0.40,'C3','zNIR')
#tmp=SRL.Draw_with_Conf_Level(MeanNIR,1.0,'C3','zNIR')
tmp=SRL.Draw_with_Conf_Level(MeanRED,0.690,'r','zRED')
#tmp=SRL.Draw_with_Conf_Level(MeanRED,1.0,'r','zRED')
tmp=SRL.Draw_with_Conf_Level(MeanGRN,0.980,'g','zGRN')
tmp=SRL.Draw_with_Conf_Level(MeanBLU,0.875,'b','zBLU')
tmp=SRL.Draw_with_Conf_Level(MeanNUV,0.016,'m','zNUV')

#TransRED=MeanRED[430:2325,:]
#TransRED[:,1]=(0.58*TransRED[:,1]/1.08)/Mean200linespermm1260mm[:,1]
#tmp=SRL.Draw_with_Conf_Level(TransRED,1.0,'r','TRED')

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Flux Calibration/"

pylab.savefig(path+'SystemResponse.png',dpi=300)

np.savetxt(path+'SystemResponseCLR.txt',Mean200linespermm1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")



tmp=SRL.Compute_EWs(path,Mean200linespermm1260mm)
