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
import SysRespLIB as SRL

#RETRIEVE ST2000 RESPONSES#####################################################
path=drive+"/Astronomy/Projects/Techniques/Flux Calibration/"
# Read and reshape spectral data files    

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("1260mm200lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean200linespermm1260mm=MeanCLR
tmp=SRL.Compute_EWs(path,"1260mm200lpm-550CLR-EW",Mean200linespermm1260mm)
np.savetxt(path+'SystemResponseCLR-1260mm200lpm.txt',Mean200linespermm1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("1260mm100lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean100linespermm1260mm=MeanCLR

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("135mm100lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean100linespermm135mm=MeanCLR
tmp=SRL.Compute_EWs(path,"135mm100lpm-550CLR-EW",Mean100linespermm135mm)
np.savetxt(path+'SystemResponseCLR-135mm100lpm.txt',Mean100linespermm135mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
CLR550_ObsList.load_select_data("135mm200lpm")
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean200linespermm135mm=MeanCLR

#RETRIEVE FILTER RESPONSES#####################################################
PlotParams=SRL.SysResp_plot_params(drive,"685NIR","TBD")
NIR685_ObsList=SRL.measurement_list(PlotParams.DataFile)
NIR685_ObsList.load_all_data()
MeanNIRtmp=SRL.Average_Spectrum("f:",NIR685_ObsList)
MeanNIR=MeanNIRtmp[430:2325,:]

PlotParams=SRL.SysResp_plot_params(drive,"650RED","TBD")
RED650_ObsList=SRL.measurement_list(PlotParams.DataFile)
RED650_ObsList.load_all_data()
MeanREDtmp=SRL.Average_Spectrum("f:",RED650_ObsList)
MeanRED=MeanREDtmp[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-650RED-EW",MeanRED)
np.savetxt(path+'SystemResponseRED-135mm100lpm.txt',MeanRED,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

PlotParams=SRL.SysResp_plot_params(drive,"550GRN","TBD")
GRN550_ObsList=SRL.measurement_list(PlotParams.DataFile)
GRN550_ObsList.load_all_data()
MeanGRNtmp=SRL.Average_Spectrum("f:",GRN550_ObsList)
MeanGRN=MeanGRNtmp[430:2325,:]

PlotParams=SRL.SysResp_plot_params(drive,"450BLU","TBD")
BLU450_ObsList=SRL.measurement_list(PlotParams.DataFile)
BLU450_ObsList.load_all_data()
MeanBLUtmp=SRL.Average_Spectrum("f:",BLU450_ObsList)
MeanBLU=MeanBLUtmp[430:2325,:]

PlotParams=SRL.SysResp_plot_params(drive,"380NUV","TBD")
NUV380_ObsList=SRL.measurement_list(PlotParams.DataFile)
NUV380_ObsList.load_all_data()
print "*********NUV380_ObsList.FileList",NUV380_ObsList.FileList
MeanNUVtmp=SRL.Average_Spectrum("f:",NUV380_ObsList)
MeanNUV=MeanNUVtmp[430:2325,:]

#MAKE ST2000 RESPONSE PLOT#####################################################
SRL.Setup_Plot("linear")
tmp=SRL.Draw_with_Conf_Level(Mean200linespermm1260mm,1.0189,'C0','1260mm200lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm1260mm,1.0,'b','1260mm100lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm135mm,1.0244,'C3','135mm100lpm')
tmp=SRL.Draw_with_Conf_Level(Mean200linespermm135mm,1.0122,'r','135mm200lpm')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'PanchromaticSystemResponse.png',dpi=300)

#MAKE 135MM-ST2000 PLOT WITH FILTERS###########################################
SRL.Setup_Plot("linear")
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm135mm,1.0244,'0.5','135mm100lpm')
tmp=SRL.Draw_with_Conf_Level(MeanNIR,0.40,'C3','NIR')
tmp=SRL.Draw_with_Conf_Level(MeanRED,0.690,'r','RED')
tmp=SRL.Draw_with_Conf_Level(MeanGRN,0.980,'g','GRN')
tmp=SRL.Draw_with_Conf_Level(MeanBLU,0.875,'b','BLU')
tmp=SRL.Draw_with_Conf_Level(MeanNUV,0.016,'m','NUV')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'FilterRelativeSystemResponse.png',dpi=300)

#MAKE TRANSMISSION PLOT########################################################
SRL.Setup_Plot("linear")
#Should be able to make a transmission computing function with two inputs
TransNIR=SRL.Compute_Transmission(MeanNIR,Mean100linespermm135mm)
tmp=SRL.Draw_with_Conf_Level(TransNIR,0.40/1.0244,'C3','NIR')
TransRED=SRL.Compute_Transmission(MeanRED,Mean100linespermm135mm)
tmp=SRL.Draw_with_Conf_Level(TransRED,0.69/1.0244,'r','RED')
TransGRN=SRL.Compute_Transmission(MeanGRN,Mean100linespermm135mm)
tmp=SRL.Draw_with_Conf_Level(TransGRN,0.98/1.0244,'g','GRN')
TransBLU=SRL.Compute_Transmission(MeanBLU,Mean100linespermm135mm)
tmp=SRL.Draw_with_Conf_Level(TransBLU,0.875/1.0244,'b','BLU')
TransNUV=SRL.Compute_Transmission(MeanNUV,Mean100linespermm135mm)
tmp=SRL.Draw_with_Conf_Level(TransNUV,0.016/01.244,'m','NUV')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'FilterTransmission.png',dpi=300)