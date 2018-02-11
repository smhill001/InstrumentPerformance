# -*- coding: utf-8 -*-
"""
Updated on Tue Mar 07 07:22:23 2017

NAME:       SystemResponse.py

PURPOSE:    Analyzes response of optics-camera-filter combinations.

@author: steven.hill
Update: 2/11/2018

This script plots the responses of optics-camera-filter combinations as well
as individual filter transmissions. It also computes the equivalent widths and
mean sensitivities of filters. Input data are text files of spectra. Program 
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
MeanNIR=MeanNIRtmp.MeanSpec[430:2325,:]

REDPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
REDPlotParams.loadplotparams(drive,"650RED","TBD")
RED650_ObsList=SRL.measurement_list(REDPlotParams.DataFile)
RED650_ObsList.load_all_data()
MeanREDtmp=SRL.SpectrumAggregation("f:",RED650_ObsList)
MeanREDtmp.ComputeAverageandStats()
MeanRED=MeanREDtmp.MeanSpec[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-650RED-EW",MeanRED,0.69/1.02424)
np.savetxt(path+'SystemResponseRED-135mm100lpm.txt',MeanRED,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

GRNPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
GRNPlotParams.loadplotparams(drive,"550GRN","TBD")
GRN550_ObsList=SRL.measurement_list(GRNPlotParams.DataFile)
GRN550_ObsList.load_all_data()
MeanGRNtmp=SRL.SpectrumAggregation("f:",GRN550_ObsList)
MeanGRNtmp.ComputeAverageandStats()
MeanGRN=MeanGRNtmp.MeanSpec[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-550GRN-EW",MeanGRN,0.98/1.0224)
np.savetxt(path+'SystemResponseGRN-135mm100lpm.txt',MeanGRN,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

BLUPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
BLUPlotParams.loadplotparams(drive,"450BLU","TBD")
BLU450_ObsList=SRL.measurement_list(BLUPlotParams.DataFile)
BLU450_ObsList.load_all_data()
MeanBLUtmp=SRL.SpectrumAggregation("f:",BLU450_ObsList)
MeanBLUtmp.ComputeAverageandStats()
MeanBLU=MeanBLUtmp.MeanSpec[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-550BLU-EW",MeanBLU,0.875/1.0224)
np.savetxt(path+'SystemResponseBLU-135mm100lpm.txt',MeanBLU,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

NUVPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
NUVPlotParams.loadplotparams(drive,"380NUV","TBD")
NUV380_ObsList=SRL.measurement_list(NUVPlotParams.DataFile)
NUV380_ObsList.load_all_data()
MeanNUVtmp=SRL.SpectrumAggregation("f:",NUV380_ObsList)
MeanNUVtmp.ComputeAverageandStats()
MeanNUV=MeanNUVtmp.MeanSpec[430:2325,:]

#MAKE ST2000 RESPONSE PLOT#####################################################
CLRPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(Mean200linespermm1260mm.MeanSpec,1.0189,'C0','1260mm200lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm1260mm.MeanSpec,1.0,'b','1260mm100lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm135mm.MeanSpec,1.0244,'C3','135mm100lpm')
tmp=SRL.Draw_with_Conf_Level(Mean200linespermm135mm.MeanSpec,1.0122,'r','135mm200lpm')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'PanchromaticSystemResponse.png',dpi=300)

#MAKE 135MM-ST2000 PLOT WITH FILTERS###########################################
CLRPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm135mm.MeanSpec,1.0244,'0.5','135mm100lpm')
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
CLRPlotParams.Setup_Plot()
#Should be able to make a transmission computing function with two inputs
print "MeanNIR.shape:",MeanNIR.shape

TransNIR=SRL.Compute_Transmission(MeanNIR,Mean100linespermm135mm.MeanSpec)
tmp=SRL.Draw_with_Conf_Level(TransNIR,0.40/1.0244,'C3','NIR')
TransRED=SRL.Compute_Transmission(MeanRED,Mean100linespermm135mm.MeanSpec)
tmp=SRL.Draw_with_Conf_Level(TransRED,0.69/1.0244,'r','RED')
TransGRN=SRL.Compute_Transmission(MeanGRN,Mean100linespermm135mm.MeanSpec)
tmp=SRL.Draw_with_Conf_Level(TransGRN,0.98/1.0244,'g','GRN')
TransBLU=SRL.Compute_Transmission(MeanBLU,Mean100linespermm135mm.MeanSpec)
tmp=SRL.Draw_with_Conf_Level(TransBLU,0.875/1.0244,'b','BLU')
TransNUV=SRL.Compute_Transmission(MeanNUV,Mean100linespermm135mm.MeanSpec)
tmp=SRL.Draw_with_Conf_Level(TransNUV,0.016/01.244,'m','NUV')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'FilterTransmission.png',dpi=300)

REDPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(TransRED,0.69/1.0244,'r','RED Transmission')
tmp=SRL.Draw_with_Conf_Level(MeanRED,0.69,'C3','RED Response')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'REDTransmission.png',dpi=300)

GRNPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(TransGRN,0.98/1.0244,'g','GRN Transmission')
tmp=SRL.Draw_with_Conf_Level(MeanGRN,0.98,'C2','GRN Response')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'GRNTransmission.png',dpi=300)

BLUPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(TransBLU,0.875/1.0244,'b','BLU Transmission')
tmp=SRL.Draw_with_Conf_Level(MeanBLU,0.875,'C0','BLU Response')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'BLUTransmission.png',dpi=300)