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
#import numpy as np
import SysRespLIB as SRL

#RETRIEVE FILTER REGION REFERENCE SPECTRA######################################
REDPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
REDPlotParams.loadplotparams(drive,"PicklesRED","TBD")

REDA_FilesList=SRL.measurement_list(REDPlotParams.DataFile)
REDA_FilesList.load_select_data("A_Stars")
MeanREDA=SRL.SpectrumAggregation("f:",REDA_FilesList)
MeanREDA.ComputeAverageandStats()

REDF_FilesList=SRL.measurement_list(REDPlotParams.DataFile)
REDF_FilesList.load_select_data("F_Stars")
MeanREDF=SRL.SpectrumAggregation("f:",REDF_FilesList)
MeanREDF.ComputeAverageandStats()

REDG_FilesList=SRL.measurement_list(REDPlotParams.DataFile)
REDG_FilesList.load_select_data("G_Stars")
MeanREDG=SRL.SpectrumAggregation("f:",REDG_FilesList)
MeanREDG.ComputeAverageandStats()

REDK_FilesList=SRL.measurement_list(REDPlotParams.DataFile)
REDK_FilesList.load_select_data("K_Stars")
MeanREDK=SRL.SpectrumAggregation("f:",REDK_FilesList)
MeanREDK.ComputeAverageandStats()

REDM_FilesList=SRL.measurement_list(REDPlotParams.DataFile)
REDM_FilesList.load_select_data("M_Stars")
MeanREDM=SRL.SpectrumAggregation("f:",REDM_FilesList)
MeanREDM.ComputeAverageandStats()

#MAKE REFERENCE PLOT#####################################################
REDPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(MeanREDF.MeanSpec,1.33,'m','F Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDG.MeanSpec,1.12,'b','G Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDK.MeanSpec,0.88,'g','K Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDM.MeanSpec,0.42,'r','M Stars')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Nebular Diagnostics/"
pylab.savefig(path+'Stellar Types for RED Spectra.png',dpi=300)

GRNPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
GRNPlotParams.loadplotparams(drive,"PicklesGRN","TBD")
GRNPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(MeanREDF.MeanSpec,1.0,'m','F Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDG.MeanSpec,1.0,'b','G Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDK.MeanSpec,1.0,'g','K Stars')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Nebular Diagnostics/"
pylab.savefig(path+'Stellar Types for GRN Spectra.png',dpi=300)

BLUPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
BLUPlotParams.loadplotparams(drive,"PicklesBLU","TBD")
BLUPlotParams.Setup_Plot()
#tmp=SRL.Draw_with_Conf_Level(MeanREDA.MeanSpec,0.72,'k','A Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDF.MeanSpec,0.72,'m','F Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDG.MeanSpec,1.0,'b','G Stars')
tmp=SRL.Draw_with_Conf_Level(MeanREDK.MeanSpec,1.33,'g','K Stars')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Nebular Diagnostics/"
pylab.savefig(path+'Stellar Types for BLU Spectra.png',dpi=300)

