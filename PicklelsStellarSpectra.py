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
sys.path.append(drive+'\\Astronomy\Python Play\Utils')

import matplotlib.pyplot as pl
import pylab
#import numpy as np
import PlotUtils as PU
import ConfigFiles as CF
import GeneralSpecUtils as GSU

#RETRIEVE FILTER REGION REFERENCE SPECTRA######################################
REDPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
REDPlotParams.loadplotparams(drive,"PicklesRED","TBD")
#REDPlotParams.loadplotparams(drive,"PicklesREDLines","TBD")

REDO_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDO_FilesList.load_records(MeasTgt="O_Stars")
MeanREDO=GSU.SpectrumAggregation("f:",REDO_FilesList)
MeanREDO.ComputeAverageandStats()

REDB_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDB_FilesList.load_records(MeasTgt="B_Stars")
MeanREDB=GSU.SpectrumAggregation("f:",REDB_FilesList)
MeanREDB.ComputeAverageandStats()

REDA_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDA_FilesList.load_records(MeasTgt="A_Stars")
MeanREDA=GSU.SpectrumAggregation("f:",REDA_FilesList)
MeanREDA.ComputeAverageandStats()

REDF_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDF_FilesList.load_records(MeasTgt="F_Stars")
MeanREDF=GSU.SpectrumAggregation("f:",REDF_FilesList)
MeanREDF.ComputeAverageandStats()

REDG_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDG_FilesList.load_records("G_Stars")
MeanREDG=GSU.SpectrumAggregation("f:",REDG_FilesList)
MeanREDG.ComputeAverageandStats()

REDK_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDK_FilesList.load_records(MeasTgt="K_Stars")
MeanREDK=GSU.SpectrumAggregation("f:",REDK_FilesList)
MeanREDK.ComputeAverageandStats()

REDM_FilesList=CF.measurement_list(REDPlotParams.DataFile)
REDM_FilesList.load_records(MeasTgt="M_Stars")
MeanREDM=GSU.SpectrumAggregation("f:",REDM_FilesList)
MeanREDM.ComputeAverageandStats()

#MAKE REFERENCE PLOT#####################################################
REDPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(MeanREDO.MeanSpec,1.82,'k','O Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDB.MeanSpec,1.66,'C9','B Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDA.MeanSpec,1.48,'C0','A Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDF.MeanSpec,1.33,'m','F Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDG.MeanSpec,1.12,'b','G Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDK.MeanSpec,0.88,'g','K Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDM.MeanSpec,0.42,'r','M Stars',step=True)

pl.plot([654.8,654.8],[0.,2.],'C1',linestyle='--',label="NII 654.8 nm",linewidth=1.0)
pl.plot([656.3,656.3],[0.,2.],'k',linestyle='--',label="HI 656.3 nm",linewidth=1.0)
pl.plot([658.3,658.3],[0.,2.],'C1',linestyle='--',linewidth=1.0)
pl.plot([671.6,671.6],[0.,2.],'C3',linestyle='--',label="SII 671.6 nm",linewidth=1.0)
pl.plot([673.1,673.1],[0.,2.],'C3',linestyle='--',linewidth=1.0)

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Nebular Diagnostics/"
pylab.savefig(path+'Stellar Types for RED Spectra.png',dpi=300)

GRNPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
GRNPlotParams.loadplotparams(drive,"PicklesGRN","TBD")
GRNPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(MeanREDF.MeanSpec,1.0,'m','F Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDG.MeanSpec,1.0,'b','G Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDK.MeanSpec,1.0,'g','K Stars',step=True)

pl.plot([486.0,486.0],[0.,2.],'b',linestyle='--',label="HI 486.0 nm",linewidth=1.0)
pl.plot([496.0,496.0],[0.,2.],'g',linestyle='--',label="OIII 496 nm",linewidth=1.0)
pl.plot([500.7,500.7],[0.,2.],'g',linestyle='--',linewidth=1.0)

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Nebular Diagnostics/"
pylab.savefig(path+'Stellar Types for GRN Spectra.png',dpi=300)

BLUPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
BLUPlotParams.loadplotparams(drive,"PicklesBLU","TBD")
BLUPlotParams.Setup_Plot()
#tmp=SRL.Draw_with_Conf_Level(MeanREDA.MeanSpec,0.72,'k','A Stars')
tmp=PU.Draw_with_Conf_Level(MeanREDF.MeanSpec,0.72,'m','F Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDG.MeanSpec,1.0,'b','G Stars',step=True)
tmp=PU.Draw_with_Conf_Level(MeanREDK.MeanSpec,1.33,'g','K Stars',step=True)

pl.plot([486.0,486.0],[0.,2.],'b',linestyle='--',label="HI 486.0 nm",linewidth=1.0)
pl.plot([496.0,496.0],[0.,2.],'g',linestyle='--',label="OIII 496 nm",linewidth=1.0)
pl.plot([500.7,500.7],[0.,2.],'g',linestyle='--',linewidth=1.0)

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Nebular Diagnostics/"
pylab.savefig(path+'Stellar Types for BLU Spectra.png',dpi=300)

