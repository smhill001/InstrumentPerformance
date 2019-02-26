# -*- coding: utf-8 -*-
"""
Updated on Tue Mar 07 07:22:23 2017

NAME:       SystemResponse.py

PURPOSE:    Analyzes response of 1260mm-656HIA-ST2000 system.

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

ManComp7Data=SRL.manufacturer_Comp7_data(path+"baader_h-alpha.txt")
ManComp7Data.load_records()
ManComp7Data.uniform_wave_grid()

ManBaaderData=SRL.manufacturer_Comp7_data(path+"h_alpha_ccd_7nm_kurve.txt")
ManBaaderData.load_records()
ManBaaderData.uniform_wave_grid()

HIAPlotParams=SRL.SysResp_plot_params("FluxCalPlotConfig.txt")
HIAPlotParams.loadplotparams(drive,"656HIA","TBD")
HIA656_ObsList=SRL.measurement_list(HIAPlotParams.DataFile)
#HIA656_ObsList.load_select_data("1260mm200lpm")
#Mean200linespermm1260mm=SRL.SpectrumAggregation("f:",CLR550_ObsList)
#Mean200linespermm1260mm.ComputeAverageandStats()
#tmp=SRL.Compute_EWs(path,"1260mm200lpm-550CLR-EW",Mean200linespermm1260mm.MeanSpec,1.0189)
#np.savetxt(path+'SystemResponseCLR-1260mm200lpm.txt',Mean200linespermm1260mm.MeanSpec,delimiter=" ",
#           fmt="%10.3F %10.7F %10.7F %10.7F")

#MAKE ST2000 MEASURED RESPONSE PLOT#####################################################
"""
HIAPlotParams.Setup_Plot()
tmp=SRL.Draw_with_Conf_Level(Mean200linespermm1260mm.MeanSpec,1.0189,'C0','1260mm200lpm')
tmp=SRL.Draw_with_Conf_Level(Mean100linespermm1260mm.MeanSpec,1.0,'b','1260mm100lpm')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'Response1260mmCLR-Measured.png',dpi=300)
"""

CLRPlotParams.Setup_Plot()

pl.plot(ManComp7Data.wavegrid,ManCamData.griddata,
           label='Company7',linewidth=1.0,color='C7')

pl.plot(ManBaaderData.wavegrid,Atmosphere_Data.griddata,
           label='Baader',linewidth=1.0,color='b')


pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'Response1260mmHIA-Modeled.png',dpi=300)

