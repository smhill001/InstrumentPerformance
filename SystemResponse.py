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

#Vega20140902UT = scipy.fromfile(file=path+"VegaResponse20140902UT.txt", dtype=float, count=-1, sep='\t')    
#Vega20140902UT=scipy.reshape(Vega20140902UT,[Vega20140902UT.size/2,2])

#Vega20140916UT = scipy.fromfile(file=path+"VegaResponse20140916UT.txt", dtype=float, count=-1, sep='\t')    
#Vega20140916UT=scipy.reshape(Vega20140916UT,[Vega20140916UT.size/2,2])


PlotParams=SRL.SysResp_plot_params(drive,"550CLR","TBD")
CLR550_ObsList=SRL.measurement_list(PlotParams.DataFile)
MeanCLR=SRL.Average_Spectrum("f:",CLR550_ObsList)
Mean200linespermm1260mm=MeanCLR
print "#####################MeanCLR.shape=",MeanCLR.shape
###############################################################################
PlotParams=SRL.SysResp_plot_params(drive,"685NIR","TBD")
NIR685_ObsList=SRL.measurement_list(PlotParams.DataFile)
MeanNIR=SRL.Average_Spectrum("f:",NIR685_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"650RED","TBD")
RED650_ObsList=SRL.measurement_list(PlotParams.DataFile)
MeanRED=SRL.Average_Spectrum("f:",RED650_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"550GRN","TBD")
GRN550_ObsList=SRL.measurement_list(PlotParams.DataFile)
MeanGRN=SRL.Average_Spectrum("f:",GRN550_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"450BLU","TBD")
BLU450_ObsList=SRL.measurement_list(PlotParams.DataFile)
MeanBLU=SRL.Average_Spectrum("f:",BLU450_ObsList)

PlotParams=SRL.SysResp_plot_params(drive,"380NUV","TBD")
NUV380_ObsList=SRL.measurement_list(PlotParams.DataFile)
MeanNUV=SRL.Average_Spectrum("f:",NUV380_ObsList)

#compute mean ignoring NANs   

#Vega20140916UT_Interp=interpolate.interp1d(Vega20140916UT[:,0],Vega20140916UT[:,1],kind='linear', copy=True,
#                         bounds_error=False, fill_value=0.0)  
#Vega20140916UT_on_20140902UT=Vega20140916UT_Interp(Vega20140902UT[:,0])

#VegaAvg_20140916and19=Vega20140916UT_on_20140919UT
#VegaAvg_20140902and16=(Vega20140916UT_on_20140902UT+Vega20140902UT[:,1])/2.
#Begin plotting 

###############################################################################
pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=350
x1=1050

xtks=15
y0=1e-3
y1=2e0

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
# Set y ticks
pl.yscale('log')
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=7)
pl.ylabel(r"$Normalized$ $Response$",fontsize=7)
pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
pl.title("Normalized Response",fontsize=9)
pl.plot(Vega20130921UT[:,0]/10.,Vega20130921UT[:,1],label='20130921UT',linewidth=1)
#pl.plot(Vega20140902UT[:,0]/10.,Vega20140902UT[:,1],label='20140902UT',linewidth=0.5)
#pl.plot(Vega20140916UT[:,0]/10.,Vega20140916UT[:,1],label='20140916UT',linewidth=0.5)
#pl.plot(Vega20140916UT[:,0]/10.,VegaAvg_20140902and16,label='20140902-16UT Avg',linewidth=1,color='k')

pl.plot(Mean200linespermm1260mm[:,0],Mean200linespermm1260mm[:,1],label='200lpm',linewidth=1.0,color='k')
pl.plot(MeanNIR[:,0],MeanNIR[:,1]*0.295*0.740,label='zNIR',linewidth=1.0,color='C3')
pl.plot(MeanRED[:,0],MeanRED[:,1]*0.580,label='zRED',linewidth=1.0,color='r')
pl.plot(MeanGRN[:,0],MeanGRN[:,1]*0.943,label='zGRN',linewidth=1.0,color='g')
pl.plot(MeanBLU[:,0],MeanBLU[:,1]*0.980,label='zBLU',linewidth=1.0,color='b')
pl.plot(MeanNUV[:,0],MeanNUV[:,1]*0.028,label='zBLU',linewidth=1.0,color='m')

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Flux Calibration/"

pylab.savefig(path+'FluxCalibrationYears.png',dpi=300)

np.savetxt(path+'FluxCalibrationYears.txt',Mean200linespermm1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

###############################################################################
#Label,Type,Start,End,Center,Avg.,SEM,WAvg.,WEM
#The Start and End wavelengths are the limits of consideration for the 
#computed values, not the actual FWHM.

List=[['380NUV','Filter',379.,381.,0.,0.,0.,0.,0.],
      ['450BLU','Filter',410.,510.,0.,0.,0.,0.,0.],
      ['486HIB','Line  ',485.,486.,0.,0.,0.,0.,0.],
      ['501OIII','Line  ',496.,502.,0.,0.,0.,0.,0.],
      ['550GRN','Filter',480.,570.,0.,0.,0.,0.,0.],
      ['650RED','Filter',610.,685.,0.,0.,0.,0.,0.],
      ['656HIA','Line  ',655.,657.,0.,0.,0.,0.,0.],
      ['672SII','Line  ',671.,673.,0.,0.,0.,0.,0.],
      ['685NIR','Filter',685.,1050.,0.,0.,0.,0.,0.],
      ['685-742','Filter',685.,742.,0.,0.,0.,0.,0.],
      ['714ArIII','Line  ',713.,715.,0.,0.,0.,0.,0.],
      ['742NIR','Filter',742.,1050.,0.,0.,0.,0.,0.],
      ['742-807','Filter',742.,807.,0.,0.,0.,0.,0.],
      ['775ArIII','Line  ',774.,775.,0.,0.,0.,0.,0.],
      ['807NIR','Filter',807.,1050.,0.,0.,0.,0.,0.],
      ['889CH4','Filter',888.,890.,0.,0.,0.,0.,0.],
      ['907SIII','Line  ',906.,908.,0.,0.,0.,0.,0.],
      ['953SIII','Line  ',952.,954.,0.,0.,0.,0.,0.]]

#print List[0][5]      
Outfile=path+'test.txt'
Append=False
for i in range(0,len(List)):       
    StartIndex=np.where(Mean200linespermm1260mm[:,0] == np.float(List[i][2]))
    EndIndex=np.where(Mean200linespermm1260mm[:,0] == np.float(List[i][3]))
    List[i][4]=np.mean([List[i][2],List[i][3]])
    #print StartIndex[0],EndIndex[0]
    #print Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1]
    List[i][5]=np.mean(Mean200linespermm1260mm[StartIndex[0][0]:EndIndex[0][0],1])
    #List[i][6]=scipy.stats.sem(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],ddof=0)
    #frac_sem=(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],3])/(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1])
    #test=1./frac_sem**2
    #List[i][7]=np.average(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],weights=test)
    List[i][8]=np.sqrt(np.sum(np.square(Mean200linespermm1260mm[StartIndex[0][0]:EndIndex[0][0],3])))/Mean200linespermm1260mm[StartIndex[0][0]:EndIndex[0][0],3].size
    tempstring=','.join(['%.3f' % num for num in List[i][2:9]])
    tempstring=List[i][0]+","+List[i][1]+","+tempstring+",\n"
    if Append:
        with open(Outfile, "a") as text_file:
            text_file.write(tempstring)
            text_file.close() 
    else:
        text_file = open(Outfile, "w")
        text_file.write("Label,Type,Start,End,Center,Avg.,SEM,WAvg.,WEM\n")
        text_file.write(tempstring)
        text_file.close()
        Append=True
print List
#print test
