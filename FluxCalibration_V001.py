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

import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate
from copy import deepcopy

path=drive+"/Astronomy/Projects/Stars/Vega/Spectral Data/1D Spectra/"
# Read and reshape spectral data files    
Vega20130921UT = scipy.fromfile(file=path+"VegaResponse20130921UT.txt", dtype=float, count=-1, sep='\t')    
Vega20130921UT=scipy.reshape(Vega20130921UT,[Vega20130921UT.size/2,2])

#Vega20140902UT = scipy.fromfile(file=path+"VegaResponse20140902UT.txt", dtype=float, count=-1, sep='\t')    
#Vega20140902UT=scipy.reshape(Vega20140902UT,[Vega20140902UT.size/2,2])

#Vega20140916UT = scipy.fromfile(file=path+"VegaResponse20140916UT.txt", dtype=float, count=-1, sep='\t')    
#Vega20140916UT=scipy.reshape(Vega20140916UT,[Vega20140916UT.size/2,2])

Vega20150913UT = scipy.fromfile(file=path+"Vega_20150913025152_Response.txt", dtype=float, count=-1, sep='\t')    
Vega20150913UT=scipy.reshape(Vega20150913UT,[Vega20150913UT.size/2,2])

ZeroIndices=np.where(Vega20150913UT[:,1] == 0.)
Vega20150913UT[ZeroIndices,1]=np.nan


Vega20150925UT = scipy.fromfile(file=path+"Vega_20150925032310_Response.txt", dtype=float, count=-1, sep='\t')    
Vega20150925UT=scipy.reshape(Vega20150925UT,[Vega20150925UT.size/2,2])

Vega20151014UT = scipy.fromfile(file=path+"Vega_20151014025716_Response.txt", dtype=float, count=-1, sep='\t')    
Vega20151014UT=scipy.reshape(Vega20151014UT,[Vega20151014UT.size/2,2])

Vega20151015UT = scipy.fromfile(file=path+"Vega20151015023120_Response.txt", dtype=float, count=-1, sep='\t')    
Vega20151015UT=scipy.reshape(Vega20151015UT,[Vega20151015UT.size/2,2])

temparray=[Vega20150913UT[:,1],Vega20150925UT[:,1],Vega20151014UT[:,1],Vega20151015UT[:,1]]    

#where is 20151015
z=np.nanmean(temparray,axis=0)
std=np.nanstd(temparray,axis=0) 
#sem=scipy.stats.sem(temparray,axis=0,ddof=0,nan_policy='omit')
Mean200linespermm1260mm=np.zeros([Vega20150913UT.size/2,4])
Mean200linespermm1260mm[:,0]=Vega20150913UT[:,0]
Mean200linespermm1260mm[:,1]=z
Mean200linespermm1260mm[:,2]=std
Mean200linespermm1260mm[:,3]=std#sem

path=drive+"/Astronomy/Projects/Stars/Aldebaran/Spectral Data/1D Spectra/"
NIRAldebaran20130117UT = scipy.fromfile(file=path+"Aldebaran-685NIR-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
NIRAldebaran20130117UT=scipy.reshape(NIRAldebaran20130117UT,[NIRAldebaran20130117UT.size/2,2])
path=drive+"/Astronomy/Projects/Stars/Rigel/Spectral Data/1D Spectra/"
NIRRigel20130117UT = scipy.fromfile(file=path+"Rigel-685NIR-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
NIRRigel20130117UT=scipy.reshape(NIRRigel20130117UT,[NIRRigel20130117UT.size/2,2])
temparray=[NIRAldebaran20130117UT[:,1],NIRRigel20130117UT[:,1]]    
zNIR=np.nanmean(temparray,axis=0)
MeanNIR=np.zeros([NIRAldebaran20130117UT.size/2,4])
MeanNIR[:,0]=NIRAldebaran20130117UT[:,0]
MeanNIR[:,1]=zNIR

path=drive+"/Astronomy/Projects/Stars/Aldebaran/Spectral Data/1D Spectra/"
REDAldebaran20130117UT = scipy.fromfile(file=path+"Aldebaran-650RED-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
REDAldebaran20130117UT=scipy.reshape(REDAldebaran20130117UT,[REDAldebaran20130117UT.size/2,2])
path=drive+"/Astronomy/Projects/Stars/Rigel/Spectral Data/1D Spectra/"
REDRigel20130117UT = scipy.fromfile(file=path+"Rigel-650RED-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
REDRigel20130117UT=scipy.reshape(REDRigel20130117UT,[REDRigel20130117UT.size/2,2])
temparray=[REDAldebaran20130117UT[:,1],REDRigel20130117UT[:,1]]    
zRED=np.nanmean(temparray,axis=0)
MeanRED=np.zeros([REDAldebaran20130117UT.size/2,4])
MeanRED[:,0]=REDAldebaran20130117UT[:,0]
MeanRED[:,1]=zRED

path=drive+"/Astronomy/Projects/Stars/Aldebaran/Spectral Data/1D Spectra/"
GRNAldebaran20130117UT = scipy.fromfile(file=path+"Aldebaran-550GRN-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
GRNAldebaran20130117UT=scipy.reshape(GRNAldebaran20130117UT,[GRNAldebaran20130117UT.size/2,2])
path=drive+"/Astronomy/Projects/Stars/Rigel/Spectral Data/1D Spectra/"
GRNRigel20130117UT = scipy.fromfile(file=path+"Rigel-550GRN-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
GRNRigel20130117UT=scipy.reshape(GRNRigel20130117UT,[GRNRigel20130117UT.size/2,2])
temparray=[GRNAldebaran20130117UT[:,1],GRNRigel20130117UT[:,1]]    
zGRN=np.nanmean(temparray,axis=0)
MeanGRN=np.zeros([GRNAldebaran20130117UT.size/2,4])
MeanGRN[:,0]=GRNAldebaran20130117UT[:,0]
MeanGRN[:,1]=zGRN

path=drive+"/Astronomy/Projects/Stars/Aldebaran/Spectral Data/1D Spectra/"
BLUAldebaran20130117UT = scipy.fromfile(file=path+"Aldebaran-450BLU-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
BLUAldebaran20130117UT=scipy.reshape(BLUAldebaran20130117UT,[BLUAldebaran20130117UT.size/2,2])
path=drive+"/Astronomy/Projects/Stars/Rigel/Spectral Data/1D Spectra/"
BLURigel20130117UT = scipy.fromfile(file=path+"Rigel-450BLU-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
BLURigel20130117UT=scipy.reshape(BLURigel20130117UT,[BLURigel20130117UT.size/2,2])
temparray=[BLUAldebaran20130117UT[:,1],BLURigel20130117UT[:,1]]    
zBLU=np.nanmean(temparray,axis=0)
MeanBLU=np.zeros([BLUAldebaran20130117UT.size/2,4])
MeanBLU[:,0]=BLUAldebaran20130117UT[:,0]
MeanBLU[:,1]=zBLU

path=drive+"/Astronomy/Projects/Stars/Rigel/Spectral Data/1D Spectra/"
NUVRigel20130117UT = scipy.fromfile(file=path+"Rigel-380NUV-20130117UT_Response.txt", dtype=float, count=-1, sep='\t')    
NUVRigel20130117UT=scipy.reshape(NUVRigel20130117UT,[NUVRigel20130117UT.size/2,2])
temparray=[NUVRigel20130117UT[:,1]]    
zNUV=np.nanmean(temparray,axis=0)
MeanNUV=np.zeros([NUVRigel20130117UT.size/2,4])
MeanNUV[:,0]=NUVRigel20130117UT[:,0]
MeanNUV[:,1]=zNUV

#compute mean ignoring NANs   

#Vega20140916UT_Interp=interpolate.interp1d(Vega20140916UT[:,0],Vega20140916UT[:,1],kind='linear', copy=True,
#                         bounds_error=False, fill_value=0.0)  
#Vega20140916UT_on_20140902UT=Vega20140916UT_Interp(Vega20140902UT[:,0])

#VegaAvg_20140916and19=Vega20140916UT_on_20140919UT
#VegaAvg_20140902and16=(Vega20140916UT_on_20140902UT+Vega20140902UT[:,1])/2.
#Begin plotting 

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
pl.ylabel(r"$Counts-s^{-1}$-$m^{-2}$-$nm^{-1}$",fontsize=7)
pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
pl.title("Vega Response - All Data",fontsize=9)
pl.plot(Vega20130921UT[:,0]/10.,Vega20130921UT[:,1],label='20130921UT',linewidth=1)
#pl.plot(Vega20140902UT[:,0]/10.,Vega20140902UT[:,1],label='20140902UT',linewidth=0.5)
#pl.plot(Vega20140916UT[:,0]/10.,Vega20140916UT[:,1],label='20140916UT',linewidth=0.5)
#pl.plot(Vega20140916UT[:,0]/10.,VegaAvg_20140902and16,label='20140902-16UT Avg',linewidth=1,color='k')

pl.plot(Vega20150913UT[:,0],Vega20150913UT[:,1],label='20150913UT',linewidth=0.5)
pl.plot(Vega20150925UT[:,0],Vega20150925UT[:,1],label='20150925UT',linewidth=0.5)
pl.plot(Vega20151014UT[:,0],Vega20151014UT[:,1],label='20151014UT',linewidth=0.5)
pl.plot(Vega20151015UT[:,0],Vega20151015UT[:,1],label='20151015UT',linewidth=0.5)

pl.plot(Mean200linespermm1260mm[:,0],Mean200linespermm1260mm[:,1],label='200lpm',linewidth=1.0)
pl.plot(MeanNIR[:,0],MeanNIR[:,1]*0.295*0.740,label='zNIR',linewidth=1.0)
pl.plot(MeanRED[:,0],MeanRED[:,1]*0.580,label='zRED',linewidth=1.0)
pl.plot(MeanGRN[:,0],MeanGRN[:,1]*0.943,label='zGRN',linewidth=1.0)
pl.plot(MeanBLU[:,0],MeanBLU[:,1]*0.980,label='zBLU',linewidth=1.0)
pl.plot(MeanNUV[:,0],MeanNUV[:,1]*0.028,label='zBLU',linewidth=1.0)


pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

path=drive+"/Astronomy/Projects/Techniques/Flux Calibration/"


pylab.savefig(path+'FluxCalibrationYears.png',dpi=300)

np.savetxt(path+'FluxCalibrationYears.txt',Mean200linespermm1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

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
    List[i][5]=np.mean(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1])
    #List[i][6]=scipy.stats.sem(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],ddof=0)
    #frac_sem=(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],3])/(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1])
    #test=1./frac_sem**2
    #List[i][7]=np.average(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],weights=test)
    List[i][8]=np.sqrt(np.sum(np.square(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],3])))/Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],3].size
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
