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
sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
sys.path.append(drive+'\\Astronomy\Python Play\Utils')
#sys.path.append(drive+'\\Astronomy\Python Play\Galaxies')

import matplotlib.pyplot as pl
import pylab
import numpy as np
import SysRespLIB as SRL
import ConfigFiles as CF
import PlotUtils as PU
import GeneralSpecUtils as GSU
#RETRIEVE ST2000 RESPONSES#####################################################
path=drive+"/Astronomy/Projects/Techniques/Flux Calibration/"
# Read and reshape spectral data files    

ManCamData=SRL.manufacturer_camera_data(path+"CameraResponse-ST2000Data.txt")
ManCamData.load_records()
ManCamData.uniform_wave_grid()

C8_StarBright_Data=SRL.manufacturer_Celestrom_data(path+"CameraResponse - StarBright.txt")
C8_StarBright_Data.load_records()
C8_StarBright_Data.uniform_wave_grid()

C8_StarBrightXLT_Data=SRL.manufacturer_Celestrom_data(path+"CameraResponse - StarBrightXLT.txt")
C8_StarBrightXLT_Data.load_records()
C8_StarBrightXLT_Data.uniform_wave_grid()

Atmosphere_Data=SRL.atmosphere_data(path+"CameraResponse - Atmosphere.txt")
Atmosphere_Data.load_records()
Atmosphere_Data.uniform_wave_grid()

CLRPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
CLRPlotParams.loadplotparams(drive,"550CLR","TBD")
CLR550_ObsList=CF.measurement_list(CLRPlotParams.DataFile)
CLR550_ObsList.load_records(MeasTgt="1260mm200lpm")
Mean200linespermm1260mm=GSU.SpectrumAggregation("f:",CLR550_ObsList)
Mean200linespermm1260mm.ComputeAverageandStats()
tmp=SRL.Compute_EWs(path,"1260mm200lpm-550CLR-EW",Mean200linespermm1260mm.MeanSpec,1.0189)
np.savetxt(path+'SystemResponseCLR-1260mm200lpm.txt',Mean200linespermm1260mm.MeanSpec,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

CLR550_ObsList.load_records(MeasTgt="1260mm100lpm")
Mean100linespermm1260mm=GSU.SpectrumAggregation("f:",CLR550_ObsList)
Mean100linespermm1260mm.ComputeAverageandStats()

CLR550_ObsList.load_records(MeasTgt="135mm100lpm")
Mean100linespermm135mm=GSU.SpectrumAggregation("f:",CLR550_ObsList)
Mean100linespermm135mm.ComputeAverageandStats()
tmp=SRL.Compute_EWs(path,"135mm100lpm-550CLR-EW",Mean100linespermm135mm.MeanSpec,1.0224)
np.savetxt(path+'SystemResponseCLR-135mm100lpm.txt',Mean100linespermm135mm.MeanSpec,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

CLR550_ObsList.load_records("135mm200lpm")
Mean200linespermm135mm=GSU.SpectrumAggregation("f:",CLR550_ObsList)
Mean200linespermm135mm.ComputeAverageandStats()

#RETRIEVE FILTER RESPONSES#####################################################
NIRPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
NIRPlotParams.loadplotparams(drive,"685NIR","TBD")
NIR685_ObsList=CF.measurement_list(NIRPlotParams.DataFile)
NIR685_ObsList.load_records()
MeanNIRtmp=GSU.SpectrumAggregation("f:",NIR685_ObsList)
MeanNIRtmp.ComputeAverageandStats()
MeanNIR135mm=MeanNIRtmp.MeanSpec[430:2325,:]

REDPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
REDPlotParams.loadplotparams(drive,"650RED","TBD")
RED650_ObsList=CF.measurement_list(REDPlotParams.DataFile)
RED650_ObsList.load_records()
MeanREDtmp=GSU.SpectrumAggregation("f:",RED650_ObsList)
MeanREDtmp.ComputeAverageandStats()
MeanRED135mm=MeanREDtmp.MeanSpec[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-650RED-EW",MeanRED135mm,0.69/1.02424)
np.savetxt(path+'SystemResponseRED-135mm100lpm.txt',MeanRED135mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

GRNPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
GRNPlotParams.loadplotparams(drive,"550GRN","TBD")
GRN550_ObsList=CF.measurement_list(GRNPlotParams.DataFile)
GRN550_ObsList.load_records()
MeanGRNtmp=GSU.SpectrumAggregation("f:",GRN550_ObsList)
MeanGRNtmp.ComputeAverageandStats()
MeanGRN135mm=MeanGRNtmp.MeanSpec[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-550GRN-EW",MeanGRN135mm,0.98/1.0224)
np.savetxt(path+'SystemResponseGRN-135mm100lpm.txt',MeanGRN135mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

BLUPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
BLUPlotParams.loadplotparams(drive,"450BLU","TBD")
BLU450_ObsList=CF.measurement_list(BLUPlotParams.DataFile)
BLU450_ObsList.load_records()
MeanBLUtmp=GSU.SpectrumAggregation("f:",BLU450_ObsList)
MeanBLUtmp.ComputeAverageandStats()
MeanBLU135mm=MeanBLUtmp.MeanSpec[430:2325,:]
tmp=SRL.Compute_EWs(path,"135mm100lpm-550BLU-EW",MeanBLU135mm,0.875/1.0224)
np.savetxt(path+'SystemResponseBLU-135mm100lpm.txt',MeanBLU135mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")

NUVPlotParams=PU.PlotSetup("FluxCalPlotConfig.txt")
NUVPlotParams.loadplotparams(drive,"380NUV","TBD")
NUV380_ObsList=CF.measurement_list(NUVPlotParams.DataFile)
NUV380_ObsList.load_records()
MeanNUVtmp=GSU.SpectrumAggregation("f:",NUV380_ObsList)
MeanNUVtmp.ComputeAverageandStats()
MeanNUV135mm=MeanNUVtmp.MeanSpec[430:2325,:]

#MAKE ST2000 RESPONSE PLOT#####################################################
CLRPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(Mean200linespermm1260mm.MeanSpec,1.0189,'C0','1260mm200lpm')
tmp=PU.Draw_with_Conf_Level(Mean100linespermm1260mm.MeanSpec,1.0,'b','1260mm100lpm')
tmp=PU.Draw_with_Conf_Level(Mean100linespermm135mm.MeanSpec,1.0244,'C3','135mm100lpm')
tmp=PU.Draw_with_Conf_Level(Mean200linespermm135mm.MeanSpec,1.0122,'r','135mm200lpm')

#Plot Mfg. QE -> response to photon counts
#pl.scatter(ManCamData.Wavelength,ManCamData.ST2000_Norm,
#           label='ST2000 SBIG QE',s=20,color='b')
#Plot Mfg. response to energy (~QE/wavelength)
EnergyResponse=ManCamData.griddata/ManCamData.wavegrid*450.

#pl.scatter(ManCamData.Wavelength,np.array(ManCamData.ST2000_Norm)/np.array(ManCamData.Wavelength)*450.,
#           label='ST2000 SBIG Ergs',s=20,color='g')

#pl.plot(ManCamData.wavegrid,EnergyResponse,
#           label='ST2000 SBIG Ergs',linewidth=1.0,color='k')

#pl.plot(Atmosphere_Data.wavegrid,Atmosphere_Data.griddata,
#           label='ST2000 SBIG Ergs',linewidth=1.0,color='b')


#pl.plot(C8_StarBrightXLT_Data.wavegrid,C8_StarBright_Data.griddata*EnergyResponse,
#           label='StarBright Response',linewidth=1.0,color='k')
#pl.plot(C8_StarBrightXLT_Data.wavegrid,C8_StarBrightXLT_Data.griddata*EnergyResponse,
#           label='StarBrightXLT Response',linewidth=1.0,color='k')

pl.plot(C8_StarBright_Data.wavegrid,
        (C8_StarBright_Data.griddata*EnergyResponse*Atmosphere_Data.griddata)/0.52,
        label='StarBright Energy Response',linewidth=1.0,color='g')
pl.plot(C8_StarBright_Data.wavegrid,
        (C8_StarBright_Data.griddata*ManCamData.griddata*Atmosphere_Data.griddata)/0.62,
        label='StarBright QE Response',linewidth=1.0,color='k')
#pl.plot(C8_StarBrightXLT_Data.wavegrid,
#        (C8_StarBrightXLT_Data.griddata*EnergyResponse*Atmosphere_Data.griddata)/0.6,
#        label='XLT Model Response',linewidth=1.0,color='r')

pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'PanchromaticSystemResponse.png',dpi=300)

#MAKE 135MM-ST2000 PLOT WITH FILTERS###########################################
CLRPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(Mean100linespermm135mm.MeanSpec,1.0244,'0.5','135mm100lpm')
tmp=PU.Draw_with_Conf_Level(MeanNIR135mm,0.40,'C3','NIR')
tmp=PU.Draw_with_Conf_Level(MeanRED135mm,0.690,'r','RED')
tmp=PU.Draw_with_Conf_Level(MeanGRN135mm,0.980,'g','GRN')
tmp=PU.Draw_with_Conf_Level(MeanBLU135mm,0.875,'b','BLU')
tmp=PU.Draw_with_Conf_Level(MeanNUV135mm,0.016,'m','NUV')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'FilterRelativeSystemResponse135mm.png',dpi=300)

#MAKE TRANSMISSION PLOT########################################################
CLRPlotParams.Setup_Plot()
#Should be able to make a transmission computing function with two inputs
print "MeanNIR.shape:",MeanNIR135mm.shape

TransNIR=GSU.SpectrumMath(MeanNIR135mm,Mean100linespermm135mm.MeanSpec,"Divide")
tmp=PU.Draw_with_Conf_Level(TransNIR,0.40/1.0244,'C3','NIR')
TransRED=GSU.SpectrumMath(MeanRED135mm,Mean100linespermm135mm.MeanSpec,"Divide")
tmp=PU.Draw_with_Conf_Level(TransRED,0.69/1.0244,'r','RED')
TransGRN=GSU.SpectrumMath(MeanGRN135mm,Mean100linespermm135mm.MeanSpec,"Divide")
tmp=PU.Draw_with_Conf_Level(TransGRN,0.98/1.0244,'g','GRN')
TransBLU=GSU.SpectrumMath(MeanBLU135mm,Mean100linespermm135mm.MeanSpec,"Divide")
tmp=PU.Draw_with_Conf_Level(TransBLU,0.875/1.0244,'b','BLU')
TransNUV=GSU.SpectrumMath(MeanNUV135mm,Mean100linespermm135mm.MeanSpec,"Divide")
tmp=PU.Draw_with_Conf_Level(TransNUV,0.016/01.244,'m','NUV')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'FilterTransmission.png',dpi=300)

MeanNIR1260mm=GSU.SpectrumMath(TransNIR,Mean100linespermm1260mm.MeanSpec,"Multiply")
MeanRED1260mm=GSU.SpectrumMath(TransRED,Mean100linespermm1260mm.MeanSpec,"Multiply")
MeanGRN1260mm=GSU.SpectrumMath(TransGRN,Mean100linespermm1260mm.MeanSpec,"Multiply")
MeanBLU1260mm=GSU.SpectrumMath(TransBLU,Mean100linespermm1260mm.MeanSpec,"Multiply")
MeanNUV1260mm=GSU.SpectrumMath(TransNUV,Mean100linespermm1260mm.MeanSpec,"Multiply")

CLRPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(Mean200linespermm1260mm.MeanSpec,1.0,'0.5','1260mm100lpm')
tmp=PU.Draw_with_Conf_Level(MeanNIR1260mm,0.32/1.0244,'C3','NIR')
tmp=PU.Draw_with_Conf_Level(MeanRED1260mm,0.62/1.0244,'r','RED')
tmp=PU.Draw_with_Conf_Level(MeanGRN1260mm,0.98*0.94/1.0244,'g','GRN')
tmp=PU.Draw_with_Conf_Level(MeanBLU1260mm,0.875/1.0244,'b','BLU')
tmp=PU.Draw_with_Conf_Level(MeanNUV1260mm,0.011/01.244,'m','NUV')
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'FilterRelativeSystemResponse1260mm.png',dpi=300)

REDPlotParams.Setup_Plot()
print "SHAPE: ",TransRED.shape
tmp=PU.Draw_with_Conf_Level(TransRED,0.69/1.0244,'k','RED Transmission')
tmp=PU.Draw_with_Conf_Level(MeanRED135mm,0.69,'r','RED 135mm Response')
tmp=PU.Draw_with_Conf_Level(MeanRED1260mm,0.62/1.0244,'C3','RED 1260mm Response')
tmp=SRL.Compute_EWs(path,"1260mm200lpm-650RED-EW",MeanRED1260mm,0.62/1.0244)
np.savetxt(path+'SystemResponseRED-1260mm200lpm.txt',MeanRED1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")
np.savetxt(path+'TransmissionRED.txt',TransRED,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'REDTransmission.png',dpi=300)

GRNPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(TransGRN,0.98/1.0244,'k','GRN Transmission')
tmp=PU.Draw_with_Conf_Level(MeanGRN135mm,0.98,'g','GRN 135mmResponse')
tmp=PU.Draw_with_Conf_Level(MeanGRN1260mm,0.98*0.94/1.0244,'C2','GRN 1260mm Response')
tmp=SRL.Compute_EWs(path,"1260mm200lpm-450GRN-EW",MeanGRN1260mm,0.98*0.94/1.0244)
np.savetxt(path+'SystemResponseGRN-1260mm200lpm.txt',MeanGRN1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'GRNTransmission.png',dpi=300)

BLUPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(TransBLU,0.875/1.0244,'k','BLU Transmission')
tmp=PU.Draw_with_Conf_Level(MeanBLU135mm,0.875,'C0','BLU135mm Response')
tmp=PU.Draw_with_Conf_Level(MeanBLU1260mm,0.875/1.0244,'C9','BLU 1260mm Response')
tmp=SRL.Compute_EWs(path,"1260mm200lpm-550BLU-EW",MeanBLU1260mm,0.875/1.0244)
np.savetxt(path+'SystemResponseBLU-1260mm200lpm.txt',MeanBLU1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'BLUTransmission.png',dpi=300)

NUVPlotParams.Setup_Plot()
tmp=PU.Draw_with_Conf_Level(MeanNUV135mm,1.0,'k','NUV Response')
#tmp=PU.Draw_with_Conf_Level(MeanBLU135mm,0.875,'C0','BLU135mm Response')
#tmp=PU.Draw_with_Conf_Level(MeanBLU1260mm,0.875/1.0244,'C9','BLU 1260mm Response')
#tmp=SRL.Compute_EWs(path,"1260mm200lpm-550BLU-EW",MeanBLU1260mm,0.875/1.0244)
np.savetxt(path+'SystemResponseNUV-1260mm200lpm.txt',MeanBLU1260mm,delimiter=" ",
           fmt="%10.3F %10.7F %10.7F %10.7F")
pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)
pylab.savefig(path+'NUVTransmission.png',dpi=300)