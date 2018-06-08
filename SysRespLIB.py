# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 09:21:01 2017

    This library is intended to read, manipulate and write spectral data,
    specifically for the purpose of analyzing observing system response. 
    
    0CLASS SysResp_plot_params
    0  INIT
    0  Setup_Plot
    0FUNCTION Draw_with_Conf_Level(Data,scl,clr,lbl)
    1CLASS measurement_list
    1  INIT
    1  Load_all_data
    1  Load_select_data
    1FUNCTION GetObsFileNames(Path,IndexFile)
    2CLASS SpectrumAggregation(drive,ObsList)
    2  INIT
    2  ComputeAverageandStats
    2FUNCTION Compute_Transmission(Response_with_Filter,Response_without_Filter)
    2FUNCTION Compute_EWs(path,outfile,Spectrum_with_Stats)
 
@author: SM Hill
Update 2/11/2018
"""
import sys
drive='f:'
sys.path.append(drive+'\\Astronomy\Python Play\Util')

import ConfigFiles as CF


def Draw_with_Conf_Level(Data,scl,clr,lbl):                
#Plot Layout Configuration
    import pylab as pl
    pl.plot(Data[:,0],Data[:,1]*scl,label=lbl,linewidth=1.0,color=clr)
    pl.plot(Data[:,0],(Data[:,1]+1.96*Data[:,3])*scl,linewidth=0.2,color=clr)
    pl.plot(Data[:,0],(Data[:,1]-1.96*Data[:,3])*scl,linewidth=0.2,color=clr)
    #ax.fill_between((Data[:,0]),(Data[:,1]+1.96*Data[:,3])*scl,(Data[:,1]-1.96*Data[:,3])*scl)
    return 0        
  
class measurement_list(CF.readtextfilelines):
    pass
    def load_all_data(self):
        
        self.MeasTarget=['']   #Keyword for star identification
        self.DataType=['']           #Target, e.g., component of a multiple star
        self.DataTarget=['']           #Target, e.g., component of a multiple star
        self.DateUT=['']           #UT Date of observation: YYYYMMDDUT
        self.Optics=['']       #Instrument code, to be used for aperture
        self.Camera=['']       #Instrument code, to be used for aperture
        self.Grating=['']    #Grating 100lpm or 200lpm or None
        self.FileList=['']         #List of observation image files (FITS)
        self.NObs=0               #Number of observatinos
        FirstTime=True

        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            if FirstTime:
                self.MeasTarget[0]=str(fields[0])
                self.DataType[0]=str(fields[1])
                self.DataTarget[0]=str(fields[2])
                self.DateUT[0]=str(fields[3])
                self.Optics[0]=str(fields[4])
                self.Camera[0]=str(fields[5])
                self.Grating[0]=str(fields[6])
                self.FileList[0]=str(fields[7])
                FirstTime=False
                self.NObs=1
            else:
                self.MeasTarget.extend([str(fields[0])])
                self.DataType.extend([str(fields[1])])
                self.DataTarget.extend([str(fields[2])])
                self.DateUT.extend([str(fields[3])])
                self.Optics.extend([str(fields[4])])
                self.Camera.extend([str(fields[5])])
                self.Grating.extend([str(fields[6])])
                self.FileList.extend([str(fields[7])])
                self.NObs=self.NObs+1

    def load_select_data(self,MeasTgt):
        
        self.MeasTarget=['']   #Keyword for star identification
        self.DataType=['']           #Target, e.g., component of a multiple star
        self.DataTarget=['']           #Target, e.g., component of a multiple star
        self.DateUT=['']           #UT Date of observation: YYYYMMDDUT
        self.Optics=['']       #Instrument code, to be used for aperture
        self.Camera=['']       #Instrument code, to be used for aperture
        self.Grating=['']    #Grating 100lpm or 200lpm or None
        self.FileList=['']         #List of observation image files (FITS)
        self.NObs=0                #Number of observatinos
        FirstTime=True

        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            if fields[0]==MeasTgt:
                if FirstTime:
                    self.MeasTarget[0]=str(fields[0])
                    self.DataType[0]=str(fields[1])
                    self.DataTarget[0]=str(fields[2])
                    self.DateUT[0]=str(fields[3])
                    self.Optics[0]=str(fields[4])
                    self.Camera[0]=str(fields[5])
                    self.Grating[0]=str(fields[6])
                    self.FileList[0]=str(fields[7])
                    FirstTime=False
                    self.NObs=1
                else:
                    self.MeasTarget.extend([str(fields[0])])
                    self.DataType.extend([str(fields[1])])
                    self.DataTarget.extend([str(fields[2])])
                    self.DateUT.extend([str(fields[3])])
                    self.Optics.extend([str(fields[4])])
                    self.Camera.extend([str(fields[5])])
                    self.Grating.extend([str(fields[6])])
                    self.FileList.extend([str(fields[7])])
                    self.NObs=self.NObs+1

def GetObsFileNames(Path,IndexFile):
    """
    A base class for reading a list of data files (observations). This should
    also be a single base class used for the photometry project.
    """
    CfgFile=open(IndexFile,'r')
    CfgLines=CfgFile.readlines()
    CfgFile.close()
    nrecords=len(CfgLines)
    #print CfgLines
    FNArray=['']
    FirstTime=True
    for recordindex in range(0,nrecords):
        fields=CfgLines[recordindex].split(',')
        if FirstTime:
            FNArray[0]=str(fields[0])
            FirstTime=False
        else:
            FNArray.extend([str(fields[0])])
            
    return FNArray

class SpectrumAggregation:
    def __init__(self,drive,ObsList):
#def Average_Spectrum(drive,ObsList):
        import numpy as np
        import scipy
        self.path=drive+"/Astronomy/Python Play/TechniquesLibrary/"
        #print path
        for i in range(0,len(ObsList.FileList)):
            print "******** i=",i,ObsList.FileList[i]
            self.FNList=GetObsFileNames(self.path,ObsList.FileList[i])
            print "********len(FNList)",len(self.FNList),self.FNList
            if ObsList.DataType[i]=="Reference":
                path=drive+"/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/"
            else:    
                path=drive+"/Astronomy/Projects/"+ObsList.DataType[i]+"/"+ObsList.DataTarget[i]+"/Spectral Data/1D Spectra/"
            
            ###Need loop over data files here!!! "j"
            for j in range(0,len(self.FNList)):
                print "****** j=",j
                if ObsList.DataType[i]=="Reference":
                    temp2 = scipy.loadtxt(path+self.FNList[j], dtype=float, usecols=(0,1))
                else:
                    temp1 = scipy.fromfile(file=path+self.FNList[j], dtype=float, count=-1, sep='\t')    
                    temp2 = scipy.reshape(temp1,[temp1.size/2,2])
                self.wave = temp2[:,0]
                tmpsig=temp2[:,1]
                print i
                print temp2.shape
        
                if i==0 and j==0:
                    self.signalarray=np.zeros([tmpsig.size,1])
                    print self.signalarray.shape
                    self.signalarray[:,0]=tmpsig
                else:
                    print "i>0 self.signalarray.shape:",self.signalarray.shape
                    print "temp2.shape: ",temp2.shape
                    print "temp2[0,0]: ",temp2[0,0]
                    self.signalarray=np.insert(self.signalarray,1,tmpsig,axis=1)
                    print "i>0 self.signalarray.shape:",self.signalarray.shape
        if ObsList.DataType[i]=="Reference":
            self.wave=self.wave/10.

    def ComputeAverageandStats(self):
        import scipy.stats as ST
        import numpy as np
        ZeroIndices=np.where(self.signalarray <= 0.)
        self.signalarray[ZeroIndices]=np.nan
        #pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")
        #pl.plot(wave,signalarray[:,0])        
        AvgSignal=np.nanmean(self.signalarray,axis=1)
        std=np.nanstd(self.signalarray,axis=1) 
        sem=ST.sem(self.signalarray,axis=1,ddof=0,nan_policy='omit')
            
        self.MeanSpec=np.zeros([self.wave.size,4])
        self.MeanSpec[:,0]=self.wave
        self.MeanSpec[:,1]=AvgSignal
        self.MeanSpec[:,2]=std
        self.MeanSpec[:,3]=sem
        
        return 0
    
def SpectrumMath(Spectrum1,Spectrum2,Operation):
    from copy import deepcopy
    import numpy as np
    ResultSpectrum=deepcopy(Spectrum1)
    if Operation == "Divide":
        ResultSpectrum[:,1]=Spectrum1[:,1]/Spectrum2[:,1]
    elif Operation == "Multiply":
        ResultSpectrum[:,1]=Spectrum1[:,1]*Spectrum2[:,1]
    #Combine uncertainty arrays in quadrature
    ResultSpectrum[:,2]=ResultSpectrum[:,1]*np.sqrt((Spectrum1[:,2]/Spectrum1[:,1])**2+
        (Spectrum2[:,2]/Spectrum2[:,1])**2,)
    ResultSpectrum[:,3]=ResultSpectrum[:,1]*np.sqrt((Spectrum1[:,3]/Spectrum1[:,1])**2+
        (Spectrum2[:,3]/Spectrum2[:,1])**2,)
        
    return ResultSpectrum

def Compute_EWs(path,outfile,Spectrum_with_Stats,Scale):
###############################################################################
#Label,Type,Start,End,Center,Avg. Response,SEM Response,WAvg.,WEM
#The Start and End wavelengths are the limits of consideration for the 
#computed values, not the actual FWHM.
    import numpy as np

    List=[['380NUV','Filter',379.,381.,0.,0.,0.,0.,0.],
          ['450BLU','Filter',410.,510.,0.,0.,0.,0.,0.],
          ['486HIB','Line  ',485.,486.,0.,0.,0.,0.,0.],
          ['501OIII','Line  ',496.,502.,0.,0.,0.,0.,0.],
          ['550GRN','Filter',480.,570.,0.,0.,0.,0.,0.],
          ['650RED','Filter',610.,685.,0.,0.,0.,0.,0.],
          ['650RED-A','Filter',620.,670.,0.,0.,0.,0.,0.],
          ['656HIA--','Line  ',636.,646.,0.,0.,0.,0.,0.],
          ['656HIA-','Line  ',649.,652.,0.,0.,0.,0.,0.],
          ['656HIA','Line  ',655.,657.,0.,0.,0.,0.,0.],
          ['656HIA+','Line  ',661.,664.,0.,0.,0.,0.,0.],
          ['656HIA++','Line  ',666.,676.,0.,0.,0.,0.,0.],
          ['658NII--','Line  ',638.,648.,0.,0.,0.,0.,0.],
          ['658NII-','Line  ',649.,652.,0.,0.,0.,0.,0.],
          ['658NII','Line  ',657.,659.,0.,0.,0.,0.,0.],
          ['658NII+','Line  ',663.,666.,0.,0.,0.,0.,0.],
          ['658NII++','Line  ',668.,678.,0.,0.,0.,0.,0.],
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
    Outfile=path+outfile+'.txt'
    Append=False
    for i in range(0,len(List)):       
        StartIndex=np.where(Spectrum_with_Stats[:,0] == np.float(List[i][2]))
        EndIndex=np.where(Spectrum_with_Stats[:,0] == np.float(List[i][3]))
        List[i][4]=np.nanmean([List[i][2],List[i][3]])
        #print StartIndex[0],EndIndex[0]
        #print Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1]
        List[i][5]=np.nanmean(Spectrum_with_Stats[StartIndex[0][0]:EndIndex[0][0],1])*Scale
        #List[i][6]=scipy.stats.sem(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],ddof=0)
        #frac_sem=(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],3])/(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1])
        #test=1./frac_sem**2
        #List[i][7]=np.average(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],weights=test)
        List[i][8]=Scale*(np.sqrt(np.sum(np.square(Spectrum_with_Stats[StartIndex[0][0]:EndIndex[0][0],3])))/Spectrum_with_Stats[StartIndex[0][0]:EndIndex[0][0],3].size)
        tempstring=','.join(['%.3f' % num for num in List[i][2:9]])
        tempstring=List[i][0]+","+List[i][1]+","+tempstring+",\n"
        if Append:
            with open(Outfile, "a") as text_file:
                text_file.write(tempstring)
                text_file.close() 
        else:
            text_file = open(Outfile, "w")
            text_file.write("Label ,Type  ,Start  ,End    ,Center ,Avg. ,SEM  ,WAvg. ,WEM  \n")
            text_file.write(tempstring)
            text_file.close()
            Append=True
    #print List
    #print test
    
class manufacturer_camera_data(CF.readtextfilelines):
    pass
    def load_all_data(self):
        
        self.Wavelength=[0.]   #Keyword for star identification
        self.ST2000_QE=[0.]           #Target, e.g., component of a multiple star
        self.ST2000_Norm=[0.]           #UT Date of observation: YYYYMMDDUT
        self.NObs=0                #Number of observatinos
        FirstTime=True

        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            print fields
            if FirstTime:
                self.Wavelength[0]=float(fields[0])
                self.ST2000_QE[0]=float(fields[1])
                self.ST2000_Norm[0]=float(fields[2])
                FirstTime=False
                self.NObs=1
            else:
                self.Wavelength.extend([float(fields[0])])
                self.ST2000_QE.extend([float(fields[1])])
                self.ST2000_Norm.extend([float(fields[2])])
                self.NObs=self.NObs+1
                
    def uniform_wave_grid(self):
        import numpy as np
        from scipy import interpolate
        self.wavegrid=np.arange(115,1062.5,0.5,dtype=float)
        Interp=interpolate.interp1d(self.Wavelength,self.ST2000_Norm,kind='linear', 
                                    copy=True,bounds_error=False, 
                                    fill_value=np.NaN,axis=0)  
        self.griddata=Interp(self.wavegrid)
                
class manufacturer_Celestrom_data(CF.readtextfilelines):
    pass
    def load_all_data(self):
        
        self.Wavelength=[0.]   #Keyword for star identification
        self.Transmission=[0.]           #Target, e.g., component of a multiple star
        self.NObs=0                #Number of observatinos
        FirstTime=True

        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            if FirstTime:
                self.Wavelength[0]=float(fields[0])
                self.Transmission[0]=float(fields[1])
                FirstTime=False
                self.NObs=1
            else:
                self.Wavelength.extend([float(fields[0])])
                self.Transmission.extend([float(fields[1])])
                self.NObs=self.NObs+1

    def uniform_wave_grid(self):
        import numpy as np
        from scipy import interpolate
        self.wavegrid=np.arange(115,1062.5,0.5,dtype=float)
        Interp=interpolate.interp1d(self.Wavelength,self.Transmission,kind='linear', 
                                    copy=True,bounds_error=False, 
                                    fill_value=np.NaN,axis=0)  
        self.griddata=Interp(self.wavegrid)

class atmosphere_data(CF.readtextfilelines):
    pass
    def load_all_data(self):
        
        self.Wavelength=[0.]   #Keyword for star identification
        self.Transmission=[0.]           #Target, e.g., component of a multiple star
        self.NormTransmission=[0.]           #Target, e.g., component of a multiple star
        self.NObs=0                #Number of observatinos
        FirstTime=True

        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            if FirstTime:
                self.Wavelength[0]=float(fields[0])
                self.Transmission[0]=float(fields[1])
                self.NormTransmission[0]=float(fields[2])
                FirstTime=False
                self.NObs=1
            else:
                self.Wavelength.extend([float(fields[0])])
                self.Transmission.extend([float(fields[1])])
                self.NormTransmission.extend([float(fields[2])])
                self.NObs=self.NObs+1

    def uniform_wave_grid(self):
        import numpy as np
        from scipy import interpolate
        self.wavegrid=np.arange(115,1062.5,0.5,dtype=float)
        Interp=interpolate.interp1d(self.Wavelength,self.NormTransmission,kind='linear', 
                                    copy=True,bounds_error=False, 
                                    fill_value=np.NaN,axis=0)  
        self.griddata=Interp(self.wavegrid)
