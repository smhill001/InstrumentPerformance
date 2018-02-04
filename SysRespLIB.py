# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 09:21:01 2017

    This library is intended to read, manipulate and write spectral data,
    specifically for the purpose of analyzing observing system response. 

@author: Astronomy
"""
import sys
drive='f:'
sys.path.append(drive+'\\Astronomy\Python Play\Galaxies')
import GalaxyLIB as GL



class SysResp_plot_params(GL.Plot_Parameters):
    """
    This class builds on the base class to add parameters specific to
    HII plots. In this case, there are no additions, just the code to
    populate the object.
    
    SMH 1/11/18
    """
    def __init__(self,drive,PlotID,PlotType):
        #View has two options: raw or flux?

        self.ID=PlotID

        CfgFile=open(drive+'/Astronomy/Python Play/TechniquesLibrary/FluxCalPlotConfig.txt','r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            #print fields[0], fields[1]
            if fields[0] == PlotID:
                if fields[1] == PlotType:
                    #print "In first if, fields[1]",fields[:]
                    self.PlotType=str(fields[1])
                    self.X0=float(fields[2])
                    self.X1=float(fields[3])
                    self.DX=float(fields[4])
                    self.Xtype=str(fields[5])
                    self.Y0=float(fields[6])
                    self.Y1=float(fields[7])
                    self.DY=float(fields[8])
                    self.Ytype=str(fields[9])
                    self.DataFile=str(fields[10])

        
class measurement_list:
    def __init__(self,MeasurementListFile):
        #The initial plan is to read ALL records in the observation list

        CfgFile=open(MeasurementListFile,'r')
        self.CfgLines=CfgFile.readlines()
        CfgFile.close()
        self.nrecords=len(self.CfgLines)
        self.MeasurmentListFile=MeasurementListFile

    def load_all_data(self):
        
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


def Average_Spectrum(drive,ObsList):
    import numpy as np
    import scipy
    import scipy.stats as ST
    path=drive+"/Astronomy/Python Play/TechniquesLibrary/"
    print path
    for i in range(0,len(ObsList.FileList)):
        print "******** i=",i
        FNList=GetObsFileNames(path,ObsList.FileList[i])
        print "********len(FNList)",len(FNList),FNList
        path=drive+"/Astronomy/Projects/"+ObsList.DataType[i]+"/"+ObsList.DataTarget[i]+"/Spectral Data/1D Spectra/"
        ###Need loop over data files here!!! "j"
        for j in range(0,len(FNList)):
            print "****** j=",j
            temp1 = scipy.fromfile(file=path+FNList[j], dtype=float, count=-1, sep='\t')    
            temp2 = scipy.reshape(temp1,[temp1.size/2,2])
            wave = temp2[:,0]
            tmpsig=temp2[:,1]
            print i
            print temp2.shape
    
            if i==0 and j==0:
                signalarray=np.zeros([tmpsig.size,1])
                print signalarray.shape
                signalarray[:,0]=tmpsig
            else:
                signalarray=np.insert(signalarray,1,tmpsig,axis=1)
                print "i>0:",signalarray.shape
            
    ZeroIndices=np.where(signalarray <= 0.)
    signalarray[ZeroIndices]=np.nan
    #pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")
    #pl.plot(wave,signalarray[:,0])        
    AvgSignal=np.nanmean(signalarray,axis=1)
    std=np.nanstd(signalarray,axis=1) 
    sem=ST.sem(signalarray,axis=1,ddof=0,nan_policy='omit')
        
    MeanSpec=np.zeros([wave.size,4])
    MeanSpec[:,0]=wave
    MeanSpec[:,1]=AvgSignal
    MeanSpec[:,2]=std
    MeanSpec[:,3]=sem
    
    return MeanSpec
    
 
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

def Draw_with_Conf_Level(Data,scl,clr,lbl):                
#Plot Layout Configuration
    import pylab as pl
    pl.plot(Data[:,0],Data[:,1]*scl,label=lbl,linewidth=1.0,color=clr)
    pl.plot(Data[:,0],(Data[:,1]+1.96*Data[:,3])*scl,linewidth=0.5,color=clr)
    pl.plot(Data[:,0],(Data[:,1]-1.96*Data[:,3])*scl,linewidth=0.5,color=clr)
    #ax.fill_between((Data[:,0]),(Data[:,1]+1.96*Data[:,3])*scl,(Data[:,1]-1.96*Data[:,3])*scl)
    return 0        


def Compute_EWs(path,Spectrum_with_Stats):
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
        StartIndex=np.where(Spectrum_with_Stats[:,0] == np.float(List[i][2]))
        EndIndex=np.where(Spectrum_with_Stats[:,0] == np.float(List[i][3]))
        List[i][4]=np.mean([List[i][2],List[i][3]])
        #print StartIndex[0],EndIndex[0]
        #print Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1]
        List[i][5]=np.mean(Spectrum_with_Stats[StartIndex[0][0]:EndIndex[0][0],1])
        #List[i][6]=scipy.stats.sem(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],ddof=0)
        #frac_sem=(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],3])/(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1])
        #test=1./frac_sem**2
        #List[i][7]=np.average(Mean200linespermm1260mm[StartIndex[0]:EndIndex[0],1],weights=test)
        List[i][8]=np.sqrt(np.sum(np.square(Spectrum_with_Stats[StartIndex[0][0]:EndIndex[0][0],3])))/Spectrum_with_Stats[StartIndex[0][0]:EndIndex[0][0],3].size
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
