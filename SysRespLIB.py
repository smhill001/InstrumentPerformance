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
    #This is the first cut at a base class. It is a truncated version
    #  of the class used by the photometry code. Thus, the photometry code
    #  should be updated to build on this base class.
    def __init__(self,MeasurementListFile):
        #The initial plan is to read ALL records in the observation list
        self.MeasTarget=['']   #Keyword for star identification
        self.DataType=['']           #Target, e.g., component of a multiple star
        self.DataTarget=['']           #Target, e.g., component of a multiple star
        self.DateUT=['']           #UT Date of observation: YYYYMMDDUT
        self.Optics=['']       #Instrument code, to be used for aperture
        self.Camera=['']       #Instrument code, to be used for aperture
        self.FileList=['']         #List of observation image files (FITS)
        self.NObs=0                #Number of observatinos
        self.FirstTime=True
        
        CfgFile=open(MeasurementListFile,'r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            if self.FirstTime:
                self.MeasTarget[0]=str(fields[0])
                self.DataType[0]=str(fields[1])
                self.DataTarget[0]=str(fields[2])
                self.DateUT[0]=str(fields[3])
                self.Optics[0]=str(fields[4])
                self.Camera[0]=str(fields[5])
                self.FileList[0]=str(fields[6])
                self.FirstTime=False
                self.NObs=1
            else:
                self.MeasTarget.extend([str(fields[0])])
                self.DataType.extend([str(fields[1])])
                self.DataTarget.extend([str(fields[2])])
                self.DateUT.extend([str(fields[3])])
                self.Optics.extend([str(fields[4])])
                self.Camera.extend([str(fields[5])])
                self.FileList.extend([str(fields[6])])
                self.NObs=self.NObs+1


def Average_Spectrum(drive,ObsList):
    import numpy as np
    import scipy
    path=drive+"/Astronomy/Python Play/TechniquesLibrary/"
    print path
    for i in range(0,len(ObsList.FileList)):
        print "******** i=",i
        FNList=GetObsFileNames(path,ObsList.FileList[i])
        print "len(FNList)",len(FNList),FNList
        path=drive+"/Astronomy/Projects/"+ObsList.DataType[i]+"/"+ObsList.DataTarget[i]+"/Spectral Data/1D Spectra/"
        ###Need loop over data files here!!! "j"
        for j in range(0,len(FNList)):
            temp1 = scipy.fromfile(file=path+FNList[j], dtype=float, count=-1, sep='\t')    
            temp2 = scipy.reshape(temp1,[temp1.size/2,2])
            wave = temp2[:,0]
            tmpsig=temp2[:,1]
            #ZeroIndices=np.where(tmpsig <= 0.)
            #tmpsig[ZeroIndices]=np.nan
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
            
    AvgSignal=np.nanmean(signalarray,axis=1)
    std=np.nanstd(signalarray,axis=1) 
        
    MeanSpec=np.zeros([wave.size,4])
    MeanSpec[:,0]=wave
    MeanSpec[:,1]=AvgSignal
    MeanSpec[:,2]=std
    #Mean200linespermm1260mm[:,3]=std#sem
    
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

def PlotHII(Target,X_data,Y_data,plotparams):                
#Plot Layout Configuration
    import pylab as pl
    import numpy as np
    
    pl.grid()
    pl.xlim(plotparams.X0,plotparams.X1)
    pl.xticks(np.arange(plotparams.X0,plotparams.X1+.000001,plotparams.DX))
    pl.ylim(plotparams.Y0,plotparams.Y1)
    pl.yticks(np.arange(plotparams.Y0,plotparams.Y1+.000001,plotparams.DY))  
    pl.tick_params(axis='both', which='major', labelsize=7)
    
    pl.title("H II Radial Distribution",fontsize=9)
    if plotparams.PlotType == "POSHist":
        pl.hist(X_data,bins=Y_data,label=Target)
        pl.ylabel(r"$H$ $II$ $Region$ $(Count)$",fontsize=7)
    elif plotparams.PlotType == "POSDens":
        p=pl.plot(X_data,Y_data,label=Target,marker='.',linewidth=1.0)
        pl.xlabel(r"$Radius$ $(arcmin)$",fontsize=7)
        pl.ylabel(r"$H$ $II$ $Region$ $(Count)$",fontsize=7)
        print "COLOR=",p[0].get_color()

    elif plotparams.PlotType == "POSLoc":
        pl.xlabel(r"$Delta$ $RA$ $(arcmin)$",fontsize=7)
        pl.ylabel(r"$Delta$ $DE$ $(arcmin)$",fontsize=7)
        pl.title("H II Location Distribution",fontsize=9)
        pl.scatter(X_data,Y_data,label=Target,marker='o',edgecolor="#1f77b4",
                   linewidth =1.0,facecolor="") 
        
    pl.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})        

    return 0        

