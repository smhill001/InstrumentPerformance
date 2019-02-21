# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 09:21:01 2017

    This library is intended to read, manipulate and write spectral data,
    specifically for the purpose of analyzing observing system response. 
    
    0CLASS measurement_list
    0  INIT
    0  Load_all_data
    0  Load_select_data
    1FUNCTION GetObsFileNames(Path,IndexFile)
    2FUNCTION Compute_EWs(path,outfile,Spectrum_with_Stats)
    3CLASS manufacturer_camera_data
    3  INIT
    3  load_all_data
    3  uniform_wave_grid
    4CLASS manufacturer_Celestrom_data
    4  INIT
    4  load_all_data
    4  uniform_wave_grid
    5CLASS atmosphere_data
    5  INIT
    5  load_all_data
    5  uniform_wave_grid
 
@author: SM Hill
"""
import sys
drive='f:'
sys.path.append(drive+'\\Astronomy\Python Play\Util')

import ConfigFiles as CF
  
class measurement_list(CF.readtextfilelines):
    pass
    def load_all_data(self):
        print "Hi in measurement_list>load_all_data"
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

    def load_select_data(self,MeasTgt,DateUTSelect="All"):
        
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
                if DateUTSelect=="All":
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
                else:
                    if DateUTSelect==fields[3]:
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
