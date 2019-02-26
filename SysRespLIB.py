# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 09:21:01 2017

    This library is intended to read, manipulate and write spectral data,
    specifically for the purpose of analyzing observing system response. 
    
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
    def load_records(self):
        
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
    def load_records(self):
        
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
    def load_records(self):
        
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
