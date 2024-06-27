#*****************************************************************************#
#---Model function------------------------------------------------------------#
#*****************************************************************************#


#-----------------------------------------------------------------------------#
#---For IV; without enterohepatic circulation---------------------------------#

IV.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                        #Abbreviation (C for concentration and then 2 first letters of tissue name)
    #CIn = AIn / VIn;                                                                #ZEN concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                                #Zen concentration in liver (umol/L)
    CKi= AKi / VKi;                                                                 #ZEN concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                                #ZEN concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                                #ZEN concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                                #ZEN concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                                #ZEN concentration in blood (umol/L)
    
    #Tissue venous concentrations (umol/L)                                          #Only for perfused organs; Equals tissue concentration/tissue:blood partition coefficient
    CvLi = CLi / PLib;                                                              #ZEN concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                              #ZEN concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                              #ZEN concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                              #ZEN concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                              #ZEN concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    # 
    # #Differential equations gastro-intestinal tract - not perfused
    # d/dt(AIn) = - (ka*AIn) - (kTsi*AIn);                                            #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    # d/dt(ACo) = (kTsi*AIn) - (kTco*ACo);                                            #(umol/h) rate of change in amount of ZEN in colon
    # d/dt(AFe) = kTco*ACo;                                                           #(umol/h) rate of change in amount of ZEN eliminated in the feces
    # 
    #Differential equations perfused tissues
    d/dt(ALi) = QLi*(CBl - CvLi) + #(ka*AIn) 
      - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi_clearance) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free))                                                           #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                                 #(umol/h) rate of change in amount in the kidney
    d/dt(AKi_clearance) = CLRe*CKi_free;                                            #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                   #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                   #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                   #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                         #(umol/h) rate of change in amount in the blood
    
  }) 
  return(ModelZEN)
}


#---End-----------------------------------------------------------------------#
#*****************************************************************************#