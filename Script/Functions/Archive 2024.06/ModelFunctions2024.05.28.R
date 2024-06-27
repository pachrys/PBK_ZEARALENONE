#*****************************************************************************#
#---Model functions-----------------------------------------------------------#
#*****************************************************************************#

#*****************************************************************************#
#---For differential dissolution in the GIT - Oral Dose---------#

Diss_EHC_PO.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                      #Abbreviation (C for concentration and then 2 first letters of tissue name)
    
    #GIT Dissolved
    CSt = ASt/VSt;
    CIn1 = AIn1/VIn1;
    CIn2 = AIn2/VIn2;
    CIn3 = AIn3/VIn3;
    CIn4 = AIn4/VIn4;
    CIn5 = AIn5/VIn5;
    CIn6 = AIn6/VIn6;
    CIn7 = AIn7/VIn7;
    CIn = CIn1 + CIn2 + CIn3 + CIn4 + CIn5 + CIn6 + CIn7;
    CCo = ACo/VCo;
    
    #CIn.glu = AIn.glu / (VIn1+VIn2+VIn3+VIn4+VIn5+VIn6+VIn7);
    CLi.glu = ALi.glu / VLi;                                                      #ZEN-glu concentration in the liver (umol/L)
    
    #CIn = AIn / VIn;                                                              #ZEN concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                              #Zen concentration in liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN concentration in blood (umol/L)
    
    #Tissue venous concentrations (umol/L)                                        #Only for perfused organs; Equals tissue concentration/tissue:blood partition coefficient
    CvLi = CLi / PLib;                                                            #ZEN concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free concentrations (umol/L)
    CKi_free = CvKi*fub;
    CLi_free = CvLi*fub;
    
    d/dt(AStundis) = - AStundis*SiSt - AStundis*kTge; #what gets dissolved
    d/dt(ASt) = AStundis*SiSt - ASt*kTge - ASt*ka;
    
    d/dt(AIn1undis) = AStundis*kTge - AIn1undis*SiIn1 - AIn1undis*kTsi; #what gets dissolved
    d/dt(AIn1) = ASt*kTge + AIn1undis*SiIn1 - AIn1*kTsi - AIn1*ka;
    
    d/dt(AIn2undis) = AIn1undis*kTsi - AIn2undis*SiIn2 - AIn2undis*kTsi; #what gets dissolved
    d/dt(AIn2) = AIn1*kTsi + AIn2undis*SiIn2 - AIn2*kTsi - AIn2*ka;
    
    d/dt(AIn3undis) = AIn2undis*kTsi - AIn3undis*SiIn3 - AIn3undis*kTsi; #what gets dissolved
    d/dt(AIn3) = AIn2*kTsi + AIn3undis*SiIn3 - AIn3*kTsi - AIn3*ka;
    
    d/dt(AIn4undis) = AIn3undis*kTsi - AIn4undis*SiIn4 - AIn4undis*kTsi; #what gets dissolved
    d/dt(AIn4) = AIn3*kTsi + AIn4undis*SiIn4 - AIn4*kTsi - AIn4*ka;
    
    d/dt(AIn5undis) = AIn4undis*kTsi - AIn5undis*SiIn5 - AIn5undis*kTsi; #what gets dissolved
    d/dt(AIn5) = AIn4*kTsi + AIn5undis*SiIn5 - AIn5*kTsi - AIn5*ka;
    
    d/dt(AIn6undis) = AIn5undis*kTsi - AIn6undis*SiIn6 - AIn6undis*kTsi; #what gets dissolved
    d/dt(AIn6) = AIn5*kTsi + AIn6undis*SiIn6 - AIn6*kTsi - AIn6*ka;
    
    d/dt(AIn7undis) = AIn6undis*kTsi - AIn7undis*SiIn7 - AIn7undis*kTsi; #what gets dissolved
    d/dt(AIn7) = AIn6*kTsi + AIn7undis*SiIn7 - AIn7*kTsi - AIn7*ka;
    
    d/dt(ACoundis) = AIn7undis*kTsi - ACoundis*SiCo - ACoundis*kTco; #what gets dissolved
    d/dt(ACo) = AIn7*kTsi + ACoundis*SiCo - ACo*kTco - ACo*ka;
    
    d/dt(AFe) = ACo*kTco + ACoundis*kTco;
    
    #Differential equations perfused tissues
    d/dt(ALi) = ka*(ASt+AIn1+AIn2+AIn3+AIn4+AIn5+AIn6+AIn7+ACo) + 
      QLi*(CBl - CvLi) - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));           #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides                                                          #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile   
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                               #(umol/h) rate of change in amount in the kidney
    d/dt(AUr) = CLRe*CKi_free;                                                    #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                 #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                 #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                 #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                       #(umol/h) rate of change in amount in the blood
    
    #Differential equations glucuronide
    d/dt(ALi.glu) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free)); #No EHC - (CLBi*CLi.glu); 
    #d/dt(AIn.glu) = (CLBi*CLi.glu) - ((Vmax_hyd * CIn.glu) / (Km_hyd + CIn.glu));
    
    
  }) 
  return(ModelZEN)
}

#*****************************************************************************#
#---For oral dose & enterohepatic circulation 
#with Clbi, Intestinal hydrolysis & 
#metabolism with the total liver concentration---------#

hydro_PO_EHC.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                      #Abbreviation (C for concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                              #ZEN concentration in intestine (umol/L)
    CIn.glu = AIn.glu / VIn;
    CLi = ALi / VLi;                                                              #Zen concentration in liver (umol/L)
    CLi.glu = ALi.glu / VLi;                                                      #ZEN-glu concentration in the liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN concentration in blood (umol/L)
    
    #Tissue venous concentrations (umol/L)                                        #Only for perfused organs; Equals tissue concentration/tissue:blood partition coefficient
    CvLi = CLi / PLib;                                                            #ZEN concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(ASt) = -(kTge*ASt);                                                      #(umol/h) rate of change of amount reaching the intestine (correcting for gastric emptying time)
    
    d/dt(AIn) = (kTge*ASt) + ((Vmax_hyd * CIn.glu) / (Km_hyd + CIn.glu)) - (ka*AIn) - (kTsi*AIn);              #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(AIn.glu) = (kBi*ALi.glu) - ((Vmax_hyd * CIn.glu) / (Km_hyd + CIn.glu))
    
    d/dt(ACo) = (kTsi*AIn) - (kTco*ACo);                                          #(umol/h) rate of change in amount of ZEN in colon
    d/dt(AFe) = kTco*ACo;                                                         #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi) / (Km_glu + CLi));                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi.glu) = ((Vmax_glu * CLi) / (Km_glu + CLi)) +
      - (kBi*ALi.glu);                                                           #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                               #(umol/h) rate of change in amount in the kidney
    d/dt(AUr) = CLRe*CKi_free;                                                                #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                 #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                 #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                 #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                       #(umol/h) rate of change in amount in the blood
    
  }) 
  return(ModelZEN)
}



#*****************************************************************************#
#---For oral dose & enterohepatic circulation - kBi---------------------------#

hydro_PO_EHC.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                      #Abbreviation (C for concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                              #ZEN concentration in intestine (umol/L)
    CIn.glu = AIn.glu / VIn;
    CLi = ALi / VLi;                                                              #Zen concentration in liver (umol/L)
    CLi.glu = ALi.glu / VLi;                                                      #ZEN-glu concentration in the liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN concentration in blood (umol/L)
    
    #Tissue venous concentrations (umol/L)                                        #Only for perfused organs; Equals tissue concentration/tissue:blood partition coefficient
    CvLi = CLi / PLib;                                                            #ZEN concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(ASt) = -(kTge*ASt);                                                      #(umol/h) rate of change of amount reaching the intestine (correcting for gastric emptying time)
    
    d/dt(AIn) = (kTge*ASt) + ((Vmax_hyd * CIn.glu) / (Km_hyd + CIn.glu)) - (ka*AIn) - (kTsi*AIn);              #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(AIn.glu) = (kBi*ALi.glu) - ((Vmax_hyd * CIn.glu) / (Km_hyd + CIn.glu))
    
    d/dt(ACo) = (kTsi*AIn) - (kTco*ACo);                                          #(umol/h) rate of change in amount of ZEN in colon
    d/dt(AFe) = kTco*ACo;                                                         #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi) / (Km_glu + CLi));                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi.glu) = ((Vmax_glu * CLi) / (Km_glu + CLi)) +
      - (kBi*ALi.glu);                                                           #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                               #(umol/h) rate of change in amount in the kidney
    d/dt(AKi_clearance) = CLRe*CKi_free;                                                                #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                 #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                 #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                 #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                       #(umol/h) rate of change in amount in the blood
    
  }) 
  return(ModelZEN)
}


#*****************************************************************************#
#---For oral dose & enterohepatic circulation---------------------------------#

PO_EHC.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                      #Abbreviation (C for concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                              #ZEN concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                              #Zen concentration in liver (umol/L)
    CLi.glu = ALi.glu / VLi;                                                      #ZEN-glu concentration in the liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN concentration in blood (umol/L)
    
    #Tissue venous concentrations (umol/L)                                        #Only for perfused organs; Equals tissue concentration/tissue:blood partition coefficient
    CvLi = CLi / PLib;                                                            #ZEN concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(ASt) = -(kTge*ASt);                                                      #(umol/h) rate of change of amount reaching the intestine (correcting for gastric emptying time)
    
    d/dt(AIn) = (kTge*ASt) + (CLBi*CLi.glu) - (ka*AIn) - (kTsi*AIn);              #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(ACo) = (kTsi*AIn) - (kTco*ACo);                                          #(umol/h) rate of change in amount of ZEN in colon
    d/dt(AFe) = kTco*ACo;                                                         #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi.glu) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free)) +
      - (CLBi*CLi.glu);                                                           #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                               #(umol/h) rate of change in amount in the kidney
    d/dt(AKi_clearance) = CLRe*CKi_free;                                                                #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                 #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                 #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                 #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                       #(umol/h) rate of change in amount in the blood
    
  }) 
  return(ModelZEN)
}



#*****************************************************************************#
#---For oral dose; without enterohepatic circulation--------------------------#

PO.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                      #Abbreviation (C for concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                              #ZEN concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                              #Zen concentration in liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN concentration in blood (umol/L)
    
    #Tissue venous concentrations (umol/L)                                        #Only for perfused organs; Equals tissue concentration/tissue:blood partition coefficient
    CvLi = CLi / PLib;                                                            #ZEN concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(ASt) = -(kTge*ASt);                                                      #(umol/h) rate of change of amount reaching the intestine (correcting for gastric emptying time)
    
    d/dt(AIn) = (kTge*ASt) - (ka*AIn) - (kTsi*AIn);                               #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(ACo) = (kTsi*AIn) - (kTco*ACo);                                          #(umol/h) rate of change in amount of ZEN in colon
    d/dt(AFe) = kTco*ACo;                                                         #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi_clearance) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free))                                                           #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                               #(umol/h) rate of change in amount in the kidney
    d/dt(AKi_clearance) = CLRe*CKi_free;                                                                #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                 #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                 #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                 #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                       #(umol/h) rate of change in amount in the blood
    
  }) 
  return(ModelZEN)
}


#*****************************************************************************#
#---For IV & enterohepatic circulation----------------------------------------#

IV_ECH.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular concentrations (umol/L)                                        #Abbreviation (C for concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                                #ZEN concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                                #Zen concentration in liver (umol/L)
    CLi.glu = ALi.glu / VLi;                                                        #ZEN-glu concentration in the liver (umol/L)
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
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(AIn) = (CLBi*CLi.glu) - (ka*AIn) - (kTsi*AIn);                             #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(ACo) = (kTsi*AIn) - (kTco*ACo);                                            #(umol/h) rate of change in amount of ZEN in colon
    d/dt(AFe) = kTco*ACo;                                                           #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));                              #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi.glu) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free)) +
      - (CLBi*CLi.glu);                                                             #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
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


#*****************************************************************************#
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
