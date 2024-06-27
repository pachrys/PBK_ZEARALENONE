#*****************************************************************************#
#---Model functions-----------------------------------------------------------#
#*****************************************************************************#

#*****************************************************************************#
#---For oral dose & enterohepatic circulation---------------------------------#

LInmplexModel <- function() {
  ModelZEN <- rxode2({  
    
    # Zearalenone
    ##Non perfused organ Concentrations (umol/L)
    CBi = ABi / VBi; #bile
    
    ##Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                              #ZEN Concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                              #Zen Concentration in liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN Concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN Concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN Concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN Concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN Concentration in blood (umol/L)
    
    ##Tissue venous Concentrations (umol/L)                                        #Only for perfused organs; Equals tissue Concentration/tissue:blood partition LInefficient
    CvIn = CIn / PInb;
    CvLi = CLi / PLib;                                                            #ZEN Concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN Concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN Concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN Concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN Concentration in rapidly perfused venous blood (umol/L)
    
    ##Tissue free Concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    # Zearalenone glucuronide
    ##Non perfused organ concentrations (umol/L)
    CInL.glu = AInL.glu / VInL; #Concentration intestinal lumen glucuronide
    CBi.glu = ABi.glu / VBi; #Concentration bile #need to find the bile volume
    
    
    ##Tissue cellular Concentrations (umol/L)                                      #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn.glu = AIn.glu / VIn; #Concentration intestinal tissue glucuronide
    CLi.glu = ALi.glu / VLi; #Concentration liver 
    CKi.glu = AKi.glu / VKi; #Concentration kidney glucuronide
    CAd.glu = AAd.glu / VAd; #Concentration adipose tissue glucuronide
    CSl.glu = ASl.glu / VSl; #Concentration slowly perfused tissue glucuronide
    CRa.glu = ARa.glu / VRa; #Concentration rapidly perfused tissue glucuronide
    CBl.glu = ABl.glu / VBl; #Concentration blood glucuronide
    
    ##Tissue venous Concentrations (umol/L)                                        #Only for perfused organs; Equals tissue Concentration/tissue:blood partition LInefficient
    CvIn.glu = CIn.glu / PInb.glu #Concentration leaving intestinal tissue glucuronide
    CvLi.glu = CLi.glu / PLib.glu #Concentration leaving liver tissue glucuronide
    CvKi.glu = CKi.glu / PKib.glu #Concentration leaving kidney glucuronide
    CvAd.glu = CAd.glu / PAdb.glu #Concentration leaving adipose tissue glucuronide
    CvSl.glu = CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue glucuronide
    CvRa.glu = CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue glucuronide
    
    # Active tnsport
    GlucExcretion.InL = VmaxExcIn.glu*CIn.glu/(KmExc.glu+CIn.glu); #excretion from the intestinal tissue to intestinal lumen
    InL.GlucUptake = VmaxAbsIn.glu*CInL.glu/(KmAbs.glu+CInL.glu); #uptake from the small intestinal lumen to the tissue
    
    # Metabolism
    InL.Gluc.Hydrolysis = (1-kTIn)*AInL.glu; #hydrolysis glucuronide to zen in small intestinal lumen
    Intestinal.GlucMet = CLint*CIn; #intestinal glucuronidation of zearalenone
    Liver.GlucMet = Vmax_glu*CLi/(Km_glu+CLi); #liver glucuronidation of zearalenone
    
    # Differential equations Zearalenone; rate of change in amounts in each tissue
    
    # Intestinal Lumen
    d/dt(AInL) = InL.Gluc.Hydrolysis - (ka*AInL) - (kTIn*AInL); 

    # Intestinal tissue
    d/dt(AIn) = (ka*AInL) + (QIn*CBl) - (QPo*CvIn) - Intestinal.GlucMet;  

    # Feces
    d/dt(AFe) = (kTIn*AInL); 
    
    # Liver
    d/dt(ALi) = (QPo*CvIn) + (QLi*(CBl-CvLi)) - Liver.GlucMet; 
    
    # Kidney
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free); #CLRe is equal to rat GER
    d/dt(AUr) = (CLRe*CKi_free); 
    
    # Adipose
    d/dt(AAd) = QAd*(CBl - CvAd); 
    
    # Slowly Perfused Tissue
    d/dt(ASl) = QSl*(CBl - CvSl); 
    
    # Rapidly Perfused Tissue
    d/dt(ARa) = QRa*(CBl - CvRa); 
    
    # Blood
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) + (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) + 
      -((QIn + QLi + QKi + QAd + QSl + QRa)*CBl); 
    
    
    ## Differential equations Glucuronide
    d/dt(AIn.glu) = Intestinal.GlucMet - GlucExcretion.InL + InL.GlucUptake + (QIn*CBl.glu) - (QPo*CvIn.glu); #(umol/h)
    #               fromIntZen           toInL               fromInL          fromblood       toportal
    
    d/dt(AInL.glu) = GlucExcretion.InL + (kTge*CBi.glu) - InL.GlucUptake - InL.Gluc.Hydrolysis - (kTIn*AInL.glu); #(umol/h) 
    #                 fromIntissue        fromBile        toSItissue       hydrolysedtoZEN        
    d/dt(AFe.glu) = (kTIn*AInL.glu)
    #               fromInLumen
    
    # Liver
    d/dt(ALi.glu) = (QPo*CvIn.glu) + QLi*(CBl-CvLi.glu) + Liver.GlucMet - (CBi*CLi.glu);
    #               fromPortal       fromToBlood          fromZEN         toBile
    d/dt(ABi.glu) = (CBi*CLi.glu) - (kTge*CBi.glu);
    #               fromLiver       toSmallIntestine
    
    # Kidney
    d/dt(AKi.glu) = QKi*(CBl.glu-CvKi.glu) - (CLRe*CKi.glu); #CLRe is equal to rat GER
    #               fromToBlood              toUrine
    d/dt(AUr.glu) = (CLRe*CKi.glu); 
    
    # Adipose
    d/dt(AAd.glu) = QAd*(CBl.glu - CvAd.glu); 
    # Slowly Perfused Tissue
    d/dt(ASl.glu) = QSl*(CBl.glu - CvSl.glu); 
    # Rapidly Perfused Tissue
    d/dt(ARa.glu) = QRa*(CBl.glu - CvRa.glu); 
    # Blood
    d/dt(ABl.glu) = (QLi*CvLi) + (QKi*CvKi.glu) + (QAd*CvAd.glu) + (QSl*CvSl.glu) + (QRa*CvRa.glu) + 
      -((QIn + QLi + QKi + QAd + QSl + QRa)*CBl.glu); 
  }) 
  return(ModelZEN)
}


#*****************************************************************************#
#---For oral dose & enterohepatic circulation---------------------------------#

PO_EHC.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular Concentrations (umol/L)                                      #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                              #ZEN Concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                              #Zen Concentration in liver (umol/L)
    CLi.glu = ALi.glu / VLi;                                                      #ZEN-glu Concentration in the liver (umol/L)
    CKi= AKi / VKi;                                                               #ZEN Concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                              #ZEN Concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                              #ZEN Concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                              #ZEN Concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                              #ZEN Concentration in blood (umol/L)
    
    #Tissue venous Concentrations (umol/L)                                        #Only for perfused organs; Equals tissue Concentration/tissue:blood partition LInefficient
    CvLi = CLi / PLib;                                                            #ZEN Concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                            #ZEN Concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                            #ZEN Concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                            #ZEN Concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                            #ZEN Concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free Concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(ASt) = -(kTge*ASt);                                                      #(umol/h) rate of change of amount reaching the intestine (LInrrecting for gastric emptying time)
    
    d/dt(AIn) = (kTge*ASt) + (kBi*ALi.glu) - (ka*AIn) - (kTSIn*AIn);              #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(ALIn) = (kTSIn*AIn) - (kTLIn*ALIn);                                          #(umol/h) rate of change in amount of ZEN in LInlon
    d/dt(AFe) = kTLIn*ALIn;                                                         #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi.glu) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free)) +
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
#---For IV & enterohepatic circulation----------------------------------------#

IV_ECH.PBKfunction <- function() {
  ModelZEN <- rxode2({  
    
    #Tissue cellular Concentrations (umol/L)                                        #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn;                                                                #ZEN Concentration in intestine (umol/L)
    CLi = ALi / VLi;                                                                #Zen Concentration in liver (umol/L)
    CLi.glu = ALi.glu / VLi;                                                        #ZEN-glu Concentration in the liver (umol/L)
    CKi = AKi / VKi;                                                                #ZEN Concentration in the kidney (umol/L)
    CAd = AAd / VAd;                                                                #ZEN Concentration in adipose (umol/L)
    CSl = ASl / VSl;                                                                #ZEN Concentration in slowly perfused tissue (umol/L)
    CRa = ARa / VRa;                                                                #ZEN Concentration in rapidly perfused tissue (umol/L)
    CBl = ABl / VBl;                                                                #ZEN Concentration in blood (umol/L)
    
    #Tissue venous Concentrations (umol/L)                                          #Only for perfused organs; Equals tissue Concentration/tissue:blood partition LInefficient
    CvLi = CLi / PLib;                                                              #ZEN Concentration in liver venous blood (umol/L)
    CvKi = CKi / PKib;                                                              #ZEN Concentration in kidney venous blood (umol/L)
    CvAd = CAd / PAdb;                                                              #ZEN Concentration in adipose venous blood (umol/L)
    CvSl = CSl / PSlb;                                                              #ZEN Concentration in slowly perfused venous blood (umol/L)
    CvRa = CRa / PRab;                                                              #ZEN Concentration in rapidly perfused venous blood (umol/L)
    
    #Tissue free Concentrations (umol/L)
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    #Differential equations gastro-intestinal tract - not perfused
    d/dt(AIn) = (CLBi*CLi.glu) - (ka*AIn) - (kTSIn*AIn);                             #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    d/dt(ALIn) = (kTSIn*AIn) - (kTLIn*ALIn);                                            #(umol/h) rate of change in amount of ZEN in LInlon
    d/dt(AFe) = kTLIn*ALIn;                                                           #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    d/dt(ALi) = (ka*AIn) + QLi*(CBl - CvLi) +
      - ((Vmax_glu * CLi_free) / (Km_glu + CLi_free));                              #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    d/dt(ALi.glu) = ((Vmax_glu * CLi_free) / (Km_glu + CLi_free)) +
      - (CLBi*CLi.glu);                                                             #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    d/dt(AKi) = QKi*(CBl - CvKi) - (CLRe*CKi_free);                                 #(umol/h) rate of change in amount in the kidney
    d/dt(AUr) = CLRe*CKi_free;                                            #(umol/h) rate of change in amount ZEN eliminated in urine
    
    d/dt(AAd) = QAd*(CBl - CvAd);                                                   #(umol/h) rate of change in amount in the adipose tissue
    d/dt(ASl) = QSl*(CBl - CvSl);                                                   #(umol/h) rate of change in amount in the slowly perfused
    d/dt(ARa) = QRa*(CBl - CvRa);                                                   #(umol/h) rate of change in amount in the rapidly perfused
    
    d/dt(ABl) = (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl);                                         #(umol/h) rate of change in amount in the blood
    
  }) 
  return(ModelZEN)
}

