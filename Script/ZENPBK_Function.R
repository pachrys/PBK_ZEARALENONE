# Libraries-------------------------------------------------------------------#
library(deSolve)

#---Model functions-----------------------------------------------------------#
PBK_ZENandGLUC <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #----Concentrations-------------------------------------------------------#
    
    # Zearalenone
    
    ## Non perfused organ Concentrations (umol/L)
    CInL = AInL / VInL #intestinal lumen
    # CBi = ABi / VBi #bile 
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn = AIn / VIn #intestine                                                             
    CLi = ALi / VLi #liver                                                              
    CKi = AKi / VKi #kidney                                                               
    CAd = AAd / VAd #adipose                                                              
    CSl = ASl / VSl #slowly perfused                                                     
    CRa = ARa / VRa #rapidly perfused                                                              
    CBl = ABl / VBl #blood                                                              
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs; Equals tissue Concentration/tissue:blood partition coeefficient
    CvIn = CIn / PInb #intestinal tissue
    CvLi = CLi / PLib #liver                                                            
    CvKi = CKi / PKib #kidney                                                            
    CvAd = CAd / PAdb #adipose                                                             
    CvSl = CSl / PSlb #slowly perfused                                                            
    CvRa = CRa / PRab #rapidly perfused                                                           
    
    ## Tissue free Concentrations (umol/L), used for clearance/transport
    CKi_free = CvKi*fub
    CLi_free = CvLi*fub
    
    
    
    # Zearalenone glucuronide
    
    ## Non perfused organ concentrations (umol/L)
    CInL.glu = AInL.glu / VInL #intestinal lumen glucuronide
    CBi.glu = ABi.glu / VBi # bile 
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn.glu = AIn.glu / VIn #intestinal tissue 
    CLi.glu = ALi.glu / VLi #liver 
    CKi.glu = AKi.glu / VKi #kidney 
    CAd.glu = AAd.glu / VAd #adipose tissue 
    CSl.glu = ASl.glu / VSl #slowly perfused tissue 
    CRa.glu = ARa.glu / VRa #rapidly perfused tissue 
    CBl.glu = ABl.glu / VBl #blood 
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs; Equals tissue Concentration/tissue:blood partition LInefficient
    CvIn.glu = CIn.glu / PInb.glu #Concentration leaving intestinal tissue 
    CvLi.glu = CLi.glu / PLib.glu #Concentration leaving liver  
    CvKi.glu = CKi.glu / PKib.glu #Concentration leaving kidney 
    CvAd.glu = CAd.glu / PAdb.glu #Concentration leaving adipose  
    CvSl.glu = CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue 
    CvRa.glu = CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue 
    
    
    
    #----Differential equations-------------------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    #dD = -D
    
    # Intestinal Lumen
    dAInL = - ka*AInL - kTIn*AInL + (1-kTIn)*AInL.glu #+D
    
    # Intestinal tissue
    dAIn = ka*AInL + QIn*CBl - QPo*CvIn - VmaxMet.In*CIn/(KmGluc.In+CIn)  
    
    # Feces
    dAFe = kTIn*AInL
    
    # Liver
    dALi = QPo*CvIn + QLi*(CBl-CvLi) - VmaxMet.Li*CLi/(KmGluc.Li+CLi) 
    
    # Kidney
    dAKi = QKi*(CBl - CvKi) - CL_Re*CKi_free #CL_Re is equal to rat GER
    dAUr = CL_Re*CKi_free
    
    # Adipose
    dAAd = QAd*(CBl - CvAd) 
    
    # Slowly Perfused Tissue
    dASl = QSl*(CBl - CvSl) 
    
    # Rapidly Perfused Tissue
    dARa = QRa*(CBl - CvRa) 
    
    # Blood
    dABl = QLi*CvLi + QKi*CvKi + QAd*CvAd + QSl*CvSl + QRa*CvRa + 
      -(QIn + QLi + QKi + QAd + QSl + QRa)*CBl 
    
    
    ## Differential equations Glucuronide
    dAIn.glu = VmaxMet.In*CIn/(KmGluc.In+CIn) - VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) + VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu) + QIn*CBl.glu - QPo*CvIn.glu #(umol/h)
    #          fromIntZen                       toInL                                       fromInL                                       fromblood     toportal
    
    dAInL.glu = VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) + kTge*CBi.glu - VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu) - (1-kTIn)*AInL.glu - kTIn*AInL.glu #(umol/h) 
    #           fromIntissue                                  fromBile       toSItissue                                    hydrolysedtoZEN        
    dAFe.glu = kTIn*AInL.glu
    #         fromInLumen
    
    # Liver
    dALi.glu = QPo*CvIn.glu + QHa*CBl -QLi*CvLi.glu + VmaxMet.Li*CLi/(KmGluc.Li+CLi) - CL_Bi*CLi.glu;
    #          fromPortal     fromToBlood             fromZEN                          toBile
    # Bile
    dABi.glu = CL_Bi*CLi.glu - kTge*CBi.glu
    #          fromLiver       toSmallIntestine
    
    # Kidney
    dAKi.glu = QKi*(CBl.glu-CvKi.glu) - CL_Re*CKi.glu; #CL_Re is equal to rat GER
    #          fromToBlood              toUrine
    dAUr.glu = CL_Re*CKi.glu
    
    # Adipose
    dAAd.glu = QAd*(CBl.glu - CvAd.glu) 
    # Slowly Perfused Tissue
    dASl.glu = QSl*(CBl.glu - CvSl.glu) 
    # Rapidly Perfused Tissue
    dARa.glu = QRa*(CBl.glu - CvRa.glu) 
    # Blood
    dABl.glu = QLi*CvLi + QKi*CvKi.glu + QAd*CvAd.glu + QSl*CvSl.glu + QRa*CvRa.glu + 
      - (QIn + QLi + QKi + QAd + QSl + QRa)*CBl.glu 
    
    #it gives bback the amount in each organ but also the total for mass balance calculations
    return(list(c(#dD,
                  dAInL, dAIn, dAFe, dALi, dAKi, dAUr, dAAd, dASl, dARa, dABl,
                  dAIn.glu, dAInL.glu, dAFe.glu, ALi.glu, dABi.glu, dAKi.glu, 
                  dAUr.glu, dAAd.glu, dASl.glu, dARa.glu, dABl.glu
           )))
  })
}
