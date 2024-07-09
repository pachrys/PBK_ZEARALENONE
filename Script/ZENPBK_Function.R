# Libraries-------------------------------------------------------------------#
library(deSolve)

#---Model functions-----------------------------------------------------------#

PBK_ZENandGLUC <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #----Concentrations-------------------------------------------------------#
    
    # Zearalenone
    
    ## Non perfused organ Concentrations (umol/L)
    CInL <- AInL / VInL #intestinal lumen
    # CBi <- ABi / VBi #bile 
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn <- AIn / VIn #intestine                                                             
    # CLi <- ALi / VLi #liver                                                              
    # CKi <- AKi / VKi #kidney                                                               
    # CAd <- AAd / VAd #adipose                                                              
    # CSl <- ASl / VSl #slowly perfused                                                     
    # CRa <- ARa / VRa #rapidly perfused                                                              
    # CBl <- ABl / VBl #blood                                                              
    # 
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition coeefficient
    # CvIn <- CIn / PInb #intestinal tissue
    # CvLi <- CLi / PLib #liver                                                            
    # CvKi <- CKi / PKib #kidney                                                            
    # CvAd <- CAd / PAdb #adipose                                                             
    # CvSl <- CSl / PSlb #slowly perfused                                                            
    # CvRa <- CRa / PRab #rapidly perfused                                                           
    # 
    # ## Tissue free Concentrations (umol/L), used for clearance/transport
    # CKi_free <- CvKi*fub
    # CLi_free <- CvLi*fub
    # 
    # 
    
    # Zearalenone glucuronide
    
    # ## Non perfused organ concentrations (umol/L)
    CInL.glu <- AInL.glu / VInL #intestinal lumen glucuronide
    # CBi.glu <- ABi.glu / VBi # bile 
    # 
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn.glu <- AIn.glu / VIn #intestinal tissue
    # CLi.glu <- ALi.glu / VLi #liver 
    # CKi.glu <- AKi.glu / VKi #kidney 
    # CAd.glu <- AAd.glu / VAd #adipose tissue 
    # CSl.glu <- ASl.glu / VSl #slowly perfused tissue
    # CRa.glu <- ARa.glu / VRa #rapidly perfused tissue 
    # CBl.glu <- ABl.glu / VBl #blood
    # 
    # ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition LInefficient
    # CvIn.glu <- CIn.glu / PInb.glu #Concentration leaving intestinal tissue 
    # CvLi.glu <- CLi.glu / PLib.glu #Concentration leaving liver  
    # CvKi.glu <- CKi.glu / PKib.glu #Concentration leaving kidney 
    # CvAd.glu <- CAd.glu / PAdb.glu #Concentration leaving adipose  
    # CvSl.glu <- CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue 
    # CvRa.glu <- CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue
    # 
    
    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    
    # Zearalenone-------------------- ##
    # Intestinal Lumen
    dAInL <- - ka*AInL - kTIn*AInL + ((1-kTIn))*AInL.glu #- VmaxAbs*CInL/(KmAbs+CInL) + VmaxSec*CIn/(KmSec+CIn) 
    #          toInt              
    
    # Intestinal tissue
    dAIn <-  + ka*AInL - VmaxMet.In*CIn/(KmGluc.In+CIn) #+ VmaxAbs*CInL/(KmAbs+CInL) - VmaxSec*CIn/(KmSec+CIn) + QIn*CBl - QPo*CvIn  
    #        fromLumenActive      toglucuronide                     fromBlood toLiver             toLumen                   fromLum                         
    
    # # Feces
    dAFe <- kTIn*AInL
    #       fromIntestine
    
    # # Liver
    # dALi <- QPo*CvIn + QHa*CBl - QLi*CvLi - VmaxMet.Li*CLi/(KmGluc.Li+CLi) + ka*AIn
    # #       fromIntest fromBlo   toBlood    toglucuronide
    # # Kidney
    # dAKi <- QKi*(CBl - CvKi) - CL_Re*CKi_free #CL_Re is equal to rat GER
    # #       fromtoblood        toUrine
    # dAUr <- CL_Re*CKi_free
    # #       fromKidney
    # # Adipose
    # dAAd <- QAd*(CBl - CvAd) 
    # #       fromtoblood
    # # Slowly Perfused Tissue
    # dASl <- QSl*(CBl - CvSl) 
    # #       fromtoslowly
    # # Rapidly Perfused Tissue
    # dARa <- QRa*(CBl - CvRa) 
    # #       fromtorapidly
    # # Blood
    # dABl <- QLi*CvLi + QKi*CvKi + QAd*CvAd + QSl*CvSl + QRa*CvRa + 
    # # from  liver      kidney     adipose    slowly     rapidly
    #   -(QIn + QHa + QKi + QAd + QSl + QRa)*CBl 
    # # to int  hepA  kid   adi   slo   rap  
    # 
    # ## Glucuronide---------------------- ##
    # 
    # Intestinal lumen
    dAInL.glu <- AIn.glu - ((1-kTIn))*AInL.glu #VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) - VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu) - (1-kTIn)*AInL.glu - kTIn*AInL.glu + kTge*CBi.glu  #(umol/h)
    #            fromIntissue                                  hydrolysedtoZEN   toSItissue                                      hydrolysedtoZEN     tofeces        # fromBile
    
    # Intestine
    dAIn.glu <- VmaxMet.In*CIn/(KmGluc.In+CIn) - AIn.glu #VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) #+ VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu) + QIn*CBl.glu - QPo*CvIn.glu #(umol/h)
    #           fromIntZen                       toInL                                         fromInL                                         fromblood     toportal
    
    # # Feces
    # dAFe.glu <- kTIn*AInL.glu
    # #           fromInLumen
    # 
    # # Blood
    # dABl.glu <- QLi*CvLi.glu + QKi*CvKi.glu + QAd*CvAd.glu + QSl*CvSl.glu + QRa*CvRa.glu +
    #   - ((QIn + QHa + QKi + QAd + QSl + QRa) * CBl.glu) 
    # 
    # # Liver
    # dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + QPo*CvIn.glu - QLi*CvLi.glu + QHa*CBl.glu - CL_Bi*CLi.glu
    # #           fromZen                          fromPortal     toBlood        fromHepaticA  toBile
    # 
    # # Bile
    # dABi.glu <- CL_Bi*CLi.glu - kTge*CBi.glu
    # #           fromLiver       toIntestLumen
    # 
    # # Kidney
    # dAKi.glu <- QKi*(CBl.glu-CvKi.glu) - CL_Re*CKi.glu
    # #           fromToBlood              toUrine
    # 
    # # Urine
    # dAUr.glu <- CL_Re*CKi.glu
    # #           fromKidney
    # 
    # # Adipose
    # dAAd.glu <- QAd*(CBl.glu - CvAd.glu) 
    # 
    # # Slowly Perfused Tissue
    # dASl.glu <- QSl*(CBl.glu - CvSl.glu)
    # 
    # # Rapidly Perfused Tissue
    # dARa.glu <- QRa*(CBl.glu - CvRa.glu)
    # 
    
    #----Output---------------------------------------------------------------#
    
    return(list(c(dAInL, dAIn, dAFe,
                  # dALi, dAKi, dAUr, dAAd, dASl, dARa, dABl,
                  dAIn.glu, dAInL.glu
                  #dAFe.glu,
                  # dABl.glu, dALi.glu, dABi.glu, dAKi.glu, dAUr.glu, dAAd.glu,
                  # dASl.glu, dARa.glu
                  )))
  })
}


# Tests Archive

# # Intestinal Lumen As Mukherjee
# dAInL <- - (ka.fa + ka.sl)*AInL + (1-kTIn)*AInL.glu
# #          toIntest  toFeces      fromGluhydro
# 
# # Fast Lumen As Mukherjee
# dAInL.f <- ka.fa*AInL - ka.fit*AInL.f
# 
# # Slow Lumen As Mukherjee
# dAInL.s <- ka.sl*AInL - kex*AInL.s
# 
# # Intestinal tissue As Mukherjee
# dAIn <- ka.fit*AInL.f + QIn*CBl - QPo*CvIn - VmaxMet.In*CIn/(KmGluc.In+CIn)
# #       fromLumen fromBlood toLiver    toglucuronide
# 
# # Feces As Mukherjee
# dAFe <- kex*AInL.s
# #       fromIntestine