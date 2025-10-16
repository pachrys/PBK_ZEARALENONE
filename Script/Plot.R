#I run the model with 7 days

# Exposure:
E.df <- data.frame(
  What = c("NOEL", "TDI", 
           
           "Single_ZEN.Diana et al.", "Single_ZEN.Ehrlich et al.", 
           "Single_ZEN.Frizzell et al.", "Single_ZEN.Molina-Molina et al.",  
           "Single_ZEN.Minervini et al.", "Single_ZEN.Malekinejat et al", 
           "Single_ZEN.Le Guevel and Pakdel",
           
           "RPF_ZEN.Diana et al.", "RPF_ZEN.Ehrlich et al.", 
           "RPF_ZEN.Frizzell et al.", "RPF_ZEN.Molina-Molina et al.",  
           "RPF_ZEN.Minervini et al.", "RPF_ZEN.Malekinejat et al", 
           "RPF_ZEN.Le Guevel and Pakdel"
           ), 
  Exposure = c(10, 0.25, 
               
               222, 69,
               226, 530, 
               45, 230, 
               9,
               
               52, 16,
               53, 123,
               10.5, 53.7,
               2.1), # ug/kg bw/dau
  Single_Zearalenone = c(0.07042453, 0.001760176, 
                  
                  1.571858, 0.4866583, 
                  1.600342, 3.781942, 
                  0.3171924, 1.628832, 
                  0.06338047,
                  
                  0.3665986, 0.1126964,
                  0.3736581, 0.8828815,
                  0.07465114, 0.3785999,
                  0.01478618), 
  # aZearalenol = c(0.004207185, 0.0001051519, 
  #                 
  #                 0.09392983, 0.02907569, 
  #                 0.09563237, 0.2260515, 
  #                 0.01895019, 0.09733526, 
  #                 0.003786364, 
  #                 
  #                 0.02190212, 0.006732585,
  #                 0.02232391, 0.05275234,
  #                 0.004459688, 0.02261918,
  #                 0.0008833205), #nM
  Zearalenone_Equivalent =  c((0.07042453 + (0.004207185*RPF)),(0.001760176 + (0.0001051519*RPF)),
                              
                              (1.571858 + (0.09392983*RPF)), (0.4866583 + (0.02907569*RPF)),
                              (1.600342 + (0.09563237*RPF)), (3.781942+(0.2260515*RPF)), 
                              (0.3171924+(0.01895019*RPF)), (1.628832+(0.09733526*RPF)), 
                              (0.06338047+(0.003786364*RPF)),
                              
                              1.571215,0.4829886,
                              1.601473, 3.78426,
                              0.319934, 1.622655,
                              0.06336881), #nM
  Group = c("NOEL", "TDI", 
            
            "EC50", "EC50", 
            "EC50", "EC50",  
            "EC50", "EC50", 
            "EC50",
            
            "EC50", "EC50", 
            "EC50", "EC50",  
            "EC50", "EC50", 
            "EC50"
    
  )
)




E.df$Group <- ifelse(grepl("^Single_ZEN", df$What), "Single_ZEN", df$What)




# Clean
E.df <- data.frame(
  What = c("Single_NOEL", "Single_TDI", 
           "RPF_NOEL", "RPF_TDI",
           
           "Single_ZEN.Diana et al.", "Single_ZEN.Ehrlich et al.", 
           "Single_ZEN.Frizzell et al.", "Single_ZEN.Molina-Molina et al.",  
           "Single_ZEN.Minervini et al.", "Single_ZEN.Malekinejat et al", 
           "Single_ZEN.Le Guevel and Pakdel",
           
           "RPF_ZEN.Diana et al.", "RPF_ZEN.Ehrlich et al.", 
           "RPF_ZEN.Frizzell et al.", "RPF_ZEN.Molina-Molina et al.",  
           "RPF_ZEN.Minervini et al.", "RPF_ZEN.Malekinejat et al", 
           "RPF_ZEN.Le Guevel and Pakdel"
  ), 
  Method = c("Simple", "Simple", 
             "RPF", "RPF",
             
             "Simple", "Simple", 
             "Simple", "Simple",   
             "Simple", "Simple", 
             "Simple",
             
             "RPF", "RPF", 
             "RPF", "RPF",  
             "RPF", "RPF", 
             "RPF"
  ),
  Exposure = c(10, 0.25, 
               10, 0.26,
               
               222, 69,
               226, 530, 
               45, 230, 
               9,
               
               52, 16,
               53, 123,
               10.5, 53.7,
               2.1), # ug/kg bw/dau
  Zearalenone = c(0.07042453, 0.001760176, 
                         0.07042453, 0.001760176,
                         
                         1.571858, 0.4866583, 
                         1.600342, 3.781942, 
                         0.3171924, 1.628832, 
                         0.06338047,
                         
                         0.3665986, 0.1126964,
                         0.3736581, 0.8828815,
                         0.07465114, 0.3785999,
                         0.01478618), 
  aZearalenol = c(0.004207185, 0.0001051519,
                  0.004207185, 0.0001051519,

                  0.09392983, 0.02907569,
                  0.09563237, 0.2260515,
                  0.01895019, 0.09733526,
                  0.003786364,

                  0.02190212, 0.006732585,
                  0.02232391, 0.05275234,
                  0.004459688, 0.02261918,
                  0.0008833205), #nM
  Zearalenone_Equivalent =  c((0.07042453 + (0.004207185*RPF)),(0.001760176 + (0.0001051519*RPF)),
                              (0.07042453 + (0.004207185*RPF)),(0.001760176 + (0.0001051519*RPF)),
                              
                              (1.571858 + (0.09392983*RPF)), (0.4866583 + (0.02907569*RPF)),
                              (1.600342 + (0.09563237*RPF)), (3.781942+(0.2260515*RPF)), 
                              (0.3171924+(0.01895019*RPF)), (1.628832+(0.09733526*RPF)), 
                              (0.06338047+(0.003786364*RPF)),
                              
                              1.571215,0.4829886,
                              1.601473, 3.78426,
                              0.319934, 1.622655,
                              0.06336881), #nM
  Concentration = c(0.07042453, 0.001760176,
                    (0.07042453 + (0.004207185*RPF)), (0.001760176 + (0.0001051519*RPF)),

                    1.571858, 0.4866583,
                    1.600342, 3.781942,
                    0.3171924, 1.628832,
                    0.06338047,

                    1.571858, 0.4866583,
                    1.600342, 3.781942,
                    0.3171924, 1.628832,
                    0.06338047),
  Value = c("NOEL", "TDI", # "HBGV", "HBGV",
            "NOEL", "TDI", # "HBGV", "HBGV", 
            
            "EC50", "EC50", 
            "EC50", "EC50",  
            "EC50", "EC50", 
            "EC50",
            
            "EC50", "EC50", 
            "EC50", "EC50",  
            "EC50", "EC50", 
            "EC50"
            
  )
)


# df$Group <- ifelse(grepl("^EC50_ZEN", df$What), "EC50_ZEN", df$What)

E.df_long <- E.df %>% pivot_longer(Zearalenone:Zearalenone_Equivalent, names_to = "Mycotoxin", values_to = "Concentration")
library(scales)


ggplot(E.df_long, aes(x = Value, y = log10(Concentration), color = Value, fill = Value)) +
  geom_boxplot() +
  facet_grid(Mycotoxin ~ ., scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
    x = "Exposure (μg/kg bw)",
    y = "Concentration (log10(nM))") +
  coord_flip() +
  # theme_bw() +
  CP_theme +
  theme(legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text.y = element_text(size = 30, angle = 0, face = "bold")  # Use strip.text.y for flipped plots
        # strip.text = element_text(size = 30, angle = 0, hjust = 1, face = "bold")
  )
ggsave(
  filename = here(OUTPUT,"RPF.or.notConcentrations.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)

ggplot(E.df_long, aes(x = Value, y = Concentration, color = Value, fill = Value)) +
  geom_boxplot() +
  facet_grid(Mycotoxin ~ ., scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
    x = "Exposure (μg/kg bw)",
    y = "Concentration (nM)") +
  coord_flip() +
  # theme_bw() +
  CP_theme +
  theme(legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 30, angle = 90, hjust = 0, face = "bold")
  )
ggsave(
  filename = here(OUTPUT,"RPF.or.not.nonlog.Concentrations.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)

Plot.df <- E.df %>% select(Method, Exposure, Concentration, Value)
# Plot.df$Method <- factor(Plot.df$Method, levels = c("TDI", "NOEL", "RPF", "Simple"))


ggplot(Plot.df, aes(x = Value, y = Exposure, colour = Value, fill = Value)) +
  geom_boxplot() +
  # facet_grid(Method ~ .) + #, scales = "free_x", space = "free_x"
  facet_grid(Method ~., labeller = labeller(Method = c("RPF" = "RPF (Zearalenone equivalent)", "Simple" = "Zearalenone"))) +
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
    x = "",
    y = "Exposure (μg/kg bw)") +
  coord_flip() +
  # theme_bw() +
  # theme(legend.position = "none") +
  CP_theme +
  theme( 
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text.y = element_text(size = 30, angle = 0, face = "bold")  # Use strip.text.y for flipped plots
  )
# print(Plot.C.Pre)
ggsave(
  filename = here(OUTPUT,"RPF.or.not.Exposure.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)


ggplot(Plot.df, aes(x = Method, y = Exposure, colour = Method, fill = Method)) +
  geom_boxplot() +
  # facet_grid(Method ~ .) + #, scales = "free_x", space = "free_x"
  # facet_grid(Method ~., labeller = labeller(Method = c("RPF" = "RPF (Zearalenone equivalent)", "Simple" = "Zearalenone"))) +
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
    x = "Exposure (μg/kg bw)",
    y = "Concentration (uM)") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")


ggplot(Plot.df, aes(x = Group, y = Exposure, colour = Group, fill = Group)) +
  geom_boxplot() +
  facet_grid(Method ~ .) + #, scales = "free_x", space = "free_x"
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
    x = "",
    y = "Concentration (uM)") +
  coord_flip() +
  theme_bw()

ggplot(Plot.df, aes(x = Method, y = Exposure, colour = Method, fill = Method)) +
  geom_boxplot() +
  # facet_grid(Group ~ .) + #, scales = "free_x", space = "free_x"
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Pastel2") +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
       x = "",
       y = "Exposure (μg/kg bw)") +
  coord_flip() +
  theme_bw()
  
  CP_theme +
  theme(legend.position = "bottom",  
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 30, angle = 90, hjust = 0, face = "bold")
  )
# print(Plot.C.Pre)
ggsave(
  filename = here(OUTPUT,"EC50.ZENandaZEL.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)



# C_ZEN_eq = C_ZEN + (C_aZEL * 55)
# print(C_ZEN_eq)
# Concentration:
RPF <- 55

df <- data.frame(
  What = c("NOEL", "TDI", "EC50_ZEN.Diana et al.", "EC50_ZEN.Ehrlich et al.", 
           "EC50_ZEN.Frizzell et al.", "EC50_ZEN.Molina-Molina et al.", "EC50_ZEN.Minervini et al.",
           "EC50_ZEN.Malekinejat et al", "EC50_ZEN.Le Guevel and Pakdel"), 
  Exposure = c(10, 0.25, 222, 69, 226, 530, 45, 230, 9), # ug/kg bw/dau
  Zearalenone = c(0.07042453, 0.001760176, 1.571858, 0.4866583, 1.600342, 3.781942, 0.3171924, 1.628832, 0.06338047), 
  aZearalenol = c(0.004207185, 0.0001051519, 0.09392983, 0.02907569, 0.09563237, 0.2260515, 0.01895019, 0.09733526, 0.003786364), #nM
  Zearalenone_Equivalent =  c((0.07042453 + (0.004207185*RPF)),(0.001760176 + (0.0001051519*RPF)), (1.571858 + (0.09392983*RPF)), (0.4866583 + (0.02907569*RPF)),
               (1.600342 + (0.09563237*RPF)), (3.781942+(0.2260515*RPF)), (0.3171924+(0.01895019*RPF)),
               (1.628832+(0.09733526*RPF)), (0.06338047+(0.003786364*RPF))) #nM
)

df$Group <- ifelse(grepl("^EC50_ZEN", df$What), "EC50_ZEN", df$What)

df_long <- df %>% pivot_longer(Zearalenone:Zearalenone_Equivalent, names_to = "Mycotoxin", values_to = "Concentration")

ggplot(df_long, aes(x = Group, y = log10(Concentration))) +
  geom_boxplot() +
  facet_grid(Mycotoxin ~ ., scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  # geom_hline(yintercept = NOEL, linetype = "dashed", color = "seagreen", linewidth = 0.5) +
  # geom_text(aes(x = 0.5, y = NOEL, label = "NOEL: 10 μg/kg bw/day"),
  #           hjust = -0.1, vjust = 0, color = "seagreen", size = 15) +
  labs(#title = "EC50 Values obtained from in vitro estrogenicity tests",
       x = "Exposure (μg/kg bw)",
       y = "Concentration (log10(nM))") +
  coord_flip() +
  theme_bw()


  CP_theme +
  theme(legend.position = "bottom",  
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 30, angle = 90, hjust = 0, face = "bold")
  )
# print(Plot.C.Pre)
ggsave(
  filename = here(OUTPUT,"EC50.ZENandaZEL.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)






# Exposure:
ZENec50_nM.df <- data.frame(
  Study = c("Diana et al.", "Ehrlich et al.", "Frizzell et al.", "Molina-Molina et al.", "Minervini et al.", "Malekinejat et al", "Le Guevel and Pakdel"), 
  EC50 =    c(1.58, 0.49, 1.60, 3.81, 0.31, 1.64, 0.06), 
  Chronic = c(222, 69, 226, 530, 45, 230, 9), #{ug/kg bw}
  Accute =  c(314, 121, 226, 929, 79, 402.5, 15.75), #{ug/kg bw}
  Mycotoxin = c("Zearalenone")
)




# C_ZEN <- max(C.Pre$Blood_ZEN) * 1000 #* 1000
# C_aZEL <- max(C.Pre$Blood_aZEL) * 1000 #* 1000

C_ZEN.TDI <- 0.0007671529
C_aZEL.TDI <- 4.313171e-05

C_ZEN.LB <- 7.365674e-06
C_aZEL.LB <- 4.133856e-07

C_ZEN.UB <- 8.898958e-05
C_aZEL.UB <- 4.994318e-06

ZENec50_nM.df <- data.frame(
  Study = c("Diana et al.", "Ehrlich et al.", "Frizzell et al.", "Molina-Molina et al.", "Minervini et al.", "Malekinejat et al", "Le Guevel and Pakdel"), 
  EC50 =              c(1.58, 0.49, 1.60, 3.81, 0.31, 1.64, 0.06), 
  Chronic = c(222, 69, 226, 530, 45, 230, 9), #{ug/kg bw}
  Accute =  c(314, 121, 226, 929, 79, 402.5, 15.75), #{ug/kg bw}
  Mycotoxin = c("Zearalenone")
)

aZELec50_nM.df <- data.frame(
  Study = c("Diana et al.", "Ehrlich et al.", "Frizzell et al.", "Molina-Molina et al.", "Minervini et al.", "Malekinejat et al", "Le Guevel and Pakdel"), 
  EC50 = c(0.03, 0.01, 0.02, 0.06, 1.4e-03, 0.05, 0.01),
  Chronic = c(77, 24, 55, 143, 3.36, 120, 24), #{ug/kg bw}
  Accute =  c(114, 38, 74, 228, 5.32, 190, 38), #{ug/kg bw}
  Mycotoxin = c("alpha-Zearalenol")
)


together = rbind(ZENec50_nM.df,aZELec50_nM.df) 
ec50_nM.df <- data.frame(
  Study = c("Diana et al.", "Ehrlich et al.", "Frizzell et al.", "Molina-Molina et al.", "Minervini et al.", "Malekinejat et al", "Le Guevel and Pakdel"),
  ZEN = c(1.58, 0.49, 1.60, 3.81, 0.31, 1.64, 0.06),
  a_ZEL = c(0.03, 0.01, 0.02, 0.06, 1.4e-03, 0.05, 0.01))%>% pivot_longer(ZEN:a_ZEL, names_to = "Compound", values_to = "EC50")

predicted <- data.frame(
  Study = c("Mechanistic-PBK", "Diana", "Mechanistic-PBK", "Mechanistic-PBK", "Diana", "Diana"), 
  HBGV = c("TDI", "TDI", "LB", "UB", "LB", "UB"),
  ZEN = c(C_ZEN.TDI, 0.0076, C_ZEN.LB, C_ZEN.UB, 0.00007284591461306996, 0.008893630661578934),
  a_ZEL = c(C_aZEL.TDI, 0.00029, C_aZEL.LB, C_aZEL.UB, 0.00000274654402958872, 0.000335320770213024)
) %>% mutate(differ =)#pivot_longer(ZEN:a_ZEL, names_to = "Compound", values_to = "EC50")

ggplot(ec50_nM.df, aes(x = Study, y = log10(EC50))) +
  geom_boxplot() +

ggplot()+
  geom_boxplot(data = predicted, aes(x = Study, y = log10(EC50), color = HBGV)) +
  facet_grid(Compound ~ .) +
  coord_flip() +
  labs(x = "Study",
       y = "Concentration (nM)") +
  theme_bw() 

ggplot(predicted, aes(x = Study, y = log10(EC50))) +
  geom_boxplot() +
  # geom_boxplot(data = predicted, aes(x = Study, y = log10(EC50), fill = HBGV)) +
  facet_grid(Compound ~ .) +
  coord_flip() +
  labs(x = "Study",
       y = "Concentration (nM)") +
  theme_bw() 

+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  facet_grid(HBGV ~ .) + #, scales = "free_x", space = "free_x"
  # geom_point(data = predicted, aes(x = Compound, y = log10(EC50), colour = Model, shape = HBGV))+
  

ZENec50_long <- ZENec50_nM.df %>% pivot_longer(Chronic:Accute, names_to = "Exposure_duration", values_to = "Exposure")

ggplot(ec50_nM.df, aes(x = EC50, y = Exposure, colour = Exposure_duration)) + 
  geom_point()

ggplot(ZENec50_long, aes(x = Study, y = EC50, fill = Exposure_duration)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  labs(title = "EC50 Values by Study and Type",
       x = "Study",
       y = "EC50 (nM)",
       fill = "Exposure duration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ZENec50_long <- data.frame(
  EC50 = c(ZENec50_nM.df$EC50, ZENec50_nM.df$EC50),
  Duration = rep(c("Chronic", "Acute"), each = nrow(ZENec50_nM.df)),
  Exposure = c(ZENec50_nM.df$Chronic, ZENec50_nM.df$Accute),
  Mycotoxin = "Zearalenone"
)

# ggplot(ZENec50_long, aes(x = Duration, y = Exposure, fill = Duration)) +
#   geom_boxplot() +
#   geom_abline(y = NOEL) +
#   labs(title = "EC50 Values obtained from in vitro estrogenicity tests",
#        x = "",
#        y = "Exposure (μg/kg bw)") +
#   CP_theme +
#   theme(legend.position = "bottom",  
#         legend.title = element_text(size = 35),
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0))
# # print(Plot.C.Pre)
# ggsave(
#   filename = here(OUTPUT,"EC50.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

ggplot(ZENec50_long, aes(x = Duration, y = Exposure, fill = Duration)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  geom_hline(yintercept = NOEL, linetype = "dashed", color = "seagreen", linewidth = 0.5) +
  geom_text(aes(x = 0.5, y = NOEL, label = "NOEL: 10 μg/kg bw/day"),
            hjust = -0.1, vjust = 0, color = "seagreen", size = 15) +
  labs(title = "EC50 Values obtained from in vitro estrogenicity tests",
       x = "",
       y = "Exposure (μg/kg bw)") +
  coord_flip() +
  CP_theme +
  theme(legend.position = "bottom",  
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))
# print(Plot.C.Pre)
ggsave(
  filename = here(OUTPUT,"EC50.2.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)

NOAEL <-  0.1 #mg/kg b.w. per day.
NOEL <- 10 ##μg/kg b.w.

aZELec50_long <- aZELec50_nM.df %>% pivot_longer(Chronic:Accute, names_to = "Exposure_duration", values_to = "Exposure")

aZENec50_long <- data.frame(
  EC50 = c(aZELec50_nM.df$EC50, aZELec50_nM.df$EC50),
  Duration = rep(c("Chronic", "Acute"), each = nrow(aZELec50_nM.df)),
  Exposure = c(aZELec50_nM.df$Chronic, aZELec50_nM.df$Accute),
  Mycotoxin = "alpha_Zearalenol"
)

long_together <- rbind(ZENec50_long, aZENec50_long)

ggplot(long_together, aes(x = Duration, y = Exposure, fill = Duration)) +
  geom_boxplot() +
  facet_grid(Mycotoxin ~ .) + #, scales = "free_x", space = "free_x"
  scale_fill_brewer(palette = "Set1") +
  geom_hline(yintercept = NOEL, linetype = "dashed", color = "seagreen", linewidth = 0.5) +
  # geom_text(aes(x = 0.5, y = NOEL, label = "NOEL: 10 μg/kg bw/day"),
  #           hjust = -0.1, vjust = 0, color = "seagreen", size = 15) +
  labs(title = "EC50 Values obtained from in vitro estrogenicity tests",
       x = "",
       y = "Exposure (μg/kg bw)") +
  coord_flip() +
  CP_theme +
  theme(legend.position = "bottom",  
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 30, angle = 90, hjust = 0, face = "bold")
        )
# print(Plot.C.Pre)
ggsave(
  filename = here(OUTPUT,"EC50.ZENandaZEL.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)

