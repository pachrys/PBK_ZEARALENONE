library(tidyverse)
library(showtext)
library(here)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()


CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    plot.title.position = "plot", 
    panel.grid = element_line(linewidth = 0.1),  
    axis.title = element_text(size = 34, face = "bold"),
    axis.text = element_text(size = 34),
    axis.line = element_line(linewidth = 0.5),
    legend.text = element_text(size = 34),
    legend.title = element_text(size = 34), 
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.position = "right"
  )


starttime <- 0
stoptime  <- 24                 # hours 
dtout     <- 1/60               # output resolution (min)
Times <- seq(from = starttime, to = stoptime, by = dtout)


# Bile ####
# Ref: Vonk et al 1978, DOI:10.1042/cs0550253
DataBileFlow <- data.frame(
  Time    = c(16.3,
              16.7,
              17.1,
              17.9,
              18.2,
              18.7,
              19.2,
              19.7,
              20.1,
              20.7,
              21.2,
              21.8,
              22.3,
              22.7,
              23.3,
              23.8,
              0.2,
              0.7,
              1.2,
              1.6,
              2.3,
              2.8,
              3.4,
              3.9,
              4.5,
              4.8,
              5.4,
              5.9,
              6.3,
              6.7,
              7.4,
              7.7,
              8.2,
              8.7,
              9.3,
              9.9,
              10.3,
              10.8,
              10.8,
              11.2,
              11.7,
              12.2,
              12.0,
              12.7,
              13.2,
              13.7,
              14.3,
              14.7,
              15.1,
              15.7
              
  ),
  Bileflow = c(13.2568,
              13.3486,
              13.3048,
              13.3288,
              13.5205,
              13.6884,
              13.3486,
              13.5925,
              13.1610,
              13.3486,
              13.6634,
              13.9041,
              13.9281,
              13.7363,
              12.2740,
              11.5308,
              11.2192,
              10.9795,
              10.7877,
              10.5240,
              9.8767,
              9.6849,
              9.4932,
              9.2534,
              9.2534,
              8.9897,
              8.8219,
              8.7740,
              8.6541,
              8.6541,
              8.8699,
              8.6781,
              9.1336,
              9.3452,
              9.5890,
              9.9486,
              10.4281,
              10.9555,
              10.9555,
              11.7466,
              12.2021,
              12.2021,
              12.0822,
              12.1781,
              12.2979,
              12.6096,
              12.8253,
              13.0890,
              13.0890,
              13.2568
  )
)
plot(DataBileFlow)

# Ref: DOI: Drummond et al. 1975 10.1152/ajplegacy.1975.229.5.1427
DataExcretionR <- data.frame(
  Time    = c(16,
              20,
              0,
              4,
              8,
              12,
              16),
  ExcretionR = c(100/100,
               85/100,
               70/100,
               60/100,
               70/100,
               90/100,
               100/100
  )
)
PlotExcretionR <- plot(DataExcretionR)

BileFlowIPf    <- approxfun(DataBileFlow$Time, DataBileFlow$Bileflow, method = "linear", rule = 2)
BileFlow <- BileFlowIPf(Times)
plot(BileFlowIPf(Times))

ExcretionRIPf    <- approxfun(DataExcretionR$Time, DataExcretionR$ExcretionR, method = "linear", rule = 2)
ExcretionR <- ExcretionRIPf(Times)
plot(ExcretionRIPf(Times))

  
df <- tibble(
  BileFlow = BileFlow,
  ExcretionR = ExcretionR,
  BileContentsExcretionR = BileFlow*ExcretionR*1e-6*60,
  Time = Times) 

df <- df %>% mutate(Diff = BileContentsExcretionR/0.000381)

plot(df$Time, df$Diff)

select <- df %>% filter(Time %in% c(0,2,4,5,7,10,13,16,20,24))

plot(select$Time, select$Diff)


write.csv(df, "Bile.csv")
df <- read.csv("Bile.csv")

area.df = df %>% 
  filter(Time >= 7) %>% 
  filter(Time <=19) %>% 
  mutate(BileContentsExcretionR = max(df$BileContentsExcretionR))

ggplot(df, aes(Time, BileContentsExcretionR))+
  geom_area(data = area.df, orientation = "x", fill = "goldenrod1", alpha = 0.2) +
  geom_vline(xintercept = 7, color = "goldenrod1") + 
  geom_vline(xintercept = 19, color = "goldenrod1") + 
  geom_point() +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  geom_vline(xintercept = 0, color = "purple") + 
  annotate("text", x = 1.8, y = 13.5, label="Exposure", 
           color = "purple",size = 10) +
  CP_theme+
  labs(x = "Time (h)", 
       y = "Bile content excretion rate (μl/min)")
ggsave(
  filename = here("BileContentExcretionRate_Rat.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)


# Lymph ####
# Lymph absorption rate (%dose/h)
# DataAbsorptionR <- data.frame(
#   Time    = c(0.5,
#               1.5,
#               2.4,
#               3.5,
#               4.5,
#               5.5,
#               6.5,
#               7.5,
#               8.5,
#               9.5,
#               24
#   ),
#   AbsorptionR = c(0.0424/100,
#                  0.2541/100,
#                  1.8494/100,
#                  5.0541/100,
#                  2.9365/100,
#                  1.8776/100,
#                  1.0165/100,
#                  0.8612/100,
#                  0.4800/100,
#                  0.6918/100,
#                  0.0424/100
#   )
# )
DataAbsorptionR <- data.frame(
  Time    = c(0,
              1.5,
              2.4,
              3.5,
              4.5,
              5.5,
              6.5,
              7.5,
              8.5,
              9.5,
              24
  ),
  AbsorptionR = c(0.0424/9.2,
                  0.2541/9.2,
                  1.8494/9.2,
                  5.0541/9.2,
                  2.9365/9.2,
                  1.8776/9.2,
                  1.0165/9.2,
                  0.8612/9.2,
                  0.4800/9.2,
                  0.6918/9.2,
                  0.0424/9.2
  )
) %>% 
  mutate(Diff = AbsorptionR/0.549358696)

AbsorptionRIPf    <- approxfun(DataAbsorptionR$Time, DataAbsorptionR$AbsorptionR, method = "linear", rule = 2)
AbsorptionR <- AbsorptionRIPf(df$Time)
plot(AbsorptionRIPf(df$Time))

df <- df %>% mutate(LymphR = AbsorptionR)
write.csv(df, here("Bile_LymphN.csv"), row.names = FALSE)

Bile_Lymph <- df
# read_delim(here("Input","Bile_Lymph.csv"),
#                          delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
#   filter(Time <= 24)

area.df <- area.df %>%
  mutate(LymphR = max(Bile_Lymph$LymphR))

ggplot(Bile_Lymph, aes(Time, LymphR))+
  geom_area(data = area.df, orientation = "x", fill = "goldenrod1", alpha = 0.2) +
  geom_vline(xintercept = 7, color = "goldenrod1") +
  geom_vline(xintercept = 19, color = "goldenrod1") +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  geom_vline(xintercept = 0, color = "purple") +
  annotate("text", x = 1.8, y = max(Bile_Lymph$LymphR), label="Exposure",
           color = "purple",size = 10) +
  CP_theme+
  labs(x = "Time (h)",
       y = "Lymphatic absorption rate (/h)")
ggsave(
  filename = here("LymhaticAbsorptionRate_Rat.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)

