#*****************************************************************************#
#------------------------------Ploting Function-------------------------------#

CP_theme <- function() {
  theme_bw()+
    theme(
      text = element_text(size = 14), 
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.line = element_line(linewidth = 0.05),
      panel.border = element_blank(), 
      panel.background = element_blank(),
      panel.grid = element_line(linewidth = 0.1)
  )
}