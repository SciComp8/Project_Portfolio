library(plotly)
library(haven)

# Data source: 
# This dataset originates from the Duke University Cardiovascular Disease Databank, encompassing information from 3,504 patients and featuring six key variables. These patients sought medical care at Duke University Medical Center due to chest pain concerns. 
acath_data <- read_dta("acath.dta")
acath_data <- acath_data |>
  mutate(sigdz_2 = recode(sigdz, "1" = "Yes", "0" = "No"))


g <- ggplot(
  data=acath_data,
  mapping=aes(cad_dur, color=sigdz_2)) +
  geom_density(alpha=0.5, position="identity") +
  labs(title = "Density of duration of symptoms of coronary artery disease by significant coronary disease",
       x = "Duration of symptoms of coronary artery disease",
       y = "Density")+
  guides(color=guide_legend("Significant coronary disease"))

ggplotly(g)
