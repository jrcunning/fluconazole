# Load libraries
library(tidyverse)
library(lubridate)

# Function to pivot PAM data to long form with column for AOI
ipam_convert <- function(data) {
  data %>% select_if(~ !any(is.na(.))) %>%
    pivot_longer(cols = starts_with("f") | starts_with("y")) %>%
    separate(name, into = c("var", "aoi"), sep = "(?<=[A-Za-z_])(?=[0-9])")
}

# Import PAM data
# List files with exported PAM data
pamfiles <- list.files(path = "data", pattern = "*.csv", full.names = T)

# Import data from each file
pam1 <- pamfiles %>%
  map_dfr(read_delim, delim = ";", .id = "file_id") %>%
  janitor::clean_names() %>%
  mutate(file_id = basename(pamfiles[as.numeric(file_id)]),
         date = as_date(date, format = "%d.%m.%y"))

# For files that have multiple sat pulses -- keep the last one only
pam1 <- pam1 %>%
  group_by(file_id, date) %>%
  filter(no == max(no)) %>%
  ungroup()

# For each source file, convert to long form data with F, FM, and YII for each AOI
pam1 <- pam1 %>%
  nest(-file_id, -date) %>%
  mutate(data2 = map(data, ipam_convert)) %>%
  unnest(data2) %>%
  group_by(file_id, date) %>%
  select(file_id, date, time, aoi, var, value) %>%
  mutate(aoi=as.numeric(aoi))


# Import metadata: order in which racks PAMmed, and order of genets on racks
pammd <- readxl::read_xlsx("Coral project part 2 data file.xlsx") %>%
  mutate(date = as_date(date)) %>%
  drop_na()

pammd<-pammd %>% 
  crossing(aoi=c(1,2,3)) %>%
  mutate(aoi=((position-1)*3)+aoi)
  



# Join PAM data with rack order information (which PAM file corresponds to which rack of corals)
pam <- full_join(pam1, pammd) %>%
  group_by(file_id, date) 

# Change treatment to a categorical variable
pam <- pam %>%
  mutate(treatment = factor(treatment),
         date = factor(date))

anti_join(pam1, pammd)

yii<-pam %>%
  filter(var=="y_ii_")

write_csv(yii,"yii.csv")  

ggplot(yii, aes(x=date,y=value))+
  geom_point()+
  facet_wrap(~species+tank+fragment)+
  theme(strip.text.x = element_text(size = 3))

# Create subset of data for just Bryopsis
bry<-yii%>%
  filter(species=="Bryopsis")

# Plot all Bryopsis fv/fm data
bryfig1<-ggplot(bry,aes(x=date, y=value, group=aoi, color=aoi))+
  geom_point()+
  geom_line()+
  facet_wrap(~tank)
bryfig1

ggsave(filename="Figures/bryfig1.png",plot=bryfig1, width=160, height=80, units="mm")



# Linear model for fv/fm in Bryopsis at end of experiment 
#install.packages("lme4")
library(lme4)
# Fit linear mixed model with treatment and date as fixed predictors, and tank as random
mod<-lmer(value~treatment*date+(1|tank)+(1|aoi),data=bry)

# Calculate the "marginal means" of each treatment at a specific date
library(emmeans)

# Means at the beginning
initial <- emmeans(mod, specs = "treatment", at = list(date = "2023-12-30"))
initial
# Statistically compare these marginal means through pairwise contrasts
contrast(initial, "pairwise")

# Means at the end (Janaury 28th)
final <- emmeans(mod, specs = "treatment", at = list(date = "2024-01-28"))
final
# Statistically compare these marginal means through pairwise contrasts
contrast(final, "pairwise")


# Means at each individual date
all <- emmeans(mod, specs = "treatment", by = "date")
all
# Statistically compare these marginal means through pairwise contrasts
contrast(all, "pairwise")

#### Look at just the tips
# Fit linear mixed model with treatment and date as fixed predictors, and tank as random
tips<-lmer(value~treatment*date+(1|tank), data=filter(bry, aoi == 3))
# Means at each individual date
all <- emmeans(tips, specs = "treatment", by = "date")
all
# Statistically compare these marginal means through pairwise contrasts
contrast(all, "pairwise")

#### Compare 2 control tanks vs. 4 "treated" tanks
bry2 <- bry %>%
  mutate(trt2 = case_when(as.numeric(as.character(treatment)) == 0 ~ "control",
                          as.numeric(as.character(treatment)) > 0  ~ "treated"),
         trt2 = factor(trt2))
# Fit linear mixed model with treatment and date as fixed predictors, and tank as random
mod<-lmer(value~trt2*date+(1|tank),data=filter(bry2, aoi == 3))
# Means at each individual date
all <- emmeans(mod, specs = "trt2", by = "date")
all
# Statistically compare these marginal means through pairwise contrasts
contrast(all, "pairwise")
# sthere are differences between control and 4xtreated tanks, when looking at Bryopsis tips only
