library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(nlme)
library(lme4)
library(ggpubr)
library(multcomp)
library(lsmeans)

#----------------------Analysis of A3 and A4 nicotine resistance-----------------------------:

response_data <- read.csv('nicotine_media_fly_survival.csv', stringsAsFactors = TRUE)

#formatting data:
str(response_data)
response_data$dose <- as.numeric(response_data$dose)
response_data$strain <- as.factor(response_data$strain)
response_data$replicate <- as.factor(response_data$replicate)

#normalize proportion survival by max survival on control food for each strain
#max A4 control percent survival: 55
#max A3 control percent survival: 70

response_data$normalized_survival[response_data$strain == "a4"] <- (response_data$percent[response_data$strain == "a4"]/55)
response_data$normalized_survival[response_data$strain == "a3"] <- (response_data$percent[response_data$strain == "a3"]/70)
response_data$norm_alive <- response_data$normalized_survival*20
response_data$norm_dead <- 20-response_data$norm_alive

#calculate mean survival and SE
data_summary <- response_data %>%
  group_by(strain, dose) %>%
  summarise(avg_survival = mean(percent), avg_normalized = mean(normalized_survival), n = n(),
            sd_raw = sd(percent), SE_raw = sd_raw/(sqrt(n)), 
            sd_norm = sd(normalized_survival), SE_norm = sd_norm/(sqrt(n)))

#generalized linear mixed model using w/normalized alive/dead values (normalized_survival*20 = alive, 20-alive = dead)

a4_AD_norm <- glmer(cbind(norm_alive, norm_dead) ~ dose + (1|replicate), family = binomial, data = subset(response_data, strain == 'a4'))
summary(a4_AD_norm)

a3_AD_norm <- glmer(cbind(norm_alive, norm_dead) ~ dose + (1|replicate), family = binomial, data = subset(response_data, strain == 'a3'))
summary(a3_AD_norm)

#generalized linear hypothesis testing on normalized (alive, dead) values:
dose_asfactor <- response_data
dose_asfactor$dose <- as.factor(dose_asfactor$dose)
str(dose_asfactor)

a4_AD_norm_factor <- glmer(cbind(norm_alive, norm_dead) ~ dose + (1|replicate), family = binomial, data = subset(dose_asfactor, strain == 'a4'))
summary(a4_AD_norm_factor)
summary(glht(a4_AD_norm_factor, mcp(dose = "Tukey")))

a3_AD_norm_factor <- glmer(cbind(norm_alive, norm_dead) ~ dose + (1|replicate), family = binomial, data = subset(dose_asfactor, strain == 'a3'))
summary(a3_AD_norm_factor)
summary(glht(a3_AD_norm_factor, mcp(dose = "Tukey")))

#define function to calculate lethal median concentration (LC50):
dose.p.glmm <-  function(obj, cf = 1:2, p = 0.5) {
  f <- family(obj)
  eta <- f$linkfun(p)
  b <- fixef(obj)[cf]
  x.p <- (eta - b[1L])/b[2L]
  names(x.p) <- paste("p = ", format(p), ":", sep = "")
  pd <- -cbind(1, x.p)/b[2L]
  SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
  res <- structure(x.p, SE = matrix(SE), p = p)
  class(res) <- "glm.dose"
  res
}

#lc50 from AD model
dose.p.glmm(a4_AD_norm, cf=1:2, p=0.5)
dose.p.glmm(a3_AD_norm, cf=1:2, p=0.5)

#lc90 from AD model
dose.p.glmm(a4_AD_norm, cf=1:2, p=0.1)
dose.p.glmm(a3_AD_norm, cf=1:2, p=0.1)

#lc95 from AD model
dose.p.glmm(a4_AD_norm, cf=1:2, p=0.05)
dose.p.glmm(a3_AD_norm, cf=1:2, p=0.05)

#generating fig 1

A1 <- ggplot(data = data_summary, aes(x = dose, y = avg_normalized, color = strain)) +
  geom_point(position=position_dodge(width=0.1), alpha = 0.8, size = 3) +
  geom_errorbar(position=position_dodge(width=0.1), aes(ymin = avg_normalized - SE_norm, ymax = avg_normalized + SE_norm), alpha = .8) +
  geom_smooth(data = response_data, aes(x = dose, y = normalized_survival, color = strain, alpha = .8),
              method = "glm", method.args = list(family = binomial(link = "probit")), se = FALSE) +
  geom_vline(xintercept = 1.05, colour = "#93E1D8", linetype = "dashed", size = 1, alpha = .6) +
  geom_vline(xintercept = 1.90, colour = "#AA7DCE", linetype = "dashed", size = 1, alpha = .6) +
  theme_cowplot()+scale_alpha(guide = 'none')+theme(legend.position = "top")+
  coord_cartesian(xlim = c(0, 5.1), ylim = c(0, 1.0)) +
  scale_x_continuous(breaks = seq(0, 5.0, by = .5)) +
  scale_y_continuous(breaks = seq(0, 1.0, by = .1)) +
  xlab("Nicotine Dose (mM)") +
  ylab("Proportion Survival to Adulthood") +
  scale_color_manual(values = c("#93E1D8", "#AA7DCE", "#FF8C00"),name="Fly Strain",labels=c("A3","A4"))

A1
#----------------------Analysis of effect of nicotine on D. melanogaster growth-----------------:

weight_data <- read.csv('adult_body_mass.csv', stringsAsFactors = TRUE)
shapiro.test(weight_data$weight_mg)

counts <- weight_data %>% group_by(strain, treatment, sex) %>% summarise(n = n())

avg_data <- weight_data %>%
  group_by(strain, treatment, ï..vial, sex) %>%
  summarise(vial_average = mean(weight_mg), num_flies = n())

#generalized linear model of a4 fly mass as a function of fly sex and treatment
a4_mass_model <- glm(vial_average ~ sex*treatment, data = subset(avg_data, strain == "a4"), family = Gamma)
summary(a4_mass_model)
lsmeans(a4_mass_model, list(pairwise ~ sex|treatment, pairwise ~ treatment|sex))

#generalized linear model of a3 fly mass as a function of fly sex and treatment
a3_mass_model <- glm(vial_average ~ sex*treatment, data = subset(avg_data, strain == "a3"), family = Gamma)
summary(a3_mass_model)
lsmeans(a3_mass_model, list(pairwise ~ sex|treatment, pairwise ~ treatment|sex))

#generalized linear model of fly mass (control vials only) as a function of fly strain and treatment
control_mass_model <- glm(vial_average ~ sex*strain, data = subset(avg_data, treatment == "control"), family = Gamma)
summary(control_mass_model)
lsmeans(control_mass_model, list(pairwise ~ sex|strain, pairwise ~ strain|sex))

#generalized linear model of fly mass (nicotine vials only) as a function of fly strain and treatment
nicotine_mass_model <- glm(vial_average ~ sex*strain, data = subset(avg_data, treatment == "nicotine"), family = Gamma)
summary(nicotine_mass_model)
lsmeans(nicotine_mass_model, list(pairwise ~ sex|strain, pairwise ~ strain|sex))

plotdata_growth <- avg_data %>%
  group_by(strain, treatment, sex) %>%
  dplyr::summarize(n = n(),
                   sur_mean = mean(vial_average),
                   sur_sd = sd(vial_average),
                   sur_se = sur_sd / sqrt(n),
                   ci = qt(0.975, df = n - 1) * sur_sd / sqrt(n)) %>% data.frame()

plotdata_growth$treatment <- recode(plotdata_growth$treatment, "control" = "0.00", "nicotine" = "1.25")

#----------------------Analysis of effect of nicotine on parasitism survival-----------------:

parasitism_data <- read.csv('parasitism_survival.csv')
parasitism_data$dose <- as.factor(parasitism_data$dose)
parasitism_data$exp_run <- as.factor(parasitism_data$exp_run)
colnames(parasitism_data)[1] <- "strain"
str(parasitism_data)

parasitism_summary <- parasitism_data %>%
  group_by(strain, parasitism_group, dose) %>%
  summarise(fly_avg_survival = mean(raw_survival), fly_avg_normalized = mean(norm_survival), n = n(),
            wasp_avg_survival = mean(wasp_percent), sd_wasp = sd(wasp_percent), SE_wasp = sd_wasp/sqrt(n),
            sd_raw = sd(raw_survival), SE_raw = sd_raw/(sqrt(n)), 
            sd_norm = sd(norm_survival), SE_norm = sd_norm/(sqrt(n)))

#analysis of a3 parasitism survival:

#model of a3 survival as a function of dose, parasitism, and dose*parasitism
a3_model <- glm(cbind(alive, dead) ~ dose*parasitism_group, family = binomial, data = subset(parasitism_data, strain == 'a3'))
summary(a3_model)

#testing for survival differences in  model of a3 survival as a function of dose, parasitism,
#and dose*parasitism using least-squared means test

lsmeans(a3_model, list(pairwise ~ dose|parasitism_group, pairwise ~ parasitism_group|dose))

#analysis of a4 fly survival:

#model of a4 survival as a function of dose, parasitism, and dose*parasitism
a4_model <- glm(cbind(alive, dead) ~ dose*parasitism_group, family = binomial, data = subset(parasitism_data, strain == 'a4'))
summary(a4_model)

#testing for survival differences in  model of a4 survival as a function of dose, parasitism,
#and dose*parasitism using least-squared means test
lsmeans(a4_model, list(pairwise ~ dose|parasitism_group, pairwise ~ parasitism_group|dose))

#analysis of wasp survival:

#modelling wasp survival as a function of dose and fly strain
wasp_model <- glm(cbind(wasp_alive, wasp_dead) ~ dose*strain, family = binomial, data = parasitism_data)
summary(wasp_model)

#least-squared means test of wasp survival
lsmeans(wasp_model, list(pairwise ~ dose|strain, pairwise ~ strain|dose))

#calculating mean and SE of wasp and fly survival for figure 2 panels A and B
plotdataf <- parasitism_data  %>% 
  group_by(strain,dose,parasitism_group) %>%
  dplyr::summarize(n = n(),
                   sur_mean = mean(raw_survival),
                   sur_sd = sd(raw_survival),
                   sur_se = sur_sd / sqrt(n),
                   ci = qt(0.975, df = n - 1) * sur_sd / sqrt(n)) %>% data.frame()

plotdataw <- parasitism_data  %>% subset(parasitism_group=="w") %>%
  group_by(strain,dose) %>%
  dplyr::summarize(n = n(),
                   sur_mean = mean(wasp_percent),
                   sur_sd = sd(wasp_percent),
                   sur_se = sur_sd / sqrt(n),
                   ci = qt(0.975, df = n - 1) * sur_sd / sqrt(n)) %>% data.frame()


#plots for fig 2 panels A, B, and C
A2 <- ggplot(plotdataf, aes(x = dose, y = sur_mean, group = interaction(strain, parasitism_group), color = strain,linetype=parasitism_group)) +
  geom_point(size = 3,show.legend = F) + geom_line(show.legend = F) +ylab("Percent Fly Survival") +ylim(-1,100)+
  geom_errorbar(aes(ymin = sur_mean - sur_sd,ymax = sur_mean + sur_sd),linetype="solid", width = .1)+
  theme_cowplot(font_size = 12)+theme(legend.position="none")+
  scale_linetype_manual(values = c("dashed", "solid"),name="") +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())+
  scale_color_manual(values = c("#93E1D8", "#AA7DCE"))


B2 <- ggplot(plotdataw, aes(x = dose, y = sur_mean, group = strain, color = strain)) +
  geom_point(size = 3,show.legend = F) + geom_line(show.legend = F) +ylab("Percent Wasp Survival") + ylim(-1,100)+
  geom_errorbar(aes(ymin = sur_mean - sur_sd,ymax = sur_mean + sur_sd), width = .1,show.legend = F)+
  theme_cowplot(font_size = 12)+theme(legend.position="none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = NULL)+
  scale_color_manual(values = c("#93E1D8", "#AA7DCE"))

C2 <- ggplot(plotdata_growth, aes(x = treatment, y = sur_mean, group = strain:sex, color = strain, shape = sex)) +
  geom_point(size = 3,show.legend = F) + geom_line(show.legend = F) +ylab("Adult Fly Mass (mg)") +ylim(.1,.3)+
  geom_errorbar(aes(ymin = sur_mean - sur_sd,ymax = sur_mean + sur_sd),linetype="solid", width = .1)+
  theme_cowplot(font_size = 12)+theme(legend.position="none")+
  scale_linetype_manual(values = c("solid"),name="") +
  xlab("Nicotine Concentration (mM)") +
  theme()+
  scale_color_manual(values = c("#93E1D8", "#AA7DCE"))

#----------------------Analysis of nicotine accumulation in A3/A4----------------------------:

#formatting data:
total_data <- read.csv("GC-MS_Nicotine_Cotinine_flies.csv")
total_data <- filter(total_data, ng.fly < 150)

total_data$stage <- factor(total_data$stage, levels = c("L3", "P1", "P3", "A1", "A3"))
total_data$stage <- recode(total_data$stage, "L3" = "3rd Instar", "P1" = "Day 1 Pupa", "P3" = "Day 3 Pupa", "A1" = "Day 1 Adult", "A3" = "Day 3 Adult")
total_data$Cotinine.Nicotine.Area.Ratio <- as.numeric(as.character(total_data$Cotinine.Nicotine.Area.Ratio))


#calculating mean nicotine accumulation and SE across developmental stages
mean_bystage <- total_data %>%
  subset(treatment == "nicotine") %>%
  group_by(strain, treatment, stage) %>%
  summarise(vial.count = n(), stage.mean = mean(ng.fly), stdev = sd(ng.fly), SE = stdev/(sqrt(vial.count))) 

total_data<-subset(total_data,treatment=="nicotine")
total_data$Cotinine.Area <- as.numeric(total_data$Cotinine.Area)

#a4 nicotine pairwise tests
wilcox.test(total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "3rd Instar"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 3 Adult"])

#a3 nicotine pairwise tests
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "3rd Instar"], total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"], total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"], total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 3 Adult"])

#a3 vs a4 nicotine tests (ng)
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "3rd Instar"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "3rd Instar"])
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$ng.fly[total_data$strain == "a3" & total_data$stage == "Day 3 Adult"], total_data$ng.fly[total_data$strain == "a4" & total_data$stage == "Day 3 Adult"])

#a3 vs a4 nicotine tests (GCMS area)
wilcox.test(total_data$Nicotine.Area[total_data$strain == "a3" & total_data$stage == "3rd Instar"], total_data$Nicotine.Area[total_data$strain == "a4" & total_data$stage == "3rd Instar"])
wilcox.test(total_data$Nicotine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"], total_data$Nicotine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$Nicotine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"], total_data$Nicotine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$Nicotine.Area[total_data$strain == "a3" & total_data$stage == "Day 3 Adult"], total_data$Nicotine.Area[total_data$strain == "a4" & total_data$stage == "Day 3 Adult"])


#a4 cotinine pairwise tests
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "3rd Instar"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 3 Adult"])

#a3 cotinine pairwise tests
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "3rd Instar"], total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"], total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"], total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 3 Adult"])

#a3 vs a4 cotinine tests
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "3rd Instar"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "3rd Instar"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Pupa"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Pupa"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 1 Adult"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 1 Adult"])
wilcox.test(total_data$Cotinine.Area[total_data$strain == "a3" & total_data$stage == "Day 3 Adult"], total_data$Cotinine.Area[total_data$strain == "a4" & total_data$stage == "Day 3 Adult"])


#plots for fig 2 panels C and D

plotdata_short <- total_data %>% subset(treatment=="nicotine"&ng.fly<160) %>%
  group_by(strain,stage) %>%
  dplyr::summarize(n = n(),
                   nic_mean = mean(ng.fly),
                   nic_sd = sd(ng.fly),
                   nic_se = nic_sd / sqrt(n),
                   ci = qt(0.975, df = n - 1) * nic_sd / sqrt(n),
                   cot_mean = mean(Cotinine.Area),
                   cot_sd = sd(Cotinine.Area),
                   cot_se = nic_sd / sqrt(n),
                   cot_nic_ratio = cot_mean/nic_mean,
                   ci = qt(0.975, df = n - 1) * cot_sd / sqrt(n))
                   

long_format <- total_data %>%
  pivot_longer(c("Nicotine.Area", "Cotinine.Area"), names_to = 'nicotine.cotinine', values_to = "GCMS_area")

long_format$GCMS_area <- long_format$GCMS_area/1000

plotdata_long <- long_format %>%
  group_by(strain, stage, nicotine.cotinine) %>%
  dplyr::summarize(n = n(),
                   mean_area = mean(GCMS_area),
                   area_sd = sd(GCMS_area),
                   area_se = area_sd / sqrt(n),
                   ci = qt(0.975, df = n - 1) * area_sd / sqrt(n)) 

#generating panels 2D and 2E
D2 <- ggplot(plotdata_short, aes(x = stage, y = nic_mean, group = interaction(strain, stage), color = strain)) +
  geom_point(position=position_dodge(width=0.4), alpha = 0.8, size = 2.5,show.legend = F, shape = 15) +
  ylab("Ng Nicotine per Fly") +
  ylim(0, 70)+
  geom_errorbar(aes(ymin = nic_mean - nic_sd ,ymax = nic_mean + nic_sd),position=position_dodge(width=0.4),
                linetype="solid", width = .25)+
  theme_cowplot(font_size = 12)+theme(legend.position="none")+
  scale_linetype_manual(values = c("solid"),name="") +
  xlab("Developmental Stage") +
  theme()+
  scale_color_manual(values = c("#93E1D8", "#AA7DCE"))

E2 <- ggplot(plotdata_long, aes(x = stage, y = mean_area, group = interaction(strain, stage, nicotine.cotinine), color = strain, shape = nicotine.cotinine)) +
  geom_point(position=position_dodge(width=0.4), size = 2.5, alpha = .8, show.legend = F) +
  scale_shape_manual(values = c(0, 15)) +
  ylab("GCMS Area") +
  scale_y_continuous(breaks=seq(0, 275, 25)) +
  geom_errorbar(aes(ymin = mean_area - area_sd,ymax = mean_area + area_sd),position=position_dodge(width=0.4),
                linetype="solid", width = .5)+
  theme_cowplot(font_size = 12)+theme(legend.position="none")+
  scale_linetype_manual(values = c("solid"),name="") +
  xlab("Developmental Stage") +
  theme()+
  scale_color_manual(values = c("#93E1D8", "#AA7DCE"))

#generating fig 2
ABC2 <- plot_grid(A2, B2, C2, nrow = 3,  labels = c("A", "B", "C"), align = "hv", hjust=-4.5,vjust=c(4.5,2), axis = "b")
DE2 <- plot_grid(D2, E2, NULL, nrow = 3, labels = c("D", "E", ""), align = "hv", hjust=-4.5,vjust=c(4.5,2), axis = "b")

plot_grid(ABC2, DE2)
