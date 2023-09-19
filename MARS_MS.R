#Code MARS-5 MS prospective study

#libraries
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(haven) 
library(foreign)
library(tidyverse)
library(AdhereR)
library(survival) 
library(survminer) 
library(excel.link) 
library(rms)
library(modEvA)
library(stargazer)
library(riskRegression)
library(CoxR2)
library(RColorBrewer)
library(pwr)
library(ggstatsplot)
library(readxl)
library(psych)
library(BlandAltmanLeh)
library(pROC)
library(caret)
library(caTools)
library(Information)
library(ISLR)
library(e1071)
library(randomForest)

#database
dfPts_mars <- xl.read.file( , password =)
dfRx_mars <- xl.read.file(, password =)

#DMF dispense dates
un_dfRx_mars <- distinct(dfRx_mars, ID)
colnames(dfRx_mars)[which(names(dfRx_mars) == 'issuedate')] <- 'DATE'
dfRx_mars <- dfRx_mars[order(dfRx_mars$IDcode, dfRx_mars$DATE),]
# df1st <- dfRx %>% group_by(IDcode) %>% filter(DATE == min(DATE)) druga opcija 
df1st <- dfRx_mars %>% group_by(IDcode) %>%
  mutate(rank = rank(DATE, ties.method= "first"))
df1st <- df1st %>% group_by(IDcode) %>% filter(rank == 1)
colnames(df1st)[which(names(df1st) == 'DATE')] <- 'DATE_1st'

#link database and dispenses
df1stF <- df1st %>% dplyr::select(DATE_1st,ID)
df1stF <- df1stF[order(df1stF$ID),]
dfPts_mars <- dfPts_mars[order(dfPts_mars$IDcode),]
dfPts_mars2 <- dfPts_mars %>% left_join(df1stF,by="ID")
colnames(dfPts_mars2)[which(names(dfPts_mars2) == 'ID')] <- 'PATIENT_ID'

dfPts_mars2<- dfPts_mars2[order(dfPts_mars$IDcode),]
dfRx_mars <- dfRx_mars[order(dfRx_mars$IDcode, dfRx_mars$DATE),]
df= dfRx_mars %>% left_join(dfPts_mars2,by="IDcode")
un_df <- distinct(df, IDcode) #patient, unique.

dbAdhere_mars <- df
un_dbAdhere_mars <- distinct(dbAdhere_mars, PATIENT_ID)
dbAdhere_mars$DATE <- as.Date(dbAdhere_mars$DATE)

dbAdhere_mars <- dbAdhere_mars[order(dbAdhere_mars$PATIENT_ID, dbAdhere_mars$DATE),]

#follow-up brain MRI
dbAdhere_mars$DATE_1st <- as.Date(dbAdhere_mars$DATE_1st)
dbAdhere_mars$d.MR.1stControl.INFO <- as.Date(dbAdhere_mars$d.MR.1stControl.INFO)
dbAdhere_mars$d.MR.2ndControl.INFO <- as.Date(dbAdhere_mars$d.MR.2ndControl.INFO)
dbAdhere_mars$follow_up_MR <- dbAdhere_mars$d.MR.2ndControl - dbAdhere_mars$DATE_1st + 1
dbAdhere_mars$follow_up_MR <- as.numeric(dbAdhere_mars$follow_up_MR)

#observation window brain MRI
dbAdhere_mars$MR.diff.t <- dbAdhere_mars$d.MR.2ndControl - dbAdhere_mars$d.MR.1stControl
dbAdhere_mars$MR.diff.t <- as.numeric(dbAdhere_mars$MR.diff.t)
difftime(dbAdhere_mars$d.MR.2ndControl, dbAdhere_mars$d.MR.1stControl, units = "weeks")
mean(dbAdhere_mars$MR.diff.t/30)
sd(dbAdhere_mars$MR.diff.t/30)

#follow-up adherence
dbAdhere_mars$follow_up_time_adh <- dbAdhere_mars$d.MR.2ndControl - dbAdhere_mars$DATE_1st + 1
dbAdhere_mars$follow_up_time_adh <- as.numeric(dbAdhere_mars$follow_up_time_adh)

dbAdhere_mars$MRstart <- dbAdhere_mars$d.MR.1stControl - dbAdhere_mars$DATE_1st
dbAdhere_mars$MRstart <- as.numeric(dbAdhere_mars$MRstart)

#CMA6
df_CMA6_episode <- CMA_per_episode(CMA="CMA6",
                                   data=dbAdhere_mars,
                                   ID.colname = "IDcode",
                                   event.date.colname = "DATE",
                                   event.duration.colname = "DURATION",
                                   event.daily.dose.colname = "PERDAY",
                                   medication.class.colname = "CATEGORY",
                                   carryover.within.obs.window = TRUE,
                                   carry.only.for.same.medication = TRUE,
                                   consider.dosage.change = TRUE,
                                   medication.change.means.new.treatment.episode = FALSE,
                                   dosage.change.means.new.treatment.episode = FALSE,
                                   maximum.permissible.gap = 60,
                                   maximum.permissible.gap.unit = "days",
                                   followup.window.start = 0,
                                   followup.window.duration = 365*1.75,
                                   followup.window.duration.unit = "days",
                                   observation.window.start = 0,
                                   observation.window.duration = 365*1.75,
                                   observation.window.duration.unit = "days",
                                   event.interval.colname = "event.interval",
                                   gap.days.colname = "gap.days",
                                   date.format = "%Y-%m-%d",
                                   summary = "Base_CMA6")


CMA6_episode <- getCMA(df_CMA6_episode)
mean(CMA6_episode$CMA)
sd(CMA6_episode$CMA)
CMA6_episode = CMA6_episode %>% group_by(IDcode) %>% mutate(weighted_CMA = weighted.mean(CMA, episode.duration))
mean(CMA6_episode$weighted_CMA)
sd(CMA6_episode$weighted_CMA)
CMA6_episode = CMA6_episode %>% filter(episode.ID == min(episode.ID))
mean(CMA6_episode$weighted_CMA)
sd(CMA6_episode$weighted_CMA)
sum(CMA6_episode$weighted_CMA>=0.9)
CMA6_episode_Pts <- CMA6_episode %>% left_join(dfPts_mars2,by="IDcode")
CMA6_episode_Pts %>% group_by(naive.pts) %>% summarise(mean(weighted_CMA))
sum(CMA6_episode_Pts$weighted_CMA>=0.9)
sum(CMA6_episode_Pts$weighted_CMA>=0.85)

#CMA7 criterion validity
dbAdhere_mars <- dbAdhere_mars[order(dbAdhere_mars$IDcode, dbAdhere_mars$DATE),]
df_cma7_episodeMR_gap60 <- CMA_per_episode(CMA="CMA7",
                                           data=dbAdhere_mars,
                                           ID.colname = "IDcode",
                                           event.date.colname = "DATE",
                                           event.duration.colname = "DURATION",
                                           event.daily.dose.colname = "PERDAY",
                                           medication.class.colname = "CATEGORY",
                                           carryover.within.obs.window = TRUE,
                                           carry.only.for.same.medication = TRUE,
                                           consider.dosage.change = TRUE,
                                           medication.change.means.new.treatment.episode = FALSE,
                                           dosage.change.means.new.treatment.episode = FALSE,
                                           maximum.permissible.gap = 60,
                                           maximum.permissible.gap.unit = "days",
                                           followup.window.start = 0,
                                           followup.window.duration = "follow_up_MR",
                                           followup.window.duration.unit ="days",
                                           observation.window.start = "MRstart",
                                           observation.window.start.unit = "days",
                                           observation.window.duration = "MR.diff.t",
                                           observation.window.duration.unit = "days",
                                           event.interval.colname = "event.interval",
                                           gap.days.colname = "gap.days",
                                           date.format = "%Y-%m-%d",
                                           summary = "Base_CMA7_MR");

CMA7_MR_gap60 <- getCMA(df_cma7_episodeMR_gap60)
mean(CMA7_MR_gap60$CMA)
sd(CMA7_MR_gap60$CMA)
sum(CMA7_MR_gap60$CMA>=0.9)
sum(CMA7_MR_gap60$CMA>=0.85)

CMA7_MR_gap60 = CMA7_MR_gap60 %>% group_by(IDcode) %>% mutate(weighted_CMA = weighted.mean(CMA, episode.duration))
mean(CMA7_MR_gap60$weighted_CMA)
sd(CMA7_MR_gap60$weighted_CMA)
CMA7_MR_gap60 = CMA7_MR_gap60 %>% filter(episode.ID == min(episode.ID))
mean(CMA7_MR_gap60$weighted_CMA)
sd(CMA7_MR_gap60$weighted_CMA)
sum(CMA7_MR_gap60$weighted_CMA>=0.9)
CMA7_MR_gap60 <- CMA7_MR_gap60 %>% left_join(dfPts_mars2,by="IDcode")
CMA7_MR_gap60 %>% group_by(naive.pts) %>% summarise(mean(weighted_CMA))
sum(CMA7_MR_gap60$weighted_CMA>=0.9)
sum(CMA7_MR_gap60$weighted_CMA>=0.85)

#episode duration
CMA7_MR_gap60$episode.duration.REAL <- CMA7_MR_gap60$episode.end - CMA7_MR_gap60$episode.start
CMA7_MR_gap60$episode.duration.REAL <- as.numeric(CMA7_MR_gap60$episode.duration.REAL)

CMA6_episode_Pts <- transform(CMA6_episode_Pts, combo.MR.2ndControl = ifelse(combo.lesion.2ndControl=="DA",1,0))
CMA7_MR_gap60 <- transform(CMA7_MR_gap60, combo.MR.2ndControl = ifelse(combo.lesion.2ndControl=="DA",1,0))

CMA6_episode_Pts$status <- ifelse(CMA6_episode_Pts$weighted_CMA>=0.85,1,0)
dfPts_mars2$status <- as.factor(dfPts_mars2$status)

dfPts_mars$status <- ifelse(CMA6_episode_Pts$CMA>=0.85,1,0)
dfPts_mars2$status <- as.factor(dfPts_mars2$status)

dfPts_mars$status2<- ifelse(CMA6_episode_Pts$CMA>=0.90,1,0)
dfPts_mars2$status2 <- as.factor(dfPts_mars2$status2)

CMA6_episode_Pts$status2<- ifelse(CMA6_episode_Pts$weighted_CMA>=0.90,1,0)
dfPts_mars2$status2 <- as.factor(dfPts_mars2$status2)

#psychometric tests

#test-retest
ICC(dfPts_mars[,c("MARS.1", "MARS.2")], missing = TRUE, alpha = .05, lmer = TRUE, check.keys = FALSE)

#Cronbach alfa
alpha(dfPts_mars_C[,c("M1_2","M2_2","M3_2","M4_2","M5_2")])

#criterion validity brain MRI
CMA7_MR_gap60_Pts_final$status2MR <- ifelse(CMA7_MR_gap60_Pts_final$CMA7.gap60.MR>=0.85,1,0)
CMA7_MR_gap60_Pts_final$status3MR <- ifelse(CMA7_MR_gap60_Pts_final$CMA7.gap60.MR>=0.90,1,0)

CMA7_MR_gap60$status2MR <- ifelse(CMA7_MR_gap60$weighted_CMA>=0.85,1,0)
CMA7_MR_gap60$status3MR <- ifelse(CMA7_MR_gap60$weighted_CMA>=0.90,1,0)

mean(CMA7_MR_gap60$MR.diff.lesion_volume)
sd(CMA7_MR_gap60$MR.diff.lesion_volume)

mean(CMA7_MR_gap60$MR.diff.lesion_count)
sd(CMA7_MR_gap60$MR.diff.lesion_count)

CMA7_MR_gap60$MARS.m <- (CMA7_MR_gap60$MARS.1 + CMA7_MR_gap60$MARS.2)/2
mean(CMA7_MR_gap60$MARS.m)
sd(CMA7_MR_gap60$MARS.m)

CMA7_MR_gap60$BMQ.gen.m <- (CMA7_MR_gap60$BMQ.gen + CMA7_MR_gap60$BMQ2.gen)/2
mean(CMA7_MR_gap60$BMQ.gen.m)
sd(CMA7_MR_gap60$BMQ.gen.m)

CMA7_MR_gap60$BMQ.spec.m <- (CMA7_MR_gap60$BMQ1.spec + CMA7_MR_gap60$BMQ2.spec)/2
mean(CMA7_MR_gap60$BMQ.spec.m)
sd(CMA7_MR_gap60$BMQ.spec.m)

#descriptive statistics
sum(CMA6_episode_Pts$sex=='F')
sum(CMA6_episode_Pts$sex=='M')
mean(CMA6_episode_Pts$dur.disease.y)
sd(CMA6_episode_Pts$dur.disease.y)
mean(CMA6_episode_Pts$EDSS.1stRxTec)
sd(CMA6_episode_Pts$EDSS.1stRxTec)
mean(CMA6_episode_Pts$age)
sd(CMA6_episode_Pts$age)
sum(CMA6_episode_Pts$naive.pts)
mean(CMA7_MR_gap60$MARS.1)
sd(CMA7_MR_gap60$MARS.1)
mean(CMA7_MR_gap60$MARS.2)
sd(CMA7_MR_gap60$MARS.2)
mean(CMA6_episode_Pts$BMQ.gen)
sd(CMA6_episode_Pts$BMQ.gen)
mean(CMA6_episode_Pts$BMQ2.gen)
sd(CMA6_episode_Pts$BMQ2.gen)
mean(CMA6_episode_Pts$BMQ1.spec)
sd(CMA6_episode_Pts$BMQ1.spec)
mean(CMA6_episode_Pts$BMQ2.spec)
sd(CMA6_episode_Pts$BMQ2.spec)

#PDC according to MARS-5 score
CMA6_episode_Pts %>% group_by(MARS.1) %>% summarise(mean(weighted_CMA))
CMA6_episode_Pts %>% group_by(MARS.2) %>% summarise(mean(weighted_CMA))
CMA6_episode_Pts %>% group_by(MARS.2) %>% summarise(sd(weighted_CMA))

CMA7_MR_gap60 %>% group_by(MARS.2) %>% summarise(mean(CMA))
CMA7_MR_gap60 %>% group_by(MARS.2) %>% summarise(mean(weighted_CMA))
CMA7_MR_gap60 %>% group_by(MARS.2) %>% summarise(sd(weighted_CMA))

#scatterplot
ggplot(data = CMA7_MR_gap60,aes(x = MARS.2 , y = weighted_CMA)) + 
  geom_point(size=3, shape=16, color ="black", fill="black")+theme_classic()+scale_x_continuous(name="MARS-5 follow-up brain MRI", breaks = seq(15,25, by=1))+theme(axis.title.x = element_text(size = 15))+theme(axis.title.x = element_text(size = 15))+theme(axis.text = element_text(size = 15))+scale_y_continuous(name="PDC adherence (weighted CMA7)", limits=c(0.3,1))
labs(x = "MARS follow-up brain MRI", y = "PDC adherence (CMA7)",
     title ="Scatterplot of CMA7 and MARS.2", 
     fun='pct',
     subtitle = "Correlation between CMA7 and MARS_2",
     caption = "MJ",
     alt = "XY")

#covariates MARS-5 --------
linear2 <- glm(MARS.2 ~ age + sex + EDSS.1stRxTec + naive.pts + dur.disease.y, data = CMA7_MR_gap60, family = "gaussian")
summary(linear2)

model.final_MARS.2 <- step(linear2)
summary(model.final_MARS.2)

#statistical power
cohen.ES(test = "r", size = "small")
r <- pwr.f2.test(u = 3, v = 40, f2 = 0.20, sig.level = 0.05, power = NULL)
pwr.f2.test(u = 1, v = 40, f2 = 0.71, sig.level = 0.05, power = NULL)
pwr.f2.test(u = 5, v = 40, f2 = 0.20, sig.level = 0.05, power = NULL)

#correlation CMA and MARS-5

#linear regression
linearMR3 <- glm(weighted_CMA ~ MARS.2 + age + sex + EDSS.1stRxTec + naive.pts + dur.disease.y, data = CMA7_MR_gap60, family = "gaussian")
summary(linearMR3)
confint(linearMR3)
model.final_linearMR3 <- step(linearMR3)
summary(model.final_linearMR3)
confint(model.final_linearMR3)
correlation2 <- cor(CMA7_MR_gap60$MARS.2, CMA7_MR_gap60$weighted_CMA, method = 'pearson')
print(correlation2)

#normal distribution test

shapiro.test(CMA6_episode_Pts$BMQ2.gen)
shapiro.test(CMA6_episode_Pts$BMQ2.spec)

shapiro.test(CMA6_episode_Pts$MARS.2)
shapiro.test(CMA6_episode_Pts$MARS.m)

shapiro.test(CMA6_episode_Pts$EDSS.1stRxTec)

#threshold

#cut-off score 24
CMA6_episode_Pts$MARS.1_24 <- ifelse(CMA6_episode_Pts$MARS.1>=24,1,0)
CMA6_episode_Pts$MARS.1_24 <- as.factor(CMA6_episode_Pts$MARS.1_24)
CMA6_episode_Pts$MARS.2_24 <- ifelse(CMA6_episode_Pts$MARS.2>=24,1,0)
CMA6_episode_Pts$MARS.2_24 <- as.factor(CMA6_episode_Pts$MARS.2_24)
CMA6_episode_Pts$MARS.m_24 <- ifelse(CMA6_episode_Pts$MARS.m>=24,1,0)
CMA6_episode_Pts$MARS.m_24 <- as.factor(CMA6_episode_Pts$MARS.m_24)

CMA7_MR_gap60$MARS.1_24MR <- ifelse(CMA7_MR_gap60$MARS.1>=24,1,0)
CMA7_MR_gap60$MARS.1_24MR <- as.factor(CMA7_MR_gap60$MARS.1_24MR)
CMA7_MR_gap60$MARS.2_24MR <- ifelse(CMA7_MR_gap60$MARS.2>=24,1,0)
CMA7_MR_gap60$MARS.2_24MR <- as.factor(CMA7_MR_gap60$MARS.2_24MR)
CMA7_MR_gap60$MARS.m_24MR <- ifelse(CMA7_MR_gap60$MARS.m>=24,1,0)
CMA7_MR_gap60$MARS.m_24MR <- as.factor(CMA7_MR_gap60$MARS.m_24MR)

CMA7_MR_gap60 %>% group_by(MARS.2_24MR) %>% summarise(mean(MARS.2))
CMA7_MR_gap60 %>% group_by(MARS.2_24MR) %>% summarise(sd(MARS.2))

sum(CMA7_MR_gap60$MARS.2_24MR==1)
sum(CMA7_MR_gap60$MARS.2_24MR==0)


#linearna regresija criterion validity BMQ
model_BMQ11 <- glm(BMQ2.gen ~ MARS.2  + age +sex + dur.disease.y + EDSS.1stRxTec + naive.pts, data = CMA7_MR_gap60, family = "gaussian")
summary(model_BMQ11)
model.final_BMQ11 <- step(model_BMQ11)
summary(model.final_BMQ11)
confint(model.final_BMQ11)
correlation <- cor(CMA7_MR_gap60$MARS.2, CMA7_MR_gap60$BMQ2.gen, method = 'pearson')
print(correlation)

model_BMQ10 <- glm(BMQ2.spec ~ MARS.2 + age +sex + dur.disease.y + EDSS.1stRxTec + naive.pts, data = CMA7_MR_gap60, family = "gaussian")
summary(model_BMQ10)
confint(model_BMQ10)
model.final_BMQ10 <- step(model_BMQ10)
summary(model.final_BMQ10)
confint(model.final_BMQ10)

#correlation analysis

#threshold score 24

# Mann-Whitney test quantification
wilcox.test(CMA7_MR_gap60$MR.diff.lesion_volume ~ CMA7_MR_gap60$MARS.2_24MR, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MR.diff.lesion_count ~ CMA7_MR_gap60$MARS.2_24MR, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MR.diff.new_appearing_lesion_count ~ CMA7_MR_gap60$MARS.2_24MR, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MR.diff.new_appearing_lesion_volume ~ CMA7_MR_gap60$MARS.2_24MR, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MR.diff.shrinking_lesion_volume ~ CMA7_MR_gap60$MARS.2_24MR, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MR.diff.shrinking_lesion_count ~ CMA7_MR_gap60$MARS.2_24MR, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MARS.2 ~ CMA7_MR_gap60$MR.diff.enlarging_lesion_volume_cat, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MARS.2 ~ CMA7_MR_gap60$MR.diff.new_appearing_lesion_volume_cat, exact=FALSE, paired=FALSE)
wilcox.test(CMA7_MR_gap60$MARS.2 ~ CMA7_MR_gap60$MR.diff.shrinking_lesion_volume_cat, exact=FALSE, paired=FALSE)

#Mann Whitney test radiologists' examination
wilcox.test(CMA7_MR_gap60$MARS.2 ~ CMA7_MR_gap60$combo.MR.2ndControl, exact=FALSE, paired=FALSE)
wilcox.test(CMA6_episode_Pts$MARS.2 ~ CMA6_episode_Pts$combo.MR.2ndControl, exact=FALSE, paired=FALSE)

MR <- table(CMA7_MR_gap60$MARS.2_23MR,CMA7_MR_gap60$combo.MR.2ndControl)
chisq.test(MR)

#logistic regression MARS.2
model_MR21 <- glm(combo.MR.2ndControl ~ MARS.2 + age +sex + dur.disease.y + EDSS.1stRxTec + naive.pts, data = CMA7_MR_gap60, family = "binomial"(link = "logit"))
summary(model_MR21)
exp(coefficients(model_MR21))
exp(confint(model_MR21))
model.final_MR21 <- step(model_MR21)
summary(model.final_MR21)
exp(coefficients(model.final_MR21))
exp(confint(model.final_MR21))

#logistic regression BMQ
model_BMQ <- glm(combo.MR.2ndControl ~ MARS_2 + age + relevel(dfPts_mars$sex, ref='M') + dur.disease.y + EDSS.1stRxTec + naive.pts + recode.cnt.relapse.Before1stRxTec, data = dfPts_mars, family = "binomial"(link = "logit"))
model.final_BMQ <- step(model_BMQ)
summary(model.final_BMQ)

# ROC curve----------------------------------------------------------------

model_CMA <- glm(status ~ MARS.2, data=dfPts_mars2, family = "binomial")
summary(model_CMA)
prediction <- predict(model_CMA, newdata=dfPts_mars, type="response")
prediction <- ifelse(prediction >0.5,1,0)

dfPts_mars2$predictions <- prediction
prediction <- as.factor(prediction)
dfPts_mars2$status <- as.factor(dfPts_mars2$status)
table1 <- table(ActualValue=dfPts_mars2$status, PredictedValue=prediction)
table1
confusionMatrix(table1)

view(table1)
optimal <- optimalCutoff(dfPts_mars2$status, prediction)
confusionMatrix(dfPts_mars2$status, prediction)

model_CMA2 <- glm(status2 ~ MARS.2, data=dfPts_mars2, family = "binomial")
summary(model_CMA2)
prediction <- predict(model_CMA, newdata=dfPts_mars2, type="response")
prediction <- ifelse(prediction >0.5,1,0)
dfPts_mars2$predictions <- prediction
table2 <- table(ActualValue=dfPts_mars2$status, PredictedValue=prediction)
view(table2)

plot(x=CMA6_episode_Pts$MARS.2, y=CMA6_episode_Pts$status)

glm.fit <- glm(status ~ MARS.2, data=CMA6_episode_Pts, family=binomial)
lines(CMA6_episode_Pts$MARS.2, glm.fit$fitted.values)

roc(CMA6_episode_Pts$status, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", print.thres.adj=c(15,25))
coord_list[[1]] <- coords(roc, x = "all")
roc.info <- roc(dfPts_mars2$status, glm.fit$fitted.values, legacy.axes=TRUE)
str(roc.info)

roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)

tail(roc.df)

roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]

roc(status, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

#accuracy
TN <- table1[1,1]
TP <- table1[2,2]
FN <- table1[2,1]
FP <- table1[1,2]
acc <- (TP+TN)/(TP+TN+FN+FP)
TPR <- TP/(TP+FN) #sensitivity
TNR <- TN/(TN+FP) #specificity
FPR <- FP/(FP+TN) #false positive
FNR <- FN/(FN+TP) #false negative
PPV <- TP/(TP+FP) #positive predictive value
NPV <- TN/(TN+FN) #negative predictive value

ROC_prediction <- predict(model_CMA, newdata=dfPts_mars2, type="response")
par(pty = "s")
#as.numeric(test_roc$auc)
ROC_plot<- roc(dfPts_mars2$status, ROC_prediction, plot=TRUE, print.auc=TRUE, main="ROC curve", col="black")
coords(ROC_plot, "best", ret=c("threshold", "specificity","sensitivity"), best.method=c("closest.topleft"))

ROC_prediction2 <- predict(model_CMA2, newdata=dfPts_mars2, type="response")
par(pty = "s")
#as.numeric(test_roc$auc)
ROC_plot2<- roc(dfPts_mars2$status2, ROC_prediction2, plot=TRUE, print.auc=TRUE, main="ROC curve", col="black")
coords(ROC_plot2, "best", ret=c("threshold", "specificity","sensitivity"), best.method=c("closest.topleft"))


