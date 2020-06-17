library(mice)
library(VIM)
library(survival)
library(cmprsk)
library(readxl)
library(survminer)
library(gmodels)

#change this directory, results and data will be saved in this directory. it will create the Figure and table folder/dicrectory
direc='~/Dropbox/MRI_Data/2020version/'
Figures=paste(direc,'Figures/',sep='')
Figure1_path=paste(direc,'Figures/Figure1/',sep='')
Figure2_path=paste(direc,'Figures/Figure2/',sep='')
Figure3_path=paste(direc,'Figures/Figure3/',sep='')
FigureS1_path=paste(direc,'Figures_raw/FigureS1/',sep='')
FigureS2_path=paste(direc,'Figures_raw/FigureS2/',sep='')
FigureS3_path=paste(direc,'Figures_raw/FigureS3/',sep='')


Report_path=paste(direc,'Report/',sep='')
dir.create(Figures)
dir.create(Figure1_path)
dir.create(Figure2_path)
dir.create(Figure3_path)
dir.create(Report_path)
dir.create(FigureTry_path)

#define a function to draw figures: 
drawFigure_l<-function(Figure_path,data,figurename){
  tiff(paste(Figure_path,figurename,sep=''), width = 1900, height = 1200, units = 'px', res = 250)
  s_e_o_s.surv <- survfit(Surv(time, status) ~ MRI, data = data,conf.type = "log-log") 
  p<-ggsurvplot(s_e_o_s.surv, data = data, pval = TRUE, risk.table = TRUE,risk.table.title='No. at risk' ,legend.labs=c("No MRI","MRI"),linetype =  "strata",
                legend.title='', xlab='Time (Years)',ylab=' Local Recurrence Free \n (Probability)',   
                font.main = c(22, "bold"),
                #tables.
                fontsize = 4,
                ylim = c(0.4, 1),
                font.x = c(18, "bold"),
                font.legend=c(16,"bold"),
                pval.coord=c(0.45,0.65),
                legend=c(0.75,0.25),
                pval.size=8,
                font.y = c(16, "bold"))
  p$table <- p$table + 
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.text = element_blank(),
      axis.line = element_blank(),
      title = element_text(size=3)
    )
  print(p)
  dev.off()
}

drawFigure_d<-function(Figure_path,data,figurename){
  tiff(paste(Figure_path,figurename,sep=''), width = 1900, height = 1200, units = 'px', res = 250)
  s_e_o_s.surv <- survfit(Surv(time, status) ~ MRI, data = data,conf.type = "log-log") 
  p<-ggsurvplot(s_e_o_s.surv, data = data, pval = TRUE, risk.table = TRUE,risk.table.title='No. at risk' ,legend.labs=c("No MRI","MRI"),linetype =  "strata",
             legend.title='', xlab='Time (Years)',ylab=' Distant Recurrence Free \n (Probability)',   
             font.main = c(22, "bold"),
             #tables.
             fontsize = 4,
             ylim = c(0.4, 1),
             font.x = c(18, "bold"),
             font.legend=c(16,"bold"),
             pval.coord=c(0.45,0.65),
             legend=c(0.75,0.25),
             pval.size=8,
             font.y = c(16, "bold"))
  p$table <- p$table + 
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.text = element_blank(),
      axis.line = element_blank(),
      title = element_text(size=3)
    )
  print(p)
  dev.off()
}

###################################################################################################################
#load the data, and remove un-validated data. product: dataframe: Merged_Table. containing the breast conservation data that annotated by Nat. 
##################################################################################################################
Big <- read_excel("~/Dropbox/MRI_Data/LSDB_2017_June18_Cleaned_No_meta_This_one_to_use_updated.xlsx",na = "NULL")
annotated<- read_excel("~/Dropbox/MRI_Data/LSDB_First_Annotation.xlsx", sheet = "Cleaned", skip = 1)
Merged_Table<-merge(annotated,Big,by.x='mrd_pt_id',by.y='mrd_pt_id' ,all.x=TRUE)
#some of the data need to removed, e.g. wrong diagnosis date, got mastectomy before
Merged_Table<-Merged_Table[which(is.na(Merged_Table$`mark 1 if the actual diagnosis is at least 1 year before this diagnosis date`) & is.na(Merged_Table$Delete)),]
#remove the ones have metastasis at diagnosis, this repsetns as no MRI data
Merged_Table<-Merged_Table[which(!is.na(Merged_Table$MRIs_60_surgery)),]
#remove non-invasive lobular, it is not cancer, per Dr. Khan
Merged_Table<-Merged_Table[which((Merged_Table$Histology2!='NC')),]
#try cut the value, only diagnoses after 2006-1-1
#put this into part two Merged_Table<-Merged_Table[which(as.POSIXct(Merged_Table$date_of_diagnosis_nat) > as.POSIXct("2006-01-01 00:00:00 UTC")),]

#let nat check again
##################################################################################################################
#retrieve the variables to be used in the MRI project, and remove the items that have followup days shorter than 90 days, product: dataframe: data 
##################################################################################################################
#make the table we need for further analysis
mrd_pt_id=Merged_Table$mrd_pt_id
name=Merged_Table$Name
birth=Merged_Table$birth_dt
diagnosis_date=Merged_Table$date_of_diagnosis_nat
MRI=Merged_Table$MRIs_60_surgery
age_at_diagnosis=(Merged_Table$date_of_diagnosis-Merged_Table$birth_dt)/365
race=ifelse(Merged_Table$race=='WHITE'& (Merged_Table$ethncty =='Not Hispanic' | is.na(Merged_Table$ethncty)) ,'Non-Hispanic whites',ifelse(Merged_Table$race=='BLACK'&(Merged_Table$ethncty =='Not Hispanic' | is.na(Merged_Table$ethncty)) ,'Non-Hispanic blacks',ifelse(Merged_Table$race=='Hispanic' | Merged_Table$ethncty=='Hispanic','Hispanic',ifelse(Merged_Table$race=='ASIAN','ASIAN',NA))))
race[is.na(race)] <- 'NR'
size=Merged_Table$SIZE
size[is.na(size)] <- mean(size,na.rm = TRUE)
grade=Merged_Table$GRADE
grade[is.na(grade)] <- 'NR'
lymph_status=Merged_Table$lymph_status
lymph_status[is.na(lymph_status)] <- 'NR'
histology=Merged_Table$Histology2
ER=Merged_Table$ER
ER=ifelse(ER=='LOWPOSITIVE','POSITIVE',ER)
ER[is.na(ER)] <- 'NR'
PR=Merged_Table$PR
PR=ifelse(PR=='LOWPOSITIVE','POSITIVE',PR)
PR[is.na(PR)] <- 'NR'
HER2=Merged_Table$HER2
HER2[is.na(HER2)] <- 'NR'
P53=Merged_Table$P53
P53=ifelse(P53=='LOWPOSITIVE','POSITIVE',P53)
P53[is.na(P53)] <- 'NR'
NEO<-Merged_Table$NEO
radiation=ifelse(Merged_Table$radiation_nat=='No','NO',ifelse(Merged_Table$radiation_nat=='YEs','YES',Merged_Table$radiation_nat))
chemo=ifelse(Merged_Table$any_chemo_nat=='No','NO',Merged_Table$any_chemo_nat)
endo=ifelse(Merged_Table$any_endo_nat=='No','NO',Merged_Table$any_endo_nat)
systematic_treatment=ifelse(chemo=='YES'|endo=='YES','YES','NO')
#get the outcome variables 
local_recurrence=Merged_Table$`local recurrence_nat`
local_recurrence=ifelse(local_recurrence=='YES',1,0)
Merged_Table$updated_local_day<-Merged_Table$`local recurrence date_nat`
Merged_Table$updated_local_day[which(is.na(Merged_Table$updated_local_day))]<-Merged_Table[which(is.na(Merged_Table$updated_local_day)),'last_contact_date']
local_days=unclass((Merged_Table$updated_local_day-Merged_Table$date_of_diagnosis))
distant_recurrence=ifelse(Merged_Table$`Distant recurrence_nat`=='No','NO',Merged_Table$`Distant recurrence_nat`)
distant_recurrence=ifelse(distant_recurrence=='YES',1,0)
Merged_Table$updated_distant_day<-Merged_Table$`Distant Recurrence Date_nat`
Merged_Table$updated_distant_day[which(is.na(Merged_Table$updated_distant_day))]<-Merged_Table[which(is.na(Merged_Table$updated_distant_day)),'last_contact_date']
distant_days=unclass((as.Date(Merged_Table$updated_distant_day)-as.Date(Merged_Table$date_of_diagnosis)))
survival=Merged_Table$Survival
survival=ifelse(survival=='Deceased',1,0)
survival_days=Merged_Table$Survival_Days
reason_noRadi=Merged_Table$reason_no_radiation
#defining a variable for distant competing for local recurrence (two scenarios, one is local did not happen, distant happens,censor at distant, aother is both happens, but distant first, censor at distnat )
local_distant_event=local_recurrence
local_distant_event[(local_recurrence==0 & distant_recurrence==1)|((local_recurrence==1&distant_recurrence==1)&(local_days>distant_days))] <-2
local_distant_days<-local_days
local_distant_days[local_distant_event==2]<-distant_days[local_distant_event==2]
#defining a variable for survival competing for local recurrence 
local_survival_event=local_recurrence
local_survival_event[local_recurrence==0 & survival==1] <-2
local_survival_days<-local_days
#defining a variable for distant recurernce, survival competing for local recurrence 
local_distant_survival_event=local_distant_event
local_distant_survival_event[local_distant_survival_event==0 & survival==1] <-3
local_distant_survival_days<-local_distant_days
#defining a variable for survival competing for distant recurrence 
distant_survival_event=distant_recurrence
distant_survival_event[distant_recurrence==0 & survival==1] <-2
distant_survival_days<-distant_days
#excision data
excision=(Merged_Table$ReExcision)
#mammography before diagnosis
mammo=Merged_Table$`total_mamography_before diagnosis`
#now, put the variables retried together into dataframe data, and retain only the items without events larger than 90 days, follow-up time less than 15 years (followup information not complete), and also less than 65 years old
data<-cbind.data.frame(mrd_pt_id,name,birth,diagnosis_date,MRI,age_at_diagnosis,race,size,grade,lymph_status,histology,NEO,ER,PR,HER2,P53,radiation,chemo,endo,systematic_treatment,local_recurrence,local_days,distant_recurrence,distant_days,survival,survival_days,local_distant_event,local_distant_days,local_survival_event,local_survival_days,local_distant_survival_event,local_distant_survival_days,distant_survival_event,distant_survival_days,excision,mammo,smoking=Merged_Table$smoking,alcohol_useage=Merged_Table$alcohol_useage,family_history=Merged_Table$family_history,reason_noRadi)
data<-data[which(data$local_days>= 90 & data$distant_days>= 90 &data$survival_days>= 90),]
data<-data[which(data$local_days<= 5475 & data$distant_days <=5475  &data$distant_days <=5475),]
data<-data[data$radiation=='YES',]
# should discuss if to include the age cut, or just adjust for it data<-data[which(data$age_at_diagnosis<=70),]

#this is for nat to recheck
data_death=data[data$race=='NR',]
write.csv(data_death,paste(direc,'grade_NR.csv',sep=''))
##################################################################################################################
#Below is to select out the sub-groups, and also to output the outcome events (local/distant recurrences, and death )
##################################################################################################################
#make dataframe for easier use later
#make sensitivity analysis, using only patients diagnosed after 2006, this is the time MRI get prevelent
sens_data_2006<-data[which(as.POSIXct(data$diagnosis_date) > as.POSIXct("2006-01-01 00:00:00 UTC")),]
sens_data_2006$age_binary<-ifelse(sens_data_2006$age_at_diagnosis<=50,'young','old')
CrossTable(sens_data_2006$MRI,sens_data_2006$age_binary,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

#select age smaller than 45
sens_data_age45<-sens_data_2006[sens_data_2006$age_at_diagnosis<=50,]
sens_data_2006=sens_data_age45
data=sens_data_2006
MRI_group<-data[data$MRI=='1',]
NO_MRI_group<-data[data$MRI=='0',]

#make the dataframe for future sensitivity analysis
sens_data_IDC<-sens_data_2006[which(sens_data_2006$histology=="IDC" | sens_data_2006$histology=="ILO"),]
#make another sensitity analysis using women with out radiation 
sens_data_no_radiation<-data[data$radiation=='NO',]
#make another sensitity analysis using women with out treatment 
sens_data_no_treatment<-data[data$radiation=='NO'& data$endo=='NO' &data$chemo=='NO',]
#make another sensitity analysis using women with followupday larger than three years 
sens_data_followup_3years<-sens_data_2006[sens_data_2006$survival_days>365*3,]
#make another sensitity analysis using women with followupday larger than three years 
sens_data_NEOs<-sens_data_2006[sens_data_2006$NEO=='NO',]

dim(sens_data_followup_3years)
dim(sens_data_IDC)
dim(sens_data_NEOs)

#reviewers: why patients not getting radiations 
table(MRI_group$radiation)
MRI_noradi=MRI_group[MRI_group$radiation=='NO',]
NOMRI_noradi=NO_MRI_group[NO_MRI_group$radiation=='NO',]
MRI_radi=MRI_group[MRI_group$radiation=='YES',]
NOMRI_radi=NO_MRI_group[NO_MRI_group$radiation=='YES',]

MRI_noradi_DCIS=MRI_group[MRI_group$radiation=='NO'&MRI_group$histology=='DCIS',]
NOMRI_noradi_DCIS=NO_MRI_group[NO_MRI_group$radiation=='NO'&NO_MRI_group$histology=='DCIS',]

Radiation=data[data$radiation=='YES',]
NoRadiation=data[data$radiation=='NO',]

fit=glm(MRI ~ age_at_diagnosis , data=data, family=binomial(link="logit"))
summary(fit)
#age
mean(Radiation$age_at_diagnosis)
sd(Radiation$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(Radiation$age_at_diagnosis[!is.na(Radiation$age_at_diagnosis)]))

mean(NoRadiation$age_at_diagnosis)
sd(NoRadiation$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(NoRadiation$age_at_diagnosis[!is.na(NoRadiation$age_at_diagnosis)]))

t.test(Radiation$age_at_diagnosis,NoRadiation$age_at_diagnosis)
#tumor size 
mean(Radiation$size)
sd(Radiation$size, na.rm=TRUE) /  
  sqrt(length(Radiation$size[!is.na(Radiation$size)]))

mean(NoRadiation$size)
sd(NoRadiation$size, na.rm=TRUE) /  
  sqrt(length(NoRadiation$size[!is.na(NoRadiation$size)]))
t.test(Radiation$size,NoRadiation$size)

#grade and lymph nodes 
table(data$radiation,data$grade)
table(data$radiation,data$lymph_status)
CrossTable(data$radiation,data$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(data$radiation,data$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

#below are for DCIS only: 
Radiation=data[data$radiation=='YES'&data$histology=='DCIS',]
NoRadiation=data[data$radiation=='NO'&data$histology=='DCIS',]
#age 
mean(Radiation$age_at_diagnosis)
sd(Radiation$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(Radiation$age_at_diagnosis[!is.na(Radiation$age_at_diagnosis)]))

mean(NoRadiation$age_at_diagnosis)
sd(NoRadiation$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(NoRadiation$age_at_diagnosis[!is.na(NoRadiation$age_at_diagnosis)]))

t.test(Radiation$age_at_diagnosis,NoRadiation$age_at_diagnosis)
#tumor size 
mean(Radiation$size)
sd(Radiation$size, na.rm=TRUE) /  
  sqrt(length(Radiation$size[!is.na(Radiation$size)]))

mean(NoRadiation$size)
sd(NoRadiation$size, na.rm=TRUE) /  
  sqrt(length(NoRadiation$size[!is.na(NoRadiation$size)]))
t.test(Radiation$size,NoRadiation$size)

#grade and lymph nodes 
DCIS=data[data$histology=='DCIS',]
table(data$radiation,data$grade)
table(data$radiation,data$lymph_status)
CrossTable(DCIS$radiation,DCIS$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(DCIS$radiation,DCIS$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

#below are the comments addressed by Dr. Khan
MRI_noradi=MRI_group[MRI_group$radiation=='NO',]
NOMRI_noradi=NO_MRI_group[NO_MRI_group$radiation=='NO',]
MRI_radi=MRI_group[MRI_group$radiation=='YES',]
NOMRI_radi=NO_MRI_group[NO_MRI_group$radiation=='YES',]

#age
mean(MRI_radi$age_at_diagnosis)
sd(MRI_radi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(MRI_radi$age_at_diagnosis[!is.na(MRI_radi$age_at_diagnosis)]))
mean(MRI_noradi$age_at_diagnosis)
sd(MRI_noradi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(MRI_noradi$age_at_diagnosis[!is.na(MRI_noradi$age_at_diagnosis)]))
mean(NOMRI_radi$age_at_diagnosis)
sd(NOMRI_radi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(NOMRI_radi$age_at_diagnosis[!is.na(NOMRI_radi$age_at_diagnosis)]))
mean(NOMRI_noradi$age_at_diagnosis)
sd(NOMRI_noradi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(NOMRI_noradi$age_at_diagnosis[!is.na(NOMRI_noradi$age_at_diagnosis)]))

t.test(MRI_noradi$age_at_diagnosis,MRI_radi$age_at_diagnosis)
t.test(NOMRI_noradi$age_at_diagnosis,NOMRI_radi$age_at_diagnosis)

t.test(MRI_radi$age_at_diagnosis,NOMRI_radi$age_at_diagnosis)
t.test(MRI_noradi$age_at_diagnosis,NOMRI_noradi$age_at_diagnosis)
#size
mean(MRI_radi$size)
sd(MRI_radi$size, na.rm=TRUE) /  
  sqrt(length(MRI_radi$size[!is.na(MRI_radi$size)]))
mean(MRI_noradi$size)
sd(MRI_noradi$size, na.rm=TRUE) /  
  sqrt(length(MRI_noradi$size[!is.na(MRI_noradi$size)]))
mean(NOMRI_radi$size)
sd(NOMRI_radi$size, na.rm=TRUE) /  
  sqrt(length(NOMRI_radi$size[!is.na(NOMRI_radi$size)]))
mean(NOMRI_noradi$size)
sd(NOMRI_noradi$size, na.rm=TRUE) /  
  sqrt(length(NOMRI_noradi$size[!is.na(NOMRI_noradi$size)]))

t.test(MRI_noradi$size,MRI_radi$size)
t.test(NOMRI_noradi$size,NOMRI_radi$size)
t.test(MRI_radi$size,NOMRI_radi$size)
t.test(MRI_noradi$size,NOMRI_noradi$size)
#follow up days
mean(MRI_radi$survival_days )/365
sd(MRI_radi$survival_days, na.rm=TRUE) /  
  sqrt(length(MRI_radi$survival_days[!is.na(MRI_radi$survival_days)]))/365
mean(MRI_noradi$survival_days)/365
sd(MRI_noradi$survival_days, na.rm=TRUE) /  
  sqrt(length(MRI_noradi$survival_days[!is.na(MRI_noradi$survival_days)]))/365
mean(NOMRI_radi$survival_days)/365
sd(NOMRI_radi$survival_days, na.rm=TRUE) /  
  sqrt(length(NOMRI_radi$survival_days[!is.na(NOMRI_radi$survival_days)]))/365
mean(NOMRI_noradi$survival_days)/365
sd(NOMRI_noradi$survival_days, na.rm=TRUE) /  
  sqrt(length(NOMRI_noradi$survival_days[!is.na(NOMRI_noradi$survival_days)]))/365

t.test(MRI_noradi$survival_days,MRI_radi$survival_days)
t.test(NOMRI_noradi$survival_days,NOMRI_radi$survival_days)
t.test(MRI_radi$survival_days,NOMRI_radi$survival_days)
t.test(MRI_noradi$survival_days,NOMRI_noradi$survival_days)
#grade and lymph nodes 
table(data$grade,data$radiation,data$MRI)
Rad=data[data$radiation =='YES',]
NORad=data[data$radiation =='NO',]

CrossTable(Rad$MRI,Rad$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(NORad$MRI,NORad$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

table(data$lymph_status,data$radiation,data$MRI)
CrossTable(Rad$MRI,Rad$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(NORad$MRI,NORad$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

#below are the comments addressed by Dr. Khan DCIS 
MRI_noradi=MRI_group[MRI_group$radiation=='NO' & MRI_group$histology=='DCIS',]
NOMRI_noradi=NO_MRI_group[NO_MRI_group$radiation=='NO'& NO_MRI_group$histology=='DCIS',]
MRI_radi=MRI_group[MRI_group$radiation=='YES'& MRI_group$histology=='DCIS',]
NOMRI_radi=NO_MRI_group[NO_MRI_group$radiation=='YES'& NO_MRI_group$histology=='DCIS',]

#age
mean(MRI_radi$age_at_diagnosis)
sd(MRI_radi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(MRI_radi$age_at_diagnosis[!is.na(MRI_radi$age_at_diagnosis)]))
mean(MRI_noradi$age_at_diagnosis)
sd(MRI_noradi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(MRI_noradi$age_at_diagnosis[!is.na(MRI_noradi$age_at_diagnosis)]))
mean(NOMRI_radi$age_at_diagnosis)
sd(NOMRI_radi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(NOMRI_radi$age_at_diagnosis[!is.na(NOMRI_radi$age_at_diagnosis)]))
mean(NOMRI_noradi$age_at_diagnosis)
sd(NOMRI_noradi$age_at_diagnosis, na.rm=TRUE) /  
  sqrt(length(NOMRI_noradi$age_at_diagnosis[!is.na(NOMRI_noradi$age_at_diagnosis)]))

t.test(MRI_noradi$age_at_diagnosis,MRI_radi$age_at_diagnosis)
t.test(NOMRI_noradi$age_at_diagnosis,NOMRI_radi$age_at_diagnosis)

t.test(MRI_radi$age_at_diagnosis,NOMRI_radi$age_at_diagnosis)
t.test(NOMRI_noradi$age_at_diagnosis,MRI_noradi$age_at_diagnosis)
#size
mean(MRI_radi$size)
sd(MRI_radi$size, na.rm=TRUE) /  
  sqrt(length(MRI_radi$size[!is.na(MRI_radi$size)]))
mean(MRI_noradi$size)
sd(MRI_noradi$size, na.rm=TRUE) /  
  sqrt(length(MRI_noradi$size[!is.na(MRI_noradi$size)]))
mean(NOMRI_radi$size)
sd(NOMRI_radi$size, na.rm=TRUE) /  
  sqrt(length(NOMRI_radi$size[!is.na(NOMRI_radi$size)]))
mean(NOMRI_noradi$size)
sd(NOMRI_noradi$size, na.rm=TRUE) /  
  sqrt(length(NOMRI_noradi$size[!is.na(NOMRI_noradi$size)]))

t.test(MRI_noradi$size,MRI_radi$size)
t.test(NOMRI_noradi$size,NOMRI_radi$size)
t.test(MRI_radi$size,NOMRI_radi$size)
t.test(NOMRI_noradi$size,MRI_noradi$size)
#follow up days
mean(MRI_radi$survival_days )/365
sd(MRI_radi$survival_days, na.rm=TRUE) /  
  sqrt(length(MRI_radi$survival_days[!is.na(MRI_radi$survival_days)]))/365
mean(MRI_noradi$survival_days)/365
sd(MRI_noradi$survival_days, na.rm=TRUE) /  
  sqrt(length(MRI_noradi$survival_days[!is.na(MRI_noradi$survival_days)]))/365
mean(NOMRI_radi$survival_days)/365
sd(NOMRI_radi$survival_days, na.rm=TRUE) /  
  sqrt(length(NOMRI_radi$survival_days[!is.na(NOMRI_radi$survival_days)]))/365
mean(NOMRI_noradi$survival_days)/365
sd(NOMRI_noradi$survival_days, na.rm=TRUE) /  
  sqrt(length(NOMRI_noradi$survival_days[!is.na(NOMRI_noradi$survival_days)]))/365

t.test(MRI_noradi$survival_days,MRI_radi$survival_days)
t.test(NOMRI_noradi$survival_days,NOMRI_radi$survival_days)
#grade and lymph nodes 
table(data$grade,data$radiation,data$MRI)
CrossTable(MRI_group$radiation,MRI_group$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(NO_MRI_group$radiation,NO_MRI_group$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

table(data$lymph_status,data$radiation,data$MRI)
CrossTable(MRI_group$radiation,MRI_group$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(NO_MRI_group$radiation,NO_MRI_group$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

#run code below to get the number events 
  #---local recurrence 
CrossTable(data$MRI,data$local_recurrence,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$local_recurrence,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
  #---distant recurrence 
CrossTable(data$MRI,data$distant_recurrence,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$distant_recurrence,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
#---survival 
CrossTable(data$MRI,data$survival,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_2006$MRI,sens_data_2006$survival,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
#---excision 
CrossTable(data$MRI,data$excision,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

#mammogragpy
mean(MRI_group$mammo, na.rm = TRUE)
mean(NO_MRI_group$mammo, na.rm = TRUE)
sd(MRI_group$mammo, na.rm = TRUE)
sd(NO_MRI_group$mammo, na.rm = TRUE)
t.test(MRI_group$mammo,NO_MRI_group$mammo)
(176)/330
(90)/182
#---followup-time
mean(MRI_group$survival_days)/365
mean(NO_MRI_group$survival_days)/365
sd(MRI_group$survival_days, na.rm = TRUE)/365
sd(NO_MRI_group$survival_days, na.rm = TRUE)/365
t.test(MRI_group$survival_days,NO_MRI_group$survival_days)
mean(sens_data_2006[which(sens_data_2006$MRI==1),]$survival_days)/365
mean(sens_data_2006[which(sens_data_2006$MRI==0),]$survival_days)/365
sd(sens_data_2006[which(sens_data_2006$MRI==1),]$survival_days, na.rm = TRUE)/365
sd(sens_data_2006[which(sens_data_2006$MRI==0),]$survival_days, na.rm = TRUE)/365
t.test(sens_data_2006[which(sens_data_2006$MRI==1),]$survival_days,sens_data_2006[which(sens_data_2006$MRI==0),]$survival_days)
##################################################################################################################
#Run code below to get table 1
###################################################################################################################
#--smoking
CrossTable(sens_data_age45$local_recurrence,sens_data_age45$family_history,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

  #--age
t.test(MRI_group$age_at_diagnosis,NO_MRI_group$age_at_diagnosis)    
sd(MRI_group$age_at_diagnosis) 
sd(NO_MRI_group$age_at_diagnosis)
t.test(sens_data_2006[which(sens_data_2006$MRI==1),]$age_at_diagnosis,sens_data_2006[which(sens_data_2006$MRI==0),]$age_at_diagnosis)    
sd(sens_data_2006[which(sens_data_2006$MRI==1),]$age_at_diagnosis) 
sd(sens_data_2006[which(sens_data_2006$MRI==0),]$age_at_diagnosis)
  #--race-ethnity
CrossTable(data$MRI,data$race,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$race,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(data[which(data$race!='NR'),]$MRI,data[which(data$race!='NR'),]$race,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

  #--tumor size
t.test(MRI_group$size,NO_MRI_group$size)   
sd(MRI_group$size,na.rm=TRUE) 
sd(NO_MRI_group$size,na.rm=TRUE)
t.test(sens_data_2006[which(sens_data_2006$MRI==1),]$size,sens_data_2006[which(sens_data_2006$MRI==0),]$size)   
sd(sens_data_2006[which(sens_data_2006$MRI==1),]$size,na.rm=TRUE) 
sd(sens_data_2006[which(sens_data_2006$MRI==0),]$size,na.rm=TRUE)
  #--grade
CrossTable(data$MRI,data$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$grade,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
#--histology
CrossTable(data$MRI,data$histology,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$histology,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

CrossTable(sens_data_2006$distant_recurrence,sens_data_2006$histology,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )

# --nodal status 
CrossTable(data$MRI,data$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$lymph_status,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --ER status 
CrossTable(data$MRI,data$ER,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$ER,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --PR status 
CrossTable(data$MRI,data$PR,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$PR,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --HER2 status 
CrossTable(data$MRI,data$HER2,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$HER2,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --P53 status 
CrossTable(data$MRI,data$P53,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$P53,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --radiation 
CrossTable(data$MRI,data$radiation,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$radiation,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --endo 
CrossTable(data$MRI,data$endo,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$endo,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
#-- chemo
CrossTable(data$MRI,data$chemo,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$chemo,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
#-- systematic
CrossTable(data$MRI,data$systematic_treatment,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(sens_data_age45$MRI,sens_data_age45$systematic_treatment,prop.t=FALSE, prop.r=TRUE, prop.c=FALSE,chisq=TRUE,prop.chisq=FALSE )
# --diagnosese date 
summary(MRI_group$diagnosis_date)
summary(NO_MRI_group$diagnosis_date)
summary(sens_data_age45[which(sens_data_age45$MRI==1),]$diagnosis_date)
summary(sens_data_age45[which(sens_data_age45$MRI==0),]$diagnosis_date)
# --fowwlow-up-time 
mean(MRI_group$survival_days, na.rm = TRUE)/365
mean(NO_MRI_group$survival_days, na.rm = TRUE)/365
sd(MRI_group$survival_days, na.rm = TRUE)/365
sd(NO_MRI_group$survival_days, na.rm = TRUE)/365
t.test(MRI_group$survival_days,NO_MRI_group$survival_days)
mean(sens_data_age45[which(sens_data_age45$MRI==1),]$survival_days, na.rm = TRUE)/365
mean(sens_data_age45[which(sens_data_age45$MRI==0),]$survival_days, na.rm = TRUE)/365
sd(sens_data_age45[which(sens_data_age45$MRI==1),]$survival_days, na.rm = TRUE)/365
sd(sens_data_age45[which(sens_data_age45$MRI==0),]$survival_days, na.rm = TRUE)/365
t.test(sens_data_age45[which(sens_data_age45$MRI==1),]$survival_days,sens_data_age45[which(sens_data_age45$MRI==0),]$survival_days)
#for the age
MRI_group_year<-cbind.data.frame(year=as.numeric(format(MRI_group$diagnosis_date,'%Y')),age=MRI_group$age_at_diagnosis,cate=1)
NO_MRI_group_year<-cbind.data.frame(year=as.numeric(format(NO_MRI_group$diagnosis_date,'%Y')),age=NO_MRI_group$age_at_diagnosis,cate=1)
MRI_group_year$cate=ifelse(MRI_group_year$year%in%c(2003,2004,2005),2, ifelse(MRI_group_year$year%in%c(2006,2007,2008),3,ifelse(MRI_group_year$year%in%c(2009,2010,2011),4,ifelse(MRI_group_year$year%in%c(2012,2013,2014),5,1))))
NO_MRI_group_year$cate=ifelse(NO_MRI_group_year$year%in%c(2003,2004,2005),2, ifelse(NO_MRI_group_year$year%in%c(2006,2007,2008),3,ifelse(NO_MRI_group_year$year%in%c(2009,2010,2011),4,ifelse(NO_MRI_group_year$year%in%c(2012,2013,2014),5,1))))

temp<-MRI_group_year[which(MRI_group_year$cate%in% c(3,4,5) ),]
res.aov <- aov(age ~ cate, data = temp)
summary(res.aov)
res.aov <- aov(age ~ cate, data = NO_MRI_group_year[which(MRI_group_year$cate%in% c(3,4,5) ),])
summary(res.aov)
temp1<-MRI_group_year[which(MRI_group_year$cate%in% c(3,3)),2]
temp2<-MRI_group_year[which(MRI_group_year$cate%in% c(4,4)),2]
t.test(temp1,temp2)

meanMRI_group_year<-aggregate(MRI_group_year$age, list(Region = MRI_group_year$cate), mean)
meanNOMRI_group_year<-aggregate(NO_MRI_group_year$age, list(Region = NO_MRI_group_year$cate), mean)

##################################################################################################################
#Run code below to get Figure 1 and Figure 2 (sensitivity analysis)
###################################################################################################################
  #--local recurrence
temp<-cbind.data.frame(status=data$local_recurrence,time=data$local_days/365,MRI=data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
#temp[3065,]$time=14.61917
summary(s_e_model)
drawFigure_l(Figure1_path,temp,'Local_recurrence_univariate.tif')
  #--distant recurrence
temp<-cbind.data.frame(status=data$distant_recurrence,time=data$distant_days/365,MRI=data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_d(Figure1_path,temp,'Distant_recurrence_univariate.tif')
####end of univariate 

######sensitivity analysis for Figure 2, this is the sensitivity using patients only after 2006. 
  #--local recurrence
temp<-cbind.data.frame(status=sens_data_age45$local_recurrence,time=sens_data_age45$local_days/365,MRI=sens_data_age45$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
head(order(temp$time,decreasing = TRUE))
temp[1140,]$time=10.10
temp[390,]$time=10.10
temp[1261,]
drawFigure_l(Figure1_path,temp,'Local_recurrence_univariate_ageLarger45.tif')
  #--distant recurrence
temp<-cbind.data.frame(status=sens_data_age45$distant_recurrence,time=sens_data_age45$distant_days/365,MRI=sens_data_age45$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_d(Figure1_path,temp,'Distant_recurrence_univariate_ageSmaller45.tif')

######sensitivity analysis for Figure S1, this is the sensitivity using patients excluding DCIS. 
#--local recurrence
temp<-cbind.data.frame(status=sens_data_IDC$local_recurrence,time=sens_data_IDC$local_days/365,MRI=sens_data_IDC$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_l(Figure1_path,temp,'Local_recurrence_IDC.tif')
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_IDC$distant_recurrence,time=sens_data_IDC$distant_days/365,MRI=sens_data_IDC$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_d(Figure1_path,temp,'Distant_recurrence_IDC.tif')

######sensitivity analysis for Figure S2, this is the sensitivity excluding less than 3 years . 
#--local recurrence
temp<-cbind.data.frame(status=sens_data_followup_3years$local_recurrence,time=sens_data_followup_3years$local_days/365,MRI=sens_data_followup_3years$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_l(Figure1_path,temp,'Local_recurrence_3year.tif')
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_followup_3years$distant_recurrence,time=sens_data_followup_3years$distant_days/365,MRI=sens_data_followup_3years$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_d(Figure1_path,temp,'Distant_recurrence_3year.tif')

######sensitivity analysis for Figure S2, this is the sensitivity excluding NEO . 
#--local recurrence
temp<-cbind.data.frame(status=sens_data_NEOs$local_recurrence,time=sens_data_NEOs$local_days/365,MRI=sens_data_NEOs$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_l(FigureS3_path,temp,'Local_recurrence_univariate_sensityvity_anlysis_NEO.tif')
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_NEOs$distant_recurrence,time=sens_data_NEOs$distant_days/365,MRI=sens_data_NEOs$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure_d(FigureS3_path,temp,'Distant_recurrence_univariate_sensityvity_anlysis_NEO.tif')

####end of univariate/sensitivity analysis 

######sensitivity analysis for Figure 2, this only uses the patient without radiation, chemo, endo. 
#--local recurrence
temp<-cbind.data.frame(status=sens_data_no_treatment$local_recurrence,time=sens_data_no_treatment$local_days/365,MRI=sens_data_no_treatment$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure(Figure2_path,temp,'Local_recurrence_univariate_sensityvity_anlysis_no_treatment.tif')
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_no_treatment$distant_recurrence,time=sens_data_no_treatment$distant_days/365,MRI=sens_data_no_treatment$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
drawFigure(Figure2_path,temp,'Distant_recurrence_univariate_sensityvity_anlysis_no_treatment.tif')
####end of univariate/sensitivity analysis 


######sensitivity analysis for Figure 2, this only uses the patient of IDC only  
#--local recurrence
temp<-cbind.data.frame(status=sens_data_IDC$local_recurrence,time=sens_data_IDC$local_days/365,MRI=sens_data_IDC$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_IDC$distant_recurrence,time=sens_data_IDC$distant_days/365,MRI=sens_data_IDC$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
######sensitivity analysis for Figure 2, this only uses the patient of followup for three years only  
#--local recurrence
temp<-cbind.data.frame(status=sens_data_followup_3years$local_recurrence,time=sens_data_followup_3years$local_days/365,MRI=sens_data_followup_3years$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_followup_3years$distant_recurrence,time=sens_data_followup_3years$distant_days/365,MRI=sens_data_followup_3years$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
######sensitivity analysis for Figure 2, this only uses the patient of followup for neo
#--local recurrence
temp<-cbind.data.frame(status=sens_data_NEOs$local_recurrence,time=sens_data_NEOs$local_days/365,MRI=sens_data_NEOs$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#--distant recurrence
temp<-cbind.data.frame(status=sens_data_NEOs$distant_recurrence,time=sens_data_NEOs$distant_days/365,MRI=sens_data_NEOs$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
####end of univariate/sensitivity analysis 
##################################################################################################################
#Run code below to get Table 2, this is for multivariate analysis, using local recurrence as outcome, and also includes the sensitiviy analysis
###################################################################################################################
#multivariate, local reucrrence
data=sens_data_age45
temp<-cbind.data.frame(status=data$local_recurrence,time=data$local_days/365,data)
#use hispanic as reference 
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
temp <- within(temp, histology <- relevel(histology, ref = 'IDC'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=temp)
summary<-cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])
summary_report<-cbind.data.frame(HR=paste(round(summary$`exp(coef)`, digits = 2),' (',round(summary$`lower .95`, digits = 2),', ',round(summary$`upper .95`, digits = 2),')',sep=''),P_value=round(summary$`summary(s_e_model)$coefficients[, 5]`, digits = 2))
rownames(summary_report)<-rownames(summary)
write.table(summary_report,paste(Report_path,'Local_recurrence_multivariate.txt',sep=''),sep = '\t')
#multivariate, local recurrence, only diagnosed after 2006
temp<-cbind.data.frame(status=sens_data_age45$local_recurrence,time=sens_data_age45$local_days/365,sens_data_age45)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
temp <- within(temp, histology <- relevel(histology, ref = 'IDC'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
summary<-cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])
summary_report<-cbind.data.frame(HR=paste(round(summary$`exp(coef)`, digits = 2),' (',round(summary$`lower .95`, digits = 2),', ',round(summary$`upper .95`, digits = 2),')',sep=''),P_value=round(summary$`summary(s_e_model)$coefficients[, 5]`, digits = 2))
rownames(summary_report)<-rownames(summary)
write.table(summary_report,paste(Report_path,'Local_recurrence_multivariate_age.txt',sep=''),sep = '\t')
##################################################################################################################
#Run code below to get Table 3, this is for multivariate analysis, using distant recurrence as outcome, and also includes the sensitiviy analysis
###################################################################################################################
#multivariate, distant reucrrence
temp<-cbind.data.frame(status=data$distant_recurrence,time=data$distant_days/365,data)
#use hispanic as reference 
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
temp <- within(temp, histology <- relevel(histology, ref = 'IDC'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=temp)
summary<-cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])
summary_report<-cbind.data.frame(HR=paste(round(summary$`exp(coef)`, digits = 2),' (',round(summary$`lower .95`, digits = 2),', ',round(summary$`upper .95`, digits = 2),')',sep=''),P_value=round(summary$`summary(s_e_model)$coefficients[, 5]`, digits = 2))
rownames(summary_report)<-rownames(summary)
write.table(summary_report,paste(Report_path,'Distant_recurrence_multivariate.txt',sep=''),sep = '\t')
#multivariate, distnat recurrence, only diagnosed after 2006
temp<-cbind.data.frame(status=sens_data_age45$distant_recurrence,time=sens_data_age45$distant_days/365,sens_data_age45)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
temp <- within(temp, histology <- relevel(histology, ref = 'IDC'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
summary<-cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])
summary_report<-cbind.data.frame(HR=paste(round(summary$`exp(coef)`, digits = 2),' (',round(summary$`lower .95`, digits = 2),', ',round(summary$`upper .95`, digits = 2),')',sep=''),P_value=round(summary$`summary(s_e_model)$coefficients[, 5]`, digits = 2))
rownames(summary_report)<-rownames(summary)
write.table(summary_report,paste(Report_path,'Distant_recurrence_multivariate_age.txt',sep=''),sep = '\t')

##################################################################################################################
#Below are the data to compare local/distant with Disease free groups 
###################################################################################################################
#--local recurrence
DF_data<-data[which(data$local_recurrence==0 & data$distant_recurrence ==0),]
dim(DF_data)
Local_data<-data[which(data$local_recurrence==1 ),]
dim(Local_data)
Local_Only_data<-data[which(data$local_recurrence==1 & data$distant_recurrence ==0),]
dim(Local_Only_data)
Distant_data<-data[which(data$distant_recurrence==1),]
dim(Distant_data)
Distant_Only_data<-data[which(data$local_recurrence==0 &data$distant_recurrence==1),]
dim(Distant_Only_data)
Local_Distnat_data<-data[which(data$local_recurrence==1 & data$distant_recurrence ==1),]
dim(Local_Distnat_data)
#start trying for different combinations
#try compare local recurrence with DF
Temp_data<-rbind.data.frame(DF_data,Local_data)
  #-local
temp<-cbind.data.frame(status=Temp_data$local_recurrence,time=Temp_data$local_days/365,MRI=Temp_data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
  #within local, now try to control for other covariates 
temp<-cbind.data.frame(status=Temp_data$local_recurrence,time=Temp_data$local_days/365,Temp_data)
temp <- within(temp, race <- relevel(race, ref = 'Hispanic'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+PR+HER2+radiation+systematic_treatment,data=temp)
summary<-cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])
summary_report<-cbind.data.frame(HR=paste(round(summary$`exp(coef)`, digits = 2),' (',round(summary$`lower .95`, digits = 2),', ',round(summary$`upper .95`, digits = 2),')',sep=''),P_value=round(summary$`summary(s_e_model)$coefficients[, 5]`, digits = 2))
rownames(summary_report)<-rownames(summary)
summary_report


#within local, but after year of 2006
DF_data<-data[which(sens_data_2006$local_recurrence==0 & sens_data_2006$distant_recurrence ==0),]
dim(DF_data)
Local_data<-data[which(sens_data_2006$local_recurrence==1 ),]
dim(Local_data)
#start trying for different combinations
#try compare local recurrence with DF
Temp_data<-rbind.data.frame(DF_data,Local_data)
#-local
temp<-cbind.data.frame(status=Temp_data$local_recurrence,time=Temp_data$local_days/365,MRI=Temp_data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#within local, now try to control for other covariates 
temp<-cbind.data.frame(status=Temp_data$local_recurrence,time=Temp_data$local_days/365,Temp_data)
temp <- within(temp, race <- relevel(race, ref = 'Hispanic'))s
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+PR+HER2+P53+radiation+systematic_treatment,data=temp)
summary<-cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])
summary_report<-cbind.data.frame(HR=paste(round(summary$`exp(coef)`, digits = 2),' (',round(summary$`lower .95`, digits = 2),', ',round(summary$`upper .95`, digits = 2),')',sep=''),P_value=round(summary$`summary(s_e_model)$coefficients[, 5]`, digits = 2))
rownames(summary_report)<-rownames(summary)
summary_report
#-Distant
temp<-cbind.data.frame(status=Temp_data$distant_recurrence,time=data$distant_days/365,MRI=Temp_data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)

drawFigure(FigureTry_path,temp,'Local_recurrence_univariate.tif')
#--distant recurrence

drawFigure(Figure1_path,temp,'Distant_recurrence_univariate.tif')
####end of univariate 

##################################################################################################################
#Below are the competing risks using distant recurrence and survival 
###################################################################################################################
#local recurrence
temp<-cbind.data.frame(status=data$local_recurrence,time=data$local_days/365,MRI=data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#local recurrence, introduce distant recurrence to compete
temp<-cbind.data.frame(status=data$local_distant_event,time=data$local_distant_days/365,MRI=data$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                         fstatus = temp$status,  # variable with distinct codes for different causes of failure
                         cov1  = MRI_input,  # estimates will calculated within groups
                         failcode = 1, # code of fstatus that denotes the failure type of interest
                         cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)

#local recurrence, introduce death  to compete
temp<-cbind.data.frame(status=data$local_survival_event,time=data$local_survival_days/365,MRI=data$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)
#local recurrence, introduce distant recurrence and death  to compete
temp<-cbind.data.frame(status=data$local_distant_survival_event,time=data$local_distant_survival_days/365,MRI=data$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)
#################below, still local recurrence, but use the women diagnosed after year of 2006
temp<-cbind.data.frame(status=sens_data_2006$local_recurrence,time=sens_data_2006$local_days/365,MRI=sens_data_2006$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#local recurrence, introduce distant recurrence to compete
temp<-cbind.data.frame(status=sens_data_2006$local_distant_event,time=sens_data_2006$local_distant_days/365,MRI=sens_data_2006$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)
#local recurrence, introduce death  to compete
temp<-cbind.data.frame(status=sens_data_2006$local_survival_event,time=sens_data_2006$local_survival_days/365,MRI=sens_data_2006$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)
#local recurrence, introduce distant recurrence and death  to compete
temp<-cbind.data.frame(status=sens_data_2006$local_distant_survival_event,time=sens_data_2006$local_distant_survival_days/365,MRI=sens_data_2006$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)

#distant recurrence
temp<-cbind.data.frame(status=data$distant_recurrence,time=data$distant_days/365,MRI=data$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
temp<-cbind.data.frame(status=sens_data_2006$distant_recurrence,time=sens_data_2006$distant_days/365,MRI=sens_data_2006$MRI)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI),data=temp)
summary(s_e_model)
#distnat reucrrence, use death as compete
temp<-cbind.data.frame(status=data$distant_survival_event,time=data$distant_survival_days/365,MRI=data$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)
#distnat reucrrence sensitivity
temp<-cbind.data.frame(status=sens_data_2006$distant_recurrence,time=sens_data_2006$local_distant_survival_days/365,MRI=sens_data_2006$MRI)
MRI_input <- matrix(as.numeric(temp$MRI == "1"))
colnames(MRI_input) <- "dis"
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
summary(resCumIncByDis)

##now start doing multivariate####################
#local recurrence 
temp<-cbind.data.frame(status=data$local_recurrence,time=data$local_days/365,data)
#use hispanic as reference 
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
temp <- within(temp, histology <- relevel(histology, ref = 'IDC'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,], digits = 2)

data <- within(data, race <- relevel(race, ref = 'Non-Hispanic whites'))
data <- within(data, histology <- relevel(histology, ref = 'IDC'))

#local recurrence, introduce distant recurrence to compete
temp<-cbind.data.frame(status=data$local_distant_event,time=data$local_distant_days/365,MRI=data$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=data)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)

#local recurrence, introduce death  to compete
temp<-cbind.data.frame(status=data$local_survival_event,time=data$local_survival_days/365,MRI=data$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=data)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)
#local recurrence, introduce distant recurrence and death  to compete
temp<-cbind.data.frame(status=data$local_distant_survival_event,time=data$local_distant_survival_days/365,MRI=data$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=data)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)
#################below, still local recurrence, but use the women diagnosed after year of 2006
temp<-cbind.data.frame(status=sens_data_2006$local_recurrence,time=sens_data_2006$local_days/365,sens_data_2006)
#use hispanic as reference 

s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,], digits = 2)

#local recurrence, introduce distant recurrence to compete
temp<-cbind.data.frame(status=sens_data_2006$local_distant_event,time=sens_data_2006$local_distant_days/365,MRI=sens_data_2006$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=sens_data_2006)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)
#local recurrence, introduce death  to compete
temp<-cbind.data.frame(status=sens_data_2006$local_survival_event,time=sens_data_2006$local_survival_days/365,MRI=sens_data_2006$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=sens_data_2006)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)
#local recurrence, introduce distant recurrence and death  to compete
temp<-cbind.data.frame(status=sens_data_2006$local_distant_survival_event,time=sens_data_2006$local_distant_survival_days/365,MRI=sens_data_2006$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=sens_data_2006)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)

#distant recurrence
temp<-cbind.data.frame(status=data$distant_recurrence,time=data$distant_days/365,data)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,],2)

temp<-cbind.data.frame(status=sens_data_2006$distant_recurrence,time=sens_data_2006$distant_days/365,sens_data_2006)
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,],2)
#distnat reucrrence, use death as compete
temp<-cbind.data.frame(status=data$distant_survival_event,time=data$distant_survival_days/365,MRI=data$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+systematic_treatment,data=data)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)
#distnat reucrrence sensitivity
temp<-cbind.data.frame(status=sens_data_2006$distant_recurrence,time=sens_data_2006$local_distant_survival_days/365,MRI=sens_data_2006$MRI)
MRI_input <- model.matrix(~as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=sens_data_2006)[, -1]
resCumIncByDis <- crr(ftime   = temp$time,  # failure time variable
                      fstatus = temp$status,  # variable with distinct codes for different causes of failure
                      cov1  = MRI_input,  # estimates will calculated within groups
                      failcode = 1, # code of fstatus that denotes the failure type of interest
                      cencode  = 0 # value of fstatus variable which indicates the failure time is censored.
)
round(cbind.data.frame(summary(resCumIncByDis)$conf[,c(1,3,4)],summary(resCumIncByDis)$coef[,5])[1,],2)

#################below, is for the other types of sensitivity analysis 
#######################follow up three years 
#local recurrence
temp<-cbind.data.frame(status=sens_data_followup_3years$local_recurrence,time=sens_data_followup_3years$local_days/365,sens_data_followup_3years)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,], digits = 2)
#distant recurrence 
temp<-cbind.data.frame(status=sens_data_followup_3years$distant_recurrence,time=sens_data_followup_3years$distant_days/365,sens_data_followup_3years)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,],2)
#######################
#######################IDC
#local recurrence
temp<-cbind.data.frame(status=sens_data_IDC$local_recurrence,time=sens_data_IDC$local_days/365,sens_data_IDC)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,], digits = 2)
#distant recurrence 
temp<-cbind.data.frame(status=sens_data_IDC$distant_recurrence,time=sens_data_IDC$distant_days/365,sens_data_IDC)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,],2)
#######################
#######################NEOs 
#local recurrence
temp<-cbind.data.frame(status=sens_data_NEOs$local_recurrence,time=sens_data_NEOs$local_days/365,sens_data_NEOs)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,], digits = 2)
#distant recurrence 
temp<-cbind.data.frame(status=sens_data_NEOs$distant_recurrence,time=sens_data_NEOs$distant_days/365,sens_data_NEOs)
temp <- within(temp, race <- relevel(race, ref = 'Non-Hispanic whites'))
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(MRI)+age_at_diagnosis+race+size+grade+lymph_status+histology+ER+HER2+P53+radiation+systematic_treatment,data=temp)
round(cbind.data.frame(summary(s_e_model)$conf.int[,c(1,3,4)],summary(s_e_model)$coefficients[,5])[1,],2)
#######################

##################################################################################################################
#Below are some old codes that are not used in the analysis any more 
###################################################################################################################

temp1<-new_data
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+age_at_diagnosis+race2+SIZE+GRADE+lymph_all+Histology2+ER+PR+HER2+Ki67_2+P53+radiation+Systematic_therapy,data=temp1)
summary(s_e_model)

#table 2 new sensitivity
temp1<-sens_new_data
temp1$status<-ifelse(temp1$Local_Recurrence=='NO',0,1)
temp1$time<-temp1$Local_Recurrence_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+age_at_diagnosis+race2+SIZE+GRADE+lymph_all+ER+PR+HER2+Ki67_2+P53+radiation+Systematic_therapy,data=temp1)
summary(s_e_model)

temp1<-sens_new_data
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~as.factor(type)+age_at_diagnosis+race2+SIZE+GRADE+lymph_all+ER+PR+HER2+Ki67_2+P53+radiation+Systematic_therapy,data=temp1)
summary(s_e_model)



##
#test interaction
temp1<-new_data
temp1$status<-ifelse(temp1$Local_Recurrence=='NO',0,1)
temp1$time<-temp1$Local_Recurrence_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)*Histology2,data=temp1)
summary(s_e_model)

temp1<-new_data
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)*Histology2,data=temp1)
summary(s_e_model)




################################################################################################################ old datathat are no more used 

#For table 7:
table(MRI,new_data$race2)

table(MRI_group$race2,MRI_group$Survival )
table(MRI_group$race2,MRI_group$Local_Recurrence )
table(MRI_group$race2,MRI_group$Distant_Recurrence )

table(NO_MRI_group$race2,NO_MRI_group$Survival )
table(NO_MRI_group$race2,NO_MRI_group$Local_Recurrence )
table(NO_MRI_group$race2,NO_MRI_group$Distant_Recurrence )
################################################################################################################

############################################################################################################################################
#for table 3 and figure 1.
####################################################################################
temp1<-new_data
temp1$status<-ifelse(temp1$Survival=='Alive',0,1)
temp1$time<-temp1$Survival_Days

temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type),new_data=temp1)
summary(s_e_model)
table(temp1$status)
table(temp1$status,temp1$MRIs_60_surgery)

temp1$time<-temp1$time/365
png(filename="/Users/zza847/Documents/MRI_Study/figures/figure1/Figure1_Survival.png",width = 700, height = 480)
s_e_o_s.surv <- survfit(Surv(time, status) ~ type, new_data = temp1,conf.type = "log-log") 
ggsurvplot(s_e_o_s.surv, new_data = temp1,pval=TRUE, risk.table = FALSE,legend.labs=c("No MRI","MRI"),
           legend.title='Class: ',title=' Survival ', xlab='Time / Years',ylab='Probability',   xlim = c(0, 11),
           font.main = c(22, "bold"),
           ylim = c(0.4, 1),
           font.x = c(18, "bold"),
           font.legend=c(16,"bold"),
           legend=c(0.35,0.15),
           pval.size=10,
           pval.coord=c(0.45,0.65),
           font.y = c(18, "bold")) 
dev.off()

#-----------------------local recurrence
temp1<-new_data
temp1$status<-ifelse(temp1$Local_Recurrence=='NO',0,1)
temp1$time<-temp1$Local_Recurrence_Days
temp1$type<-temp1$MRIs_60_surgery

s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type),new_data=temp1)
summary(s_e_model)
table(temp1$status)
table(temp1$MRIs_60_surgery,temp1$status)

temp1$time<-temp1$time/365
png(filename="/Users/zza847/Documents/MRI_Study/figures/figure1/Figure1_Local_recurrence.png",width = 700, height = 480)
s_e_o_s.surv <- survfit(Surv(time, status) ~ type, new_data = temp1,conf.type = "log-log") 
ggsurvplot(s_e_o_s.surv, new_data = temp1,pval=TRUE, risk.table = FALSE,legend.labs=c("No MRI","MRI"),
           legend.title='Class: ',title=' Local Recurrence ', xlab='Time / Years',ylab='Probability',   xlim = c(0, 11),
           font.main = c(22, "bold"),
           ylim = c(0.4, 1),
           font.x = c(18, "bold"),
           font.legend=c(16,"bold"),
           legend=c(0.35,0.15),
           pval.size=10,
           pval.coord=c(0.45,0.65),
           font.y = c(18, "bold")) 
dev.off()

#-----------------------distant recurrence
temp1<-new_data
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type),new_data=temp1)
summary(s_e_model)
table(temp1$status)
table(temp1$MRIs_60_surgery,temp1$status)

temp1$time<-temp1$time/365
png(filename="/Users/zza847/Documents/MRI_Study/figures/figure1/Figure1_distant_recurrence.png",width = 700, height = 480)
s_e_o_s.surv <- survfit(Surv(time, status) ~ type, new_data = temp1,conf.type = "log-log") 
ggsurvplot(s_e_o_s.surv, new_data = temp1,pval=TRUE, risk.table = FALSE,legend.labs=c("No MRI","MRI"),
           legend.title='Class: ',title=' Distant Recurrence ', xlab='Time / Years',ylab='Probability',   xlim = c(0, 11),
           font.main = c(22, "bold"),
           ylim = c(0.4, 1),
           font.x = c(18, "bold"),
           font.legend=c(16,"bold"),
           legend=c(0.35,0.15),
           pval.size=10,
           pval.coord=c(0.45,0.65),
           font.y = c(18, "bold")) 
dev.off()


################################################################################################################
# free try : 
#use only breast conservation 
temp1<-new_data[which(new_data$Conservation_Surgery=='Breast Conservation Surgery' ),]
temp1<-new_data[which(new_data$Conservation_Surgery=='MASTECTOMY' ),]
temp1<-new_data[which(new_data$Conservation_Surgery=='MASTECTOMY' & new_data$race2=="Other"),]
temp1<-new_data[which(new_data$Conservation_Surgery=='Breast Conservation Surgery' & new_data$race2=="Other"),]
temp1<-new_data[which( new_data$race2=="WHITE"),]

temp1$status<-ifelse(temp1$Survival=='Alive',0,1)
temp1$time<-temp1$Survival_Days
temp1$status<-ifelse(temp1$Local_Recurrence=='NO',0,1)
temp1$time<-temp1$Local_Recurrence_Days
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days

temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+GRADE+SIZE+age_at_diagnosis+Histology2+radiation+NEO+TNEG,new_data=temp1)   
summary(s_e_model)
table(temp1$status)
table(temp1$MRIs_60_surgery,temp1$status)

#use only mastectomy 
temp1<-new_data[which(new_data$Conservation_Surgery=='MASTECTOMY'),]

temp1$status<-ifelse(temp1$Survival=='Alive',0,1)
temp1$time<-temp1$Survival_Days
temp1$status<-ifelse(temp1$Local_Recurrence=='NO',0,1)
temp1$time<-temp1$Local_Recurrence_Days
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days

temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+GRADE+SIZE+age_at_diagnosis+Histology2+radiation+NEO,new_data=temp1)   
summary(s_e_model)
table(temp1$status)
table(temp1$MRIs_60_surgery,temp1$status)


####################################################################################
############################################################################################################################################
#for new table 4 
####################################################################################
temp1<-new_data
temp1$status<-ifelse(temp1$Survival=='Alive',0,1)
temp1$time<-temp1$Survival_Days

temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+GRADE+SIZE+age_at_diagnosis+Histology2+radiation+NEO+TNEG+race2+Conservation_Surgery,new_data=temp1)
summary(s_e_model)
table(temp1$status)
table(temp1$status,temp1$MRIs_60_surgery)


#-----------------------local recurrence
temp1<-new_data
temp1$status<-ifelse(temp1$Local_Recurrence=='NO',0,1)
temp1$time<-temp1$Local_Recurrence_Days
temp1$type<-temp1$MRIs_60_surgery

s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+GRADE+SIZE+age_at_diagnosis+Histology2+radiation+NEO+TNEG+race2+Conservation_Surgery,new_data=temp1)
summary(s_e_model)
table(temp1$status)
table(temp1$MRIs_60_surgery,temp1$status)

#-----------------------distant recurrence
temp1<-new_data
temp1$status<-ifelse(temp1$Distant_Recurrence=='NO',0,1)
temp1$time<-temp1$D_Days
temp1$type<-temp1$MRIs_60_surgery
s_e_model<-coxph(formula = Surv(time,status) ~ as.factor(type)+GRADE+SIZE+age_at_diagnosis+Histology2+radiation+NEO+TNEG+race2+Conservation_Surgery,new_data=temp1)
summary(s_e_model)
table(temp1$status)
table(temp1$MRIs_60_surgery,temp1$status)

####################################################################################
#####################################################################################en
#end

#follow up length: 
hist(new_data$Survival_Days,main='Follow-Length',xlab='Follow-Up Length / days',ylab='Number of patients',breaks = 10)
summary(new_data$Survival_Days)

#curve for survival
status_S<- ifelse(new_data$Survival=="Alive", 0, 1)
recu <- survfit(Surv(new_data$Survival_Days, status_S) ~ 1,  conf.type="log-log")
plot(recu, xlab="Time/Days", ylab="Survival",main="Survival Rate",lty=2:3 ,conf.int=TRUE,cex=0.2)

#Cureve for local recurrence
status_R<- ifelse(L_recur=="Yes", TRUE, FALSE)
recuR <- survfit(Surv(Local_recurrence_Days_Anki, status_R) ~1,  conf.type="log-log")
plot(recuR, xlab="Time/Days", ylab="No Local Recurrence",main="Local Recurrence Rate",lty=2:3 ,conf.int=TRUE,cex=0.2)

#Cureve for distant recurrence
status_M<- ifelse(D_recur=="Yes", TRUE, FALSE)
recuM <- survfit(Surv(Metastatis_Days_ANKI, status_M) ~1,  conf.type="log-log")
plot(recuM, xlab="Time/Days", ylab="No Distant Recurrence",main="Distant Recurrence Rate",lty=2:3 ,conf.int=TRUE,cex=0.2)

#MRI events and survival events together
model<- glm(status_S~MRI,family = "binomial")
summary(model)
#MRI event and lip
table(MRI, status_S)
table(status_R,MRI)
table(status_M,MRI)


chisq.test(table(status_R[Conservation_Surgery=='Breast Conservation Surgery'],MRI[Conservation_Surgery=='Breast Conservation Surgery']),correct=T)
chisq.test(table(status_R,MRI),correct=T)
chisq.test(table(MRI, status_M),correct=T)
chisq.test(table(MRI, status_S),correct=T)


# MRI events and surgery type
summary(Conservation_Surgery)
table(MRI, Conservation_Surgery)
Surgery<-as.factor(ifelse(Conservation_Surgery=='Breast Conservation Surgery',0,ifelse(Conservation_Surgery=='MASTECTOMY',1,NA)))
chisq.test(table(MRI, Surgery)[apply(table(MRI, Surgery), 1, function(x) all(x !=0 )),],correct=T)
model<- glm(Conservation_Surgery~MRI,family = "binomial")
summary(model)


#make table 2
## Calculate the overall CSH for mortality
#for survival 
resCshMortByDis <- survfit(formula   = Surv(Survival_Days, status_S) ~ 1, new_data= new_data, type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
resCshMortByDis <- survfit(formula   = Surv(Survival_Days, status_S) ~MRI, new_data= new_data, type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
summary(resCshMortByDis, times =c(0,1825))
survdiff(Surv(Survival_Days, status_S) ~ MRI,  rho=1) 

#for Ipsilateral recurrence event
resCshMortByDis <- survfit(formula= Surv(Local_recurrence_Days_Anki[Conservation_Surgery=='Breast Conservation Surgery'], status_R[Conservation_Surgery=='Breast Conservation Surgery']) ~ 1,  type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
resCshMortByDis <- survfit(formula= Surv(Local_recurrence_Days_Anki[Conservation_Surgery=='Breast Conservation Surgery'], status_R[Conservation_Surgery=='Breast Conservation Surgery']) ~ MRI[Conservation_Surgery=='Breast Conservation Surgery'],  type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
summary(status_R[Conservation_Surgery=='Breast Conservation Surgery'])
summary(resCshMortByDis, times =c(0,1825))
survdiff(Surv(Local_recurrence_Days_Anki[Conservation_Surgery=='Breast Conservation Surgery'], status_R[Conservation_Surgery=='Breast Conservation Surgery']) ~ MRI[Conservation_Surgery=='Breast Conservation Surgery'],  rho=1) 

#for chest wall recurrence event
resCshMortByDis <- survfit(formula= Surv(Local_recurrence_Days_Anki, status_R) ~ 1,  type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
resCshMortByDis <- survfit(formula= Surv(Local_recurrence_Days_Anki, status_R) ~ MRI,  type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
summary(status_R[Conservation_Surgery=='MASTECTOMY'])
summary(resCshMortByDis, times =c(0,1825))
survdiff(Surv(Local_recurrence_Days_Anki[Conservation_Surgery=='MASTECTOMY'], status_R[Conservation_Surgery=='MASTECTOMY']) ~ MRI[Conservation_Surgery=='MASTECTOMY'],  rho=1) 

#for distant recurrence event
resCshMortByDis <- survfit(formula   = Surv(Metastatis_Days_ANKI, status_M) ~ 1, new_data= new_data, type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
resCshMortByDis <- survfit(formula   = Surv(Metastatis_Days_ANKI, status_M) ~MRI, new_data= new_data, type = "kaplan-meier",error= "greenwood",conf.type = "log-log")
summary(status_M)
summary(resCshMortByDis, times =c(0,1825))
survdiff(Surv(Metastatis_Days_ANKI, status_M) ~ MRI,  rho=1) 


#Experiment 1:  Check if MRI predicts survival using survival model
#survival Rate:
status_S<- ifelse(Survival=="Alive", FALSE, TRUE)
recu <- survfit(Surv(Survival_Days, status_S) ~ MRI,  conf.type="log-log")
plot(recu, xlab="Time/Days", ylab="Survival Rate",main="Survival_MRI",lty=2:3 ,conf.int=TRUE,cex=0.2,col=2:3)
legend(100, .4, c("No", "Yes"), col=2:3, lty = 2:3) 
#survdiff(Surv(Survival_Days, status_S) ~ MRI,  rho=1) 
fitmodel<-coxph(Surv(Survival_Days, status_S) ~ MRI) 
summary(fitmodel)

#local recurrence survival study
status_R<- ifelse(L_recur=="Yes", TRUE, FALSE)
fitmodel<-coxph(Surv(Local_recurrence_Days_Anki, status_R) ~ MRI) 
summary(fitmodel)

#distant recurrence survival study
status_M<- ifelse(D_recur=="Yes", TRUE, FALSE)
fitmodel<-coxph(Surv(Metastatis_Days_ANKI, status_M) ~ MRI) 
summary(fitmodel)


#Experiment 2:  Check if MRI predicts survival using survival model, control other variables
#Recurrence Rate:
#Size_C<-as.numeric(levels(SIZE))[SIZE]
RACE_C<- as.factor(ifelse(race=='BLACK', 0,ifelse(race=='WHITE',1,2)))
Lymph_node_positive_F<-as.factor(ifelse(Lymph_node_positive==0,0, ifelse (0<Lymph_node_positive & Lymph_node_positive<=3,1, ifelse(4<=Lymph_node_positive&Lymph_node_positive<7,2, ifelse (8<=Lymph_node_positive&Lymph_node_positive<=10 ,3,ifelse(Lymph_node_positive>10,4, NA))))))
Lymph_node_positive_F<-as.factor(ifelse(Lymph_node_positive==0,0, ifelse (0<Lymph_node_positive , 1 ,NA)))
ER_F<-as.factor(ifelse(ER=='NEGATIVE',0, ifelse (ER=='LOWPOSITIVE',0, ifelse(ER=='POSITIVE',1, NA))))
ERPER_C<-as.numeric(as.numeric(sub("%","",ER_percent))/100)
PRPER_C<-as.numeric(as.numeric(sub("%","",PR_percent))/100)
PR_F<-as.factor(ifelse(PR=='NEGATIVE',0, ifelse (PR=='LOWPOSITIVE',1, ifelse(PR=='POSITIVE',1, NA))))
P53_F<-as.factor(ifelse(P53=='NEGATIVE',0, ifelse (P53=='LOWPOSITIVE',1, ifelse(P53=='POSITIVE',1, NA))))
HER2_F<-as.factor(ifelse(HER2=='NEGATIVE',0, ifelse (HER2=='POSITIVE',1,  NA)))
Ki67_F<-as.factor(ifelse(Ki67=='NEGATIVE',0, ifelse (Ki67=='LOW',1, ifelse(Ki67=='INTERMEDIATE',2,ifelse(Ki67=='HIGH', 3,NA)))))
Histology_F<-as.factor(ifelse(Histology=='DUCT',0, ifelse (Histology=='LOBULAR',1, ifelse(Histology=='MIXED DUCT AND LOBULAR       ',2,NA))))
Histology_2<-as.factor(ifelse(Histology2=='DCIS',0, ifelse (Histology2=='IDC',1, ifelse(Histology2=='ILO',1,NA))))
Invasive_F<-as.factor(ifelse(Invasive =='NO',0, ifelse(Invasive =='YES',1,  NA)))
GRADE_F<-as.factor(ifelse(GRADE==1 | GRADE==2, 0,ifelse(GRADE==3,1,NA)))
Surgery<-as.factor(ifelse(Conservation_Surgery=='Breast Conservation Surgery',0,ifelse(Conservation_Surgery=='MASTECTOMY',1,NA)))
status_S<- ifelse(Survival=="Alive", FALSE, TRUE)
#survdiff(Surv(Survival_Days, status_S) ~ MRI,  rho=1)
Radiation<-as.factor(ifelse(is.na(any_radiation), 0, 1))
breast_radi_F<-as.factor(ifelse(is.na(breast_radi),0, ifelse(breast_radi=='YES',1,NA)))
Chest_Wall_radi_F<-as.factor(ifelse(is.na(Chest_Wall_radi),0, ifelse(Chest_Wall_radi=='YES',1,NA)))
Nodal_Fields_radi_F<-as.factor(ifelse(is.na(Nodal_Fields_radi),0, ifelse(Nodal_Fields_radi=='YES',1,NA)))
Radiation_Unspecified_F<-as.factor(ifelse(is.na(Radiation_Unspecified),0, ifelse(Radiation_Unspecified=='YES',1,NA)))

Chemo<-as.factor(ifelse(any_chemo=='YES', 1, 0))
Endo<-as.factor(ifelse(any_endo=='YES', 1, 0))
NEOA<-as.factor(ifelse(NEO=='YES', 1, 0))
HER_2_Inhibitors_F<-as.factor(ifelse(is.na(HER_2_Inhibitors),0, ifelse(HER_2_Inhibitors=='YES',1,NA)))
Gonadotrophin_releasing_hormone_agonists_F<-as.factor(ifelse(is.na(Gonadotrophin_releasing_hormone_agonists),0, ifelse(Gonadotrophin_releasing_hormone_agonists=='YES',1,NA)))
targeted_therapy_F<-as.factor(ifelse(is.na(targeted_therapy),0, ifelse(targeted_therapy=='YES',1,NA)))
antihormonal_agents_F<-as.factor(ifelse(is.na(antihormonal_agents),0, ifelse(antihormonal_agents=='YES',1,NA)))
Alkylating_agents_F<-as.factor(ifelse(is.na(Alkylating_agents),0, ifelse(Alkylating_agents=='YES',1,NA)))
anthracyclines_F<-as.factor(ifelse(is.na(anthracyclines),0, ifelse(anthracyclines=='YES',1,NA)))
Antimetabolites_F<-as.factor(ifelse(is.na(Antimetabolites),0, ifelse(Antimetabolites=='YES',1,NA)))
anti_tubulin_F<-as.factor(ifelse(is.na(anti_tubulin),0, ifelse(anti_tubulin=='YES',1,NA)))

#study survival
fitmodel<-coxph(Surv(Survival_Days, status_S) ~
                  MRI+age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+HER2_F+Histology_2+
                  Radiation+Chemo+Endo
) 
summary(fitmodel)
#study recurrence
status_R<- ifelse(L_recur=="Yes", TRUE, FALSE)
fitmodel<-coxph(Surv(Local_recurrence_Days_Anki, status_R) ~ 
                  MRI+age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+HER2_F+Histology_2+
                  Radiation+Chemo+Endo
                ) 
summary(fitmodel)
#study distant recurrence
status_M<- ifelse(D_recur=="Yes", TRUE, FALSE)
fitmodel<-coxph(Surv(Metastatis_Days_ANKI, status_M) ~ 
                  MRI+age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+HER2_F+Histology_2+
                  Radiation+Chemo+Endo
                ) 
summary(fitmodel)

#Experiment 3: imputate the covariates using MICE
covariates<-cbind.new_data.frame(
 MRI,age_at_diagnosis,RACE_C,SIZE,Lymph_node_positive_F,GRADE_F,ER_F,PR_F,HER2_F,Histology_2,
    Radiation,Chemo,Endo
)
md.pattern(covariates)
imp_new_data<-mice(covariates)
head(complete(imp_new_data))
#study survival
fitmodel<-coxph(Surv(Survival_Days, status_S) ~ 
                  MRI+age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+HER2_F+Histology_2+
                  Radiation+Chemo+Endo
                ,new_data.frame(complete(imp_new_data,5))) 
summary(fitmodel)
#study recurrence
status_R<- ifelse(L_recur=="Yes", TRUE, FALSE)
fitmodel<-coxph(Surv(Local_recurrence_Days_Anki, status_R) ~ 
                  MRI+age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+HER2_F+Histology_2+
                  Radiation+Chemo+Endo
                ,new_data.frame(complete(imp_new_data,5))) 
summary(fitmodel)

#study distant recurrence
fitmodel<-coxph(Surv(Metastatis_Days_ANKI, status_M) ~ 
                  MRI+age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+HER2_F+Histology_2+
                  Radiation+Chemo+Endo
                ,new_data.frame(complete(imp_new_data,5))) 
summary(fitmodel)


#Experiment 6: use MRI as output variables, use all other variables as input variables. Study what is driving to MRI
covariates_6<-cbind.new_data.frame(
  age_at_diagnosis,RACE_C,SIZE,Lymph_node_positive_F,GRADE_F,ER_F,ERPER_C,PR_F,PRPER_C,P53_F,HER2_F,Histology_2
  
)
imp_new_data6<-mice(covariates_6)
fit<-with(imp_new_data6,glm(as.factor(MRI)~
                          age_at_diagnosis+RACE_C+SIZE+Lymph_node_positive_F+GRADE_F+ER_F+PR_F+P53_F+HER2_F+Histology_2
                        ,family = "binomial"))
summary(pool(fit))


#Experiment 7 Study how MRI interact with age to affect survival 
age_C<-as.factor(ifelse(age_at_diagnosis<40,0, ifelse (40<=age_at_diagnosis & age_at_diagnosis<50,1, ifelse(50<=age_at_diagnosis&age_at_diagnosis<60,2, ifelse (60<=age_at_diagnosis &age_at_diagnosis<70 ,3,ifelse (70<=age_at_diagnosis&age_at_diagnosis<80 ,4,ifelse (80<=age_at_diagnosis,5,NA)))))))
age_C<-as.factor(ifelse(age_at_diagnosis<45,0, ifelse (45<=age_at_diagnosis & age_at_diagnosis<70,1, ifelse(70<=age_at_diagnosis,2,3))))

#how age predict survival 
fitmodel<-coxph(Surv(Survival_Days, status_S) ~ age_at_diagnosis )
summary(fitmodel)

fitmodel<-coxph(Surv(Survival_Days, status_S) ~ age_C) 
summary(fitmodel)

fitmodel<-coxph(Surv(Survival_Days, status_S) ~ age_C+MRI+MRI:age_C) 
summary(fitmodel)


#experiment 8 
fitmodel<-coxph(Surv(Survival_Days, status_S) ~ GRADE_F+MRI+MRI:GRADE_F) 
summary(fitmodel)

fitmodel<-coxph(Surv(Local_recurrence_Days_Anki, status_R) ~ GRADE_F+MRI+MRI:GRADE_F) 
summary(fitmodel)


fitmodel<- glm(status_R~
                 GRADE_F+MRI+MRI:GRADE_F
               ,family = "binomial")
summary(fitmodel)


#experiment 9 test age under 62: 
#Experiment 2: In the age smaller than 62 group, test the association between MRI and survival (in survival model)
#using survival as output
new_data_small<-new_data[which(age_at_diagnosis)<=62,]
status_S<- ifelse(new_data_small$Survival=="Alive", FALSE, TRUE)
recu <- survfit(Surv(new_data_conservation$Survival_Days, status_S) ~ as.factor(new_data_conservation$MRIs_60_diagnosis), new_data = new_data_conservation, conf.type="log-log")
plot(recu, xlab="Time/Days", ylab="Survival Rate",main="Survival_MRI",lty=2:3 ,conf.int=TRUE,cex=0.2,col=2:3)
legend(100, .4, c("No", "Yes"), col=2:3, lty = 2:3) 
survdiff(Surv(new_data_conservation$Survival_Days, status_S) ~ as.factor(new_data_conservation$MRIs_60_diagnosis),  rho=1) 
fitmodel<-coxph(Surv(new_data_conservation$Survival_Days, status_S) ~ as.factor(new_data_conservation$MRIs_60_diagnosis)) 
summary(fitmodel)

# use recurrence as output, in this case, it is modeling the 
new_data_small<-new_data[which(age_at_diagnosis<=55 & Conservation_Surgery=='Breast Conservation Surgery')  ,]
status_R<- ifelse(new_data_small$L_recur=="Yes", TRUE, FALSE)
fitmodel<-coxph(Surv(new_data_small$Local_recurrence_Days_Anki, status_R) ~ as.factor(new_data_small$MRIs_60_diagnosis)) 
summary(fitmodel)


#experiemtn 9: 
library(reshape2)
library(plyr)
library(ggplot2)

age_C<-as.factor(ifelse(age_at_diagnosis<40,1, ifelse (40<=age_at_diagnosis & age_at_diagnosis<50,2, ifelse(50<=age_at_diagnosis&age_at_diagnosis<60,3, ifelse (60<=age_at_diagnosis &age_at_diagnosis<70 ,4,ifelse (70<=age_at_diagnosis,5,6))))))
status_SS<- ifelse(new_data$Survival=="Alive", 'Survival', 'Death')
x<-as.new_data.frame.matrix(table(age_C,status_SS))
id<-c('0~40','40~50','50~60','60~70','70~100')
x1<-cbind(id,x)
xlong1 <- ddply(melt(x1, id.vars = 'id'), .(id), mutate, prop = value / sum(value))
ggplot(xlong1, aes(x = id, y =prop, fill = variable)) + geom_bar(stat = 'identity')

status_SS<- ifelse(new_data$L_recur=="Yes", 'Chest wall recurrence', 'No recurrence')
x<-as.new_data.frame.matrix(table(age_C,status_SS))
id<-c('0~40','40~50','50~60','60~70','70~100')
x1<-cbind(id,x)
xlong1 <- ddply(melt(x1, id.vars = 'id'), .(id), mutate, prop = value / sum(value))
ggplot(xlong1, aes(x = id, y =prop, fill = variable)) + geom_bar(stat = 'identity')

status_SS<- ifelse(new_data$D_recur=="Yes", 'Distant recurrence', 'No recurrence')
x<-as.new_data.frame.matrix(table(age_C,status_SS))
id<-c('0~40','40~50','50~60','60~70','70~100')
x1<-cbind(id,x)
xlong1 <- ddply(melt(x1, id.vars = 'id'), .(id), mutate, prop = value / sum(value))
ggplot(xlong1, aes(x = id, y =prop, fill = variable)) + geom_bar(stat = 'identity')


