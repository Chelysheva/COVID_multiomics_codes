
# missing vaccine:final samples and TD vs D7 final sample


design.options<-c("design.timepoint",
                  "design.timepoint.diagnosis",
                  #"design.timepoint.diagnosis_paired",
                  "design.timepoint.diagnosis.vaccine"
                  # ,"design.timepoint.diagnosis.vaccine_paired"
)

cm.options<-c("cm.timepoint",
              "cm.timepoint.diagnosis",
              #  "cm.timepoint.diagnosis_paired",
              "cm.timepoint.diagnosis.vaccine"
              #  ,"cm.timepoint.diagnosis.vaccine_paired"
)



#### cm.timepoint #####

#.28/D28 vs Dremove, .remove, final samp, TD, D0
timepoint <- relevel(timepoint, "CT")


design.timepoint<- model.matrix(~0+Red.cell.count+timepoint)
# colnames(design.timepoint)
colnames(design.timepoint)<-gsub("timepoint","",colnames(design.timepoint))
colnames(design.timepoint)<-make.names(colnames(design.timepoint))
colnames(design.timepoint)


cm.timepoint <- makeContrasts(
  
  # D0 vs CT
  # CT vs D7
  # D0 vs D7
  
  
  

  
  D7vsCT=
    CT-D7,


  
  levels = design.timepoint
)


######  design.timepoint.diagnosis #####
timepoint.diagnosis<-as.factor(timepoint.diagnosis)
timepoint.diagnosis<-relevel(timepoint.diagnosis,"D7_had_covid")

design.timepoint.diagnosis<- model.matrix(~0+Red.cell.count+timepoint.diagnosis)


colnames(design.timepoint.diagnosis)<-gsub("timepoint.diagnosis","",colnames(design.timepoint.diagnosis))
colnames(design.timepoint.diagnosis)<-make.names(colnames(design.timepoint.diagnosis))
colnames(design.timepoint.diagnosis)



cm.timepoint.diagnosis <- makeContrasts(
  
  #cross section
  
  CT_COVID_vs_CT_not_COVID=
    CT_had_covid-CT_not_covid,
  
  # not doing D0_COVID_vs_D0_not_COVID as only one person had D0_not_covid, and I have renamed this timepoint anyway
  
  D7_COVID_vs_D7_not_COVID=
    D7_had_covid-D7_not_covid,
  
  # longitudinal
  
  # CT vs baseline
  

  
  # CT vs D7
  
  D7_all_vs_CT_COVID=
    CT_had_covid-(D7_not_covid+D7_had_covid)/2,
  
  D7_all_vs_CT_not_COVID=
    CT_not_covid-(D7_not_covid+D7_had_covid)/2,
  
  D7_had_COVID_vs_CT_COVID=
    CT_had_covid-D7_had_covid,
  
  D7_not_COVID_vs_CT_not_COVID=
    CT_not_covid-D7_not_covid,
  
  

  
  # CT -ve vs all non covid samples

  
  

  levels=design.timepoint.diagnosis
)

# 
# # # ##### design.timepoint.diagnosis_paired #####
# 
# timepoint.diagnosis<-as.factor(timepoint.diagnosis)
# timepoint.diagnosis<-relevel(timepoint.diagnosis,"D7_not_covid")
# 
# design.timepoint.diagnosis_paired<- model.matrix(~0+Red.cell.count+participant + timepoint.diagnosis)
# 
# 
# 
# colnames(design.timepoint.diagnosis_paired)<-gsub("timepoint.diagnosis","",colnames(design.timepoint.diagnosis_paired))
# colnames(design.timepoint.diagnosis_paired)<-make.names(colnames(design.timepoint.diagnosis_paired))
# colnames(design.timepoint.diagnosis_paired)
# 
# 
# 
# cm.timepoint.diagnosis_paired <- makeContrasts(
# 
#   #cross section
# 
#   CT_COVID_vs_CT_not_COVID_paired=
#     CT_had_covid-CT_not_covid,
# 
#   # not doing D0_COVID_vs_D0_not_COVID in paired analysis as only one person had D0_not_covid, and I have renamed this timepoint anyway
# 
#   D7_COVID_vs_D7_not_COVID_paired=
#     D7_had_covid-D7_not_covid,
# 
#   # longitudinal
# 
#   # CT vs baseline
# 
# #  D0_all_vs_CT_COVID # not doing at CT +ve didn't have baseline samples
# 
# 
#   # CT vs D7
# 
#   D7_all_vs_CT_COVID_paired=
#     CT_had_covid-(D7_not_covid+D7_had_covid)/2,
# 
#   D7_all_vs_CT_not_COVID_paired=
#     CT_not_covid-(D7_not_covid+D7_had_covid)/2,
# 
#   D7_CT_COVID_vs_CT_COVID_paired=
#     CT_had_covid-D7_had_covid,
# 
#   D7_not_COVID_vs_CT_not_COVID_paired=
#     CT_not_covid-D7_not_covid/2,
# 
# ## not doing any paired analysis which includes baseline as there wasn't really any D0s who went onto have covid
# 
#   # # CT +ve vs all non covid samples
#   # CT_COVID_vs_all_non_COVID_samples_paired=
#   #   CT_had_covid-
#   #   (D0_not_covid+D7_not_covid+CT_not_covid)/3,
#   #
#   # CT_COVID_vs_all_other_samples_incl.D7_had_covid_paired=
#   #   CT_had_covid-
#   #   (D0_not_covid+CT_not_covid+D7_not_covid+D7_had_covid)/4,
# 
# # # CT -ve vs all non covid samples
# # CT_not_COVID_vs_own_non_symptomatic_samples_paired=
# #   CT_had_covid-
# #   (D0_not_covid+D7_not_covid)/2,
# #
# # CT_not_COVID_vs_all_other_samples_paired=
# #   CT_had_covid-
# #   (D0_not_covid+CT_had_covid+D7_not_covid+D7_had_covid)/4,
# #
# #
# #   acute_covid_vs_baseline_paired=
# #     (CT_had_covid+D7_had_covid)/2-(D0_not_covid),
# 
#   levels=design.timepoint.diagnosis_paired
# )

##### design.timepoint.vaccine unpaired #####

timepoint.diagnosis.vaccine<-as.factor(timepoint.diagnosis.vaccine)
timepoint.diagnosis.vaccine<-relevel(timepoint.diagnosis.vaccine,"CT_not_covid")

design.timepoint.diagnosis.vaccine<- model.matrix(~0+Red.cell.count+timepoint.diagnosis.vaccine)


colnames(design.timepoint.diagnosis.vaccine)<-gsub("timepoint.diagnosis.vaccine","",colnames(design.timepoint.diagnosis.vaccine))
colnames(design.timepoint.diagnosis.vaccine)<-make.names(colnames(design.timepoint.diagnosis.vaccine))
colnames(design.timepoint.diagnosis.vaccine)


cm.timepoint.diagnosis.vaccine <- makeContrasts(
  
  # cross section
  
  # does covid induce changes at CT when you have no vaccine?
  
  CT.had.COVID.MEN_vs_CT.had.COVID.ChAdOx___does.vaccination.ameliorate.acute.covid=
    CT_had_covid_MenACWY-CT_had_covid_ChAdOx1,
  
  CT.had.COVID.MEN_vs_CT.any___do.cov.naive.people.show.DE.compared.with.all.other.CTs=
    CT_had_covid_MenACWY-(CT_had_covid_ChAdOx1+CT_not_covid)/2,
  
  CT.had.COVID.MEN_vs_CT.not.covid___do.cov.naive.people.show.DE.compared.with.all.other.CTs=
    CT_had_covid_MenACWY-CT_not_covid,
  
  CT.had.COVID.CHADOX_vs_CT.not.covid___do.cov.naive.people.show.DE.compared.with.all.other.CTs=
    CT_had_covid_ChAdOx1-CT_not_covid,
  
  # does covid induce changes compare with baseline when you have no vaccine?
 
  
  # does covid induce changes at CT compared with D7 when you have no vaccine?
  
  CT.had.COVID.MEN_vs_D7_had_covid_MEN=
    CT_had_covid_MenACWY-D7_had_covid_MenACWY,
  
  CT.had.COVID.MEN_vs_D7_had_covid_any_vacc=
    CT_had_covid_MenACWY-(D7_had_covid_MenACWY+D7_had_covid_ChAdOx1)/2,
  
  CT.had.COVID.MEN_vs_D7_not_COVID=
    CT_had_covid_MenACWY-(D7_not_covid),
  
  CT.had.COVID.MEN_vs_D7_all=
    CT_had_covid_MenACWY-(D7_not_covid+D7_had_covid_ChAdOx1+D7_had_covid_MenACWY)/3,
  
  CT.had.COVID.MEN_vs_D7_had_cov_chad_and_D7_not_covid=
    CT_had_covid_MenACWY-(D7_not_covid+D7_had_covid_ChAdOx1)/2,
  
  # does covid induce changes at CT compared with any other sample?
  


  
  levels=design.timepoint.diagnosis.vaccine
)


# ##### design.timepoint.vaccine paired #####
# 
# timepoint.diagnosis.vaccine<-as.factor(timepoint.diagnosis.vaccine)
# timepoint.diagnosis.vaccine<-relevel(timepoint.diagnosis.vaccine,"D7_not_covid")
# 
# design.timepoint.diagnosis.vaccine_paired<- model.matrix(~0+Red.cell.count+participant+timepoint.diagnosis.vaccine)
# 
# 
# colnames(design.timepoint.diagnosis.vaccine_paired)<-gsub("timepoint.diagnosis.vaccine","",colnames(design.timepoint.diagnosis.vaccine_paired))
# colnames(design.timepoint.diagnosis.vaccine_paired)<-make.names(colnames(design.timepoint.diagnosis.vaccine_paired))
# colnames(design.timepoint.diagnosis.vaccine_paired)
# 
# 
# cm.timepoint.diagnosis.vaccine_paired <- makeContrasts(
# 
#   # cross section
# 
#   # does covid induce changes at CT when you have no vaccine?
# 
#   CT.had.COVID.MEN_vs_CT.had.COVID.ChAdOx___does.vaccination.ameliorate.acute.covid=
#     -CT_had_covid_ChAdOx1,
# 
#   CT.had.COVID.MEN_vs_CT.any___do.cov.naive.people.show.DE.compared.with.all.other.CTs=
#     -(CT_had_covid_ChAdOx1+CT_not_covid)/2,
# 
#   CT.had.COVID.MEN_vs_CT.not.covid___do.cov.naive.people.show.DE.compared.with.all.other.CTs=
#     -CT_not_covid,
# 
#   # does covid induce changes compare with baseline when you have no vaccine?
# 
#   CT.had.COVID.MEN_vs_D0_any=
#     -D0_not_covid,
# 
#   # does covid induce changes at CT compared with D7 when you have no vaccine?
# 
#   # CT.had.COVID.MEN_vs_D7_had_covid_MEN=
#   #   -D7_had_covid_MenACWY,
#   #
#   # CT.had.COVID.MEN_vs_D7_had_covid_any_vacc=
#   #   -(D7_had_covid_MenACWY+D7_had_covid_ChAdOx1)/2,
# 
#   CT.had.COVID.MEN_vs_D7_not_COVID=
#     -(D7_not_covid),
# 
#   CT.had.COVID.MEN_vs_D7_all=
#     -(D7_not_covid+D7_had_covid_ChAdOx1+D7_had_covid_MenACWY)/2,
# 
#   CT.had.COVID.MEN_vs_D7_had_cov_chad_and_D7_not_covid=
#     -(D7_not_covid+D7_had_covid_ChAdOx1)/2,
# 
#   # does covid induce changes at CT compared with any other sample?
# 
#   CT.had.COVID.MEN_vs_D0_and_D7_all=
#     -(D7_not_covid+
#                             D7_had_covid_ChAdOx1+
#                             D7_had_covid_MenACWY+
#                             D0_not_covid
#     )/4,
# 
#   CT.had.COVID.MEN_vs_everything_except_any_CT_positive=
#     -(D7_not_covid+
#                             D7_had_covid_ChAdOx1+
#                             D7_had_covid_MenACWY+
#                             D0_not_covid+
#                             CT_not_covid
#     )/5,
# 
#   CT.had.COVID.MEN_vs_all_else=
#     -(D7_not_covid+
#                             D7_had_covid_ChAdOx1+
#                             D7_had_covid_MenACWY+
#                             D0_not_covid+
#                             CT_not_covid+
#                             CT_had_covid_ChAdOx1
#     )/6,
# 
#   levels=design.timepoint.diagnosis.vaccine_paired
# )
# 

## paired analyses are possible for CT cov -ve vs baseline
## CT.+ve.vaccine vs CT.+ve.vaccine
## theoretically possible for CT

# need to think about this a bit more. 
