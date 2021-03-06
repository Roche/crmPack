##Construct the DualResponsesDesign for simulations
##The design comprises the DLE and efficacy models, the escalation rule, starting data, 
##a cohort size and a starting dose
##Define your data set first using an empty data set 
## with dose levels from 25 to 300 with increments 25
data <- DataDual(doseGrid=seq(25,300,25),placebo=FALSE)
##First for the DLE model 
##The DLE model must be of 'ModelTox' (e.g 'LogisticIndepBeta') class 
DLEmodel <- LogisticIndepBeta(binDLE=c(1.05,1.8),
                              DLEweights=c(3,3),
                              DLEdose=c(25,300),
                              data=data)

##The efficacy model of 'ModelEff' (e.g 'Effloglog') class 
Effmodel<-Effloglog(Eff=c(1.223,2.513),Effdose=c(25,300),
                    nu=c(a=1,b=0.025),data=data,c=0)

##The escalation rule using the 'NextBestMaxGain' class
mynextbest<-NextBestMaxGain(DLEDuringTrialtarget=0.35,
                            DLEEndOfTrialtarget=0.3)


##The increments (see Increments class examples) 
## 200% allowable increase for dose below 300 and 200% increase for dose above 300
myIncrements<-IncrementsRelative(intervals=c(25,300),
                                 increments=c(2,2))
##cohort size of 3
mySize<-CohortSizeConst(size=3)
##Stop only when 36 subjects are treated
myStopping <- StoppingMinPatients(nPatients=36)
##Now specified the design with all the above information and starting with a dose of 25

design <- DualResponsesDesign(nextBest=mynextbest,
                              model=DLEmodel,
                              Effmodel=Effmodel,
                              stopping=myStopping,
                              increments=myIncrements,
                              cohortSize=mySize,
                              data=data,startingDose=25)