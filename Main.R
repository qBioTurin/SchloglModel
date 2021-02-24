library(epimod)
library(ggplot2)
library(readr)

# model_generation(net_fname = "./Net/ShlogelModel.PNPRO")
# 
# model_generation(net_fname = "./Net/ShlogelModel.PNPRO",
#                  functions_fname = "cpp/transitions.cpp")

system(paste('mv', 
             sprintf("ShlogelModel.*"),
             sprintf("NetFn/")) )


Exp.mngr=function(Gen,sim="SSA",teps=""){
  
  if(Gen){
    file_to_read<-"Input/Functions_list_FNgeneral.csv"
    solver <- "NetFn/ShlogelModel.solver"
    folder<-paste0("GeneralShlogel",sim,teps)
  }else{
    #file_to_read<-"Input/Functions_list_ModelAnalysis.csv"
    solver <- "Net/ShlogelModel.solver"
    folder<- paste0("Shlogel",sim,teps)
  }
  

  if(Gen){
    model_analysis(out_fname = "model_analysis",
                   solver_fname = solver,
                   parameters_fname = file_to_read,
                   solver_type = sim,
                   taueps = teps,
                   n_run = 1000,
                   f_time = 20, # weeks
                   s_time = 1
    )
  }else{
    model_analysis(out_fname = "model_analysis",
                   solver_fname = solver,
                   solver_type = sim,
                   taueps = teps,
                   n_run = 1000,
                   f_time = 20, # weeks
                   s_time = 1
    )
  }
  
  
  
  if(file.exists(folder)) { 
    system(paste('rm -rd ', sprintf(folder)) )
  }
  
  system(paste('mv', 
               sprintf("results_model_analysis"),
               sprintf(folder)) )
}

Exp.mngr(Gen=T,sim="SSA")
Exp.mngr(Gen=F,sim="SSA")

Exp.mngr(Gen=T,sim="TAUG",teps=0.01)



trace=as.data.frame(read.csv( "./ShlogelSSA/model_analysis-1.trace", sep = ""))
Gentrace=as.data.frame(read.csv( "./GeneralShlogelSSA/model_analysis-1.trace", sep = ""))

trace$type <- "No gen fun"
Gentrace$type <- "With gen fun"

TracesAll <-rbind(trace,Gentrace)
  
SummaryTrace<-function(trace){
  n_sim_tot<-table(trace$Time)
  n_sim <- n_sim_tot[1]
  time_delete<-as.numeric(names(n_sim_tot[n_sim_tot!=n_sim_tot[1]]))
  
  if(length(time_delete)!=0) trace = trace[which(trace$Time!=time_delete),]
  
  trace$ID <- rep(1:n_sim[1],each = length(unique(trace$Time)) )
  
  traceNew <-  data.frame(Time = rep(trace$Time,3),
                          Values = c( trace$B1, trace$X1, trace$B2),
                          place= rep(c("B1","X1","B2"),each = length(trace$Time)))
  
  
  Mean<-aggregate(trace[, 2:4], list(Time=trace$Time), mean)
  sd<- aggregate(trace[, 2:4], list(Time=trace$Time), sd)
  
  NewTrace<-data.frame(Time = rep(Mean$Time,3),
                       Mean = c( Mean$B1, Mean$X1, Mean$B2),
                       place= rep(c("B1","X1","B2"),each = length(Mean$Time)),
                       sd = c( sd$B1,sd$X1,sd$B2))
  
  return(list(Summary=NewTrace, trace = trace))
}

SummaryTrace(trace) -> traceNew
SummaryTrace(Gentrace) -> GentraceNew
traceNew$Summary$type <- traceNew$trace$type <- "No gen fun"
GentraceNew$Summary$type <- GentraceNew$trace$type <- "With gen fun"

AlltracesSummary <- rbind(traceNew$Summary, GentraceNew$Summary)
Alltraces <- rbind(traceNew$trace, GentraceNew$trace)

Alltraces$ID <- rep(1:(table(Alltraces$Time)[1]),
                    each = length(unique(Alltraces$Time)) )

ggplot(Alltraces,aes(x=Time)) +
  geom_line(aes(y=X1, group=ID, col=type, linetype=type))+
  facet_wrap(~type)


ggplot(AlltracesSummary,aes(x=Time)) +
  geom_line(aes(y=Mean, col=type, linetype=type))+
  geom_point(aes(y=Mean, col=type, shape=type))
  facet_wrap(~place,scales = "free_y")
