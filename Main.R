
library(epimod)

# downloadContainers()

start_time <- Sys.time()
model.generation(net_fname = "./Net/Schlogl_general.PNPRO",
                 transitions_fname = "./Cpp/transitions.cpp")

## Notice that model.generation() might take as input parameter a C++ 
## file defining the functions characterizing the behavior of general transitions

end_time <- Sys.time()-start_time

source("Rfunction/ModelAnalysisPlot.R")

### Model Analysis

##############################
## Deterministic setting
##############################

## Two scenarios:

## 1) permutations are approximated
## 2) permutations are no approximated

model.generation(net_fname = "./Net/Schlogl_general.PNPRO",
                 transitions_fname = "./Cpp/transitions.cpp")

model.analysis(solver_fname = "./Schlogl_general.solver",
               parameters_fname = "./Input/Functions_list_ModelAnalysis.csv",
               f_time = 100, # days
               s_time = 1
)

source("Rfunction/ModelAnalysisPlot.R")

AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./results_model.analysis")

p1 <- AnalysisPlot$list.plX1$plX1.1

ggarrange(p1, p2, ncol = 2, nrow = 1)

# approximated case

model.generation(net_fname = "./Net/Schlogl_general.PNPRO",
                 transitions_fname = "./Cpp/transitions_apprx.cpp")

model.analysis(solver_fname = "./Schlogl_general.solver",
               parameters_fname = "./Input/Functions_list_ModelAnalysis.csv",
               f_time = 100, # days
               s_time = 1
)

source("Rfunction/ModelAnalysisPlot.R")
AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./results_model.analysis")
p2 <- AnalysisPlot$list.plX1$plX1.1

library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1)

## Stochastic setting

model.analysis(solver_fname = "./Schlogl_general.solver",
               parameters_fname = "./Input/Functions_list_ModelAnalysis.csv",
               functions_fname = "./Rfunction/Functions.R",
               solver_type = "SSA",
               n_run = 100,
               parallel_processors = 2,
               f_time = 100, # days
               s_time = 1
)

AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./results_model.analysis",
                                 Stoch = T)
AnalysisPlot$list.plX1

## parameters configurations

model.analysis(solver_fname = "./Schlogl_general.solver",
               parameters_fname = "./Input/Functions_list2.csv",
               functions_fname = "./Rfunction/Functions.R",
               solver_type = "SSA",
               n_run = 100,
               n_config = 8,
               parallel_processors = 2,
               f_time = 100, # days
               s_time = 1
)

AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./results_model.analysis",
                                 Stoch = T)

### Sensitivity analysis

##########################################################
## Simple version where only the transition rates vary. ##
##########################################################

start_time <- Sys.time()
sensitivity<-sensitivity_analysis(n_config = 100,
                                  parameters_fname = "Input/Functions_list.csv", 
                                  solver_fname = "Schlogl_general.solver",
                                  reference_data = "Input/reference_data.csv",
                                  distance_measure_fname = "Rfunction/msqd.R" ,
                                  target_value_fname = "Rfunction/Target.R" ,
                                  f_time = 100, # days
                                  s_time = 1, # days      
                                  parallel_processors = 2
)

end_time <- Sys.time()-start_time

## Let draw the trajectories
source("./Rfunction/SensitivityPlot.R")
plX1

## Version where only the PRCC is calculated
# sensitivity<-sensitivity_analysis(n_config = 100,
#                                   parameters_fname = "Input/Functions_list.csv", 
#                                   functions_fname = "Rfunction/Functions.R",
#                                   solver_fname = "Schlogl_general.solver",
#                                   target_value_fname = "Rfunction/Target.R" ,
#                                   parallel_processors = 1,
#                                   f_time = 100 # days
#                                   s_time = 1 # days
#                                   )

## Version where only the ranking is calculated
# sensitivity<-sensitivity_analysis(n_config = 100,
#                                   parameters_fname = "Input/Functions_list.csv", 
#                                   functions_fname = "Rfunction/Functions.R",
#                                   solver_fname = "Schlogl_general.solver",
#                                   reference_data = "Input/reference_data.csv",
#                                   distance_measure_fname = "Rfunction/msqd.R" ,
#                                   parallel_processors = 1,
#                                   f_time = 100 # days
#                                   s_time = 1 # days
#                                   )

## Complete and more complex version where all the parameters for calculating
## the PRCC and the ranking are considered, and the initial conditions vary too.

start_time <- Sys.time()

sensitivity<-sensitivity_analysis(n_config = 100,
                                   parameters_fname = "Input/Functions_list2.csv", 
                                   functions_fname = "Rfunction/Functions.R",
                                   solver_fname = "Schlogl_general.solver",
                                   reference_data = "Input/reference_data.csv",
                                   distance_measure_fname = "Rfunction/msqd.R" ,
                                   target_value_fname = "Rfunction/Target.R" ,
                                   parallel_processors = 2,
                                   f_time = 30, # days
                                   s_time = 1 # days
                                   )

end_time <- Sys.time() - start_time

source("./Rfunction/SensitivityPlot.R")

plX1

### Calibration analysis

start_time <- Sys.time()

model.calibration(parameters_fname = "Input/Functions_list_Calibration.csv",
                  functions_fname = "Rfunction/FunctionCalibration.R",
                  solver_fname = "Schlogl_general.solver",
                  reference_data = "Input/reference_data.csv",
                  distance_measure_fname = "Rfunction/msqd.R" ,
                  f_time = 30, # days
                  s_time = 1, # days
                  # Vectors to control the optimization
                  ini_v = c(247, 0.06, 0.0005),
                  ub_v = c(249, 0.08, 0.0002),
                  lb_v = c(248, 0.028, 0.00009),
                  max.time = 1
)

end_time <- Sys.time()-start_time

source("Rfunction/CalibrationPlot.R")

plX1

################################
## from Rmd to pdf generation ##
################################

rmarkdown::render("ReadME.Rmd", "pdf_document")
