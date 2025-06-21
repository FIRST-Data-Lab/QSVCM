# QSVCM
Estimation and Inference of Quantile Spatially Varying Coefficient Models Over Complicated Domains\
(1) R code of the QBPST method (CodeData.zip) contains the following:
 - Source codes(fit.qsvcm.R, QRWboot.R, rhotau.R, stoper.R, cv.pred.qsvcm.R, QmixNorm.R, CQRPI.R). See details for description in each source function.
 - In Application folder, the main applications for mortality (Main_Application_Mortality.R) and PM2.5 (Main_Application_PM25.R).
 - In Application folder, mortality data (Data_mortality.csv) and PM2.5 data (Data_PM25.csv) 
 - In Application folder, triangulations (Tr_usa.csv, V_usa.csv).
 - In Application folder, population grid points for USA map (pop_location_d025.csv)
 - In Simulation Examples folder, Simulation Example 2 M1 (Main_Sim_M1.R) and M2 (Main_Sim_M2.R).
 - In Simulation Examples folder, simulation data (eg1_pop_mnorm.csv, eg1_pop_t2.csv, eg2_pop_mnorm.csv, eg2_pop_t2.csv) and triangulations (Tr_1.csv, Tr_2.csv, V_1.csv, V_2.csv).
 

(2) Details on data:
- The data file "Data_mortality.csv" contains 3,104 counties across 48 states (Alaska and Hawaii were excluded) and the District of Columbia in the United States with the following variables:
 > "Affluence": Social affluence is a metric that gauges economically privileged regions and takes into account the following factors: (i) the percentage of households with income over $75,000, (ii) the percentage of adults who have obtained a bachelor's degree or higher, (iii) the percentage of employed individuals in management, professional, and related occupations, and (iv) the median value of owner-occupied housing units.
 > "Disadvantage": The concentrated disadvantage which is a metric that assesses economic disadvantage and includes the following factors: (i) the percentage of households with public assistance income, (ii) the percentage of households with a female householder and no husband present, and (iii) the civilian labor force unemployment rate.
 > "Violent": Violent crime rate per 1,000 population. 
 > "UrbanRate": Urban rate.
 > "Latitude": Latitude of the county.
 > "Longtitude": Longtitude of the county.
 > "avemort": The response variable is the average age-standardized mortality rates per 100 population based on county level over the period of 1998-2002. 

 - The data file "Data_PM25.scv" contains 1,247 air quality monitoring stations throughout the United States with the following variables:
 > "winterpm": the average daily PM2.5 concentrations during winter for 2011.
 > "Latitude": Latitude of the county.
 > "Longtitude": Longtitude of the county.
 > "prec": daily total gridded precipitation.
 > "wind": surface wind speed.
 > "tmin": daily minimum air temperature.
 > "tmax": daily maximum air temperature.
 > "rhum": relative humidity.
 > "tcdc": total column cloud cover. 

- The data file "Tr_mort.csv" contains information on 459 triangles used in the application.

- The data file "V_mort.csv" contains information on 283 vertices used in the application.

- The data file "pop_location_d025.csv" corresponds to 22,4087 population grid locations to evaluate the prediction of the method at those in the application. 

- The data files "eq1_pop_mnorm.csv", "eq1_pop_t2.csv", "eg2_pop_mnorm.csv", "eg2_pop_t2.csv" are data files in Simulation Example 2 M1 and M2 for mixture normal (mnorm) and t2 distribution. Here, "eq1" and "eq2" correspond to M1 and M2, respectively.
  > "y": response value.
  > "m1", "m2" are the spatially varying functions "beta_0" and "beta_1", respectively, in the manuscript.
  > "x1", "x2" are 1 and explanatory variables, respectively, in the manuscript.
  > "u", "v" are spatial locations. 
  > "m3", "m4" are the spatially varying functions "zeta_0" and "zeta_1", respectively, in the manuscript. They only appear in the M2 model.

- The data files "Tr_1.csv", "Tr_2.csv" are triangulations used in Simulation Example 2 M1 and M2. 
- The data files "V_1.csv", "V_2.csv" are vertices used in Simulation Example 2 M1 and M2.


(3) Implementations:
(SIMULATION) 
  - 1. Open either "Simulation_Examples/Main_Sim_M1.R" or "Simulation_Examples/Main_Sim_M2.R" for Simulation Example2 M1 or M2, respectively.
  - 2. If necessary, install the following library: `mgcv', `MGLM', `BPST', `quantreg'. For installing `BPST' package, use the commands as follows.
	install.packages("devtools",repos="http://cran.r-project.org")
	library(devtools)
	install_github("funstatpackages/BPST")
  - 3. Run the code. As in our manuscript, other scenarios (E.g., n = 1000, 2000; tau = 0.50, 0.75; dist = "MN", "t2") can be chosen. 
  - 4. By default, the PCI using wild bootstrap is not implemented due to computational time. If one wants to implement it, select 1 for "run.bootstrap". 
  - 5. By default, the prediction interval for the response at the new location using the CQR (Conformalized quantile regression) is not implemented due to computational time. If one wants to implement it, select 1 for "run.CQRPI". 


(APPLICATION)
  - 1. Open "Application/Main_Application.R" for the application.
  - 2. If necessary, install the following library: `mgcv', `MGLM', `BPST', `quantreg'. For installing `BPST' package, use the commands as follows.
	install.packages("devtools",repos="http://cran.r-project.org")
	library(devtools)
	install_github("funstatpackages/BPST")
  - 3. Run the code. As in our manuscript, other scenarios (E.g., tau = 0.50, 0.75) can be chosen. 
  - 4. By default, the PCI using wild bootstrap is not implemented due to computational time. If one wants to implement it, select 1 for "run.bootstrap". 
  - 5. By default, the prediction interval for the response at the new location using the CQR is not implemented due to computational time. If one wants to implement it, select 1 for "run.CQRPI". 
