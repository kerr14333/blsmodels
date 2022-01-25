# BLS report level model (microdata)

Created by Julie Gershunskaya. Programs used to generate simulated survey data and to run the Stan model. The program of most interest `step2_run_pathfinder.R` and the associated log from running the program `step2_pathfinder.log`.

* **/reportlevel**
    * **/data** 
        * `pop.rds` - Frame (population) we select our sample from.
        * `smpl.rds` - Subsample from `pop.rds` we use as our sample.
        * `stan_data.rds` - Inputs created for `step2\*.R` programs to run stan.
    * **/fit**
    * **/stan**
        * `y_countimp_ispbase_010522.stan` - Stan file containing our MRP model.
    * `step0_smpl.R` - Creates the simulated population (`pop.rds`) and the sample (`smpl.rds`). Does some other checks.
    * `step1_stan_data.R` - Creates the data for Stan from the population and sample files.
    * `step2_run_hmc.R` - Performs some last formatting on the data for stan then runs HMC. Saves model fit in
    * `step2_run_pathfinder.R` - Performs some last formatting on the data for stan and then runs pathfinder. Saves the model fit
    * `step2_run_pathfinder.log` - Log from R program run.
    * `step3_summary.R` - Creates summaries from ouputs.    
