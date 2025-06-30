library(here)
library(tidyverse)
library(tools)
library(ncdf4)
library(FLAREr)
working_directory <- here()
setwd(working_directory)
source("get_glm_nc_var.R")
source("get_glm_nc_cont.R")

system2("rm", args = "glm3.nml")
system2("rm", args = "glm3_full.nml")
system2("rm", args = "output.nc")
system2("rm", args = "txt/*")
system2("rm", args = "*.txt")
system2("rm", args = "*.nc")
system2("rm", args = "nc/*")
system2("rm", args = "nml/*")
system2("rm", args = "Rplots.pdf")

states <- c("temp", "salt")
#states <- c("temp", "salt", "OXY_oxy")
#states <- c("temp","salt","OXY_oxy","CAR_dic","CAR_ch4","SIL_rsi","NIT_amm","NIT_nit","PHS_frp",
#            "OGM_doc","OGM_docr","OGM_poc","OGM_don","OGM_donr","OGM_pon","OGM_dop","OGM_dopr","OGM_pop",
#            "PHY_cyano","PHY_green","PHY_diatom")
round_num <- 20
file.copy("glm3_initial_spark.nml", "glm3.nml", overwrite = TRUE)
file.copy("glm3_initial_spark.nml", "glm3_full.nml", overwrite = TRUE)
full_time <- seq(as_datetime("1997-04-16 00:00:00"), as_datetime("1997-08-20 00:00:00"), "1 day")
#full_time <- seq(as_datetime("2021-01-01 00:00:00"), as_datetime("2021-06-15 00:00:00"), "1 day")
full_time_string <- strftime(full_time,
                             format="%Y-%m-%d %H:%M",tz = "UTC")
glm_list_index <- 1
update_glm_nml_list <- list()
update_glm_nml_names <- list()
update_glm_nml_list[[glm_list_index]] <- full_time_string[1]
update_glm_nml_names[glm_list_index] <- "start"
glm_list_index <- 1
update_glm_nml_list[[glm_list_index]] <- full_time_string[2]
update_glm_nml_names[glm_list_index] <- "stop"
glm_list_index <- 2
wq_states <- length(states) - 2
update_glm_nml_list[[glm_list_index]] <- wq_states
update_glm_nml_names[glm_list_index] <- "num_wq_vars"
glm_list_index <- 3
wq_states <- length(states) - 2
update_glm_nml_list[[glm_list_index]] <- wq_states
update_glm_nml_names[glm_list_index] <- "num_wq_vars"
glm_list_index <- 3
if(wq_states > 0){
  update_glm_nml_list[[glm_list_index]] <- states[3:length(states)]
  update_glm_nml_names[glm_list_index] <- "wq_names"
  glm_list_index <- 4
}
FLAREr:::update_nml(var_list = update_glm_nml_list,
                    var_name_list = update_glm_nml_names,
                    working_directory = working_directory,
                    nml = "glm3.nml")
glm_list_index <- 1
update_glm_nml_list <- list()
update_glm_nml_names <- list()
update_glm_nml_list[[glm_list_index]] <- full_time_string[1]
update_glm_nml_names[glm_list_index] <- "start"
glm_list_index <- 1
update_glm_nml_list[[glm_list_index]] <- full_time_string[length(full_time_string)]
update_glm_nml_names[glm_list_index] <- "stop"
glm_list_index <- 2
wq_states <- length(states) - 2
update_glm_nml_list[[glm_list_index]] <- wq_states
update_glm_nml_names[glm_list_index] <- "num_wq_vars"
glm_list_index <- 3
wq_states <- length(states) - 2
update_glm_nml_list[[glm_list_index]] <- wq_states
update_glm_nml_names[glm_list_index] <- "num_wq_vars"
glm_list_index <- 3
if(wq_states > 0){
  update_glm_nml_list[[glm_list_index]] <- states[3:length(states)]
  update_glm_nml_names[glm_list_index] <- "wq_names"
  glm_list_index <- 4
}
FLAREr:::update_nml(var_list = update_glm_nml_list,
                    var_name_list = update_glm_nml_names,
                    working_directory = working_directory,
                    nml = "glm3_full.nml")
nsteps <- length(full_time_string) - 1
one_step_output <- array(NA, dim = c(nsteps, 500, length(states)))
lake_depth <- array(NA, dim = c(nsteps))
z <- array(NA, dim = c(nsteps, 500))
white_ice_thickness <- array(NA, dim = c(nsteps))
blue_ice_thickness <- array(NA, dim = c(nsteps))
avg_surf_temp <- array(NA, dim = c(nsteps))
mixing_vars <-  array(NA, dim = c(nsteps, 17))
chla <- array(NA, dim = c(nsteps, 500))
mixer_count <- array(NA, dim = c(nsteps))
for(i in 1:nsteps){
  if(i > 1){
    glm_list_index <- 1
    update_glm_nml_list <- list()
    update_glm_nml_names <- list()
    curr_start <- as.character(full_time_string[i])
    curr_stop <- as.character(full_time_string[i+1])
    update_glm_nml_list[[glm_list_index]] <- curr_start
    update_glm_nml_names[glm_list_index] <- "start"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- curr_stop
    update_glm_nml_names[glm_list_index] <- "stop"
    glm_list_index <- glm_list_index + 1
    native_depth_index <- which(!is.na(z[i-1, ]))
    the_temps_glm <- rev(one_step_output[i-1, native_depth_index, 1])
    the_temps_glm <- one_step_output[i-1, native_depth_index, 1]
    update_glm_nml_list[[glm_list_index]] <- round(the_temps_glm, round_num)
    update_glm_nml_names[glm_list_index] <- "the_temps"
    glm_list_index <- glm_list_index + 1
    the_sals_glm <- rev(one_step_output[i-1, native_depth_index, 2])
    the_sals_glm <- one_step_output[i-1, native_depth_index, 2]
    update_glm_nml_list[[glm_list_index]] <- round(the_sals_glm, round_num)
    update_glm_nml_names[glm_list_index] <- "the_sals"
    glm_list_index <- glm_list_index + 1
    #the_depths <- rev(lake_depth[i-1] - z[i-1, native_depth_index])
    #update_glm_nml_list[[glm_list_index]] <- round(the_depths, round_num)
    #update_glm_nml_names[glm_list_index] <- "the_depths"
    #glm_list_index <- glm_list_index + 1
    #update_glm_nml_list[[glm_list_index]] <- length(the_depths)
    #update_glm_nml_names[glm_list_index] <- "num_depths"
    #glm_list_index <- glm_list_index + 1
    the_heights <- z[i-1, native_depth_index]
    update_glm_nml_list[[glm_list_index]] <- round(the_heights, round_num)
    update_glm_nml_names[glm_list_index] <- "the_heights"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- length(the_heights)
    update_glm_nml_names[glm_list_index] <- "num_heights"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- round(lake_depth[i-1], round_num)
    update_glm_nml_names[glm_list_index] <- "lake_depth"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- 0.0
    update_glm_nml_names[glm_list_index] <- "snow_thickness"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- round(white_ice_thickness[i-1], round_num)
    update_glm_nml_names[glm_list_index] <- "white_ice_thickness"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- round(blue_ice_thickness[i-1], round_num)
    update_glm_nml_names[glm_list_index] <- "blue_ice_thickness"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- round(avg_surf_temp[i-1], round_num)
    update_glm_nml_names[glm_list_index] <- "avg_surf_temp"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- mixing_vars[i-1, ]
    update_glm_nml_names[glm_list_index] <- "restart_variables"
    glm_list_index <- glm_list_index + 1
    update_glm_nml_list[[glm_list_index]] <- mixer_count[i-1]
    update_glm_nml_names[glm_list_index] <- "restart_mixer_count"
    glm_list_index <- glm_list_index + 1
    if((length(states)-2) > 0){
      wq_init_vals <- c()
      start_index <- 2
      for(wq in 1:(length(states)-2)){
        wq_tmp <- one_step_output[i-1, native_depth_index, start_index + wq]
        wq_init_vals <- c(wq_init_vals, wq_tmp)
      }
      update_glm_nml_list[[glm_list_index]] <- round(wq_init_vals, round_num)
      update_glm_nml_names[glm_list_index] <- "wq_init_vals"
    }
    FLAREr:::update_nml(var_list = update_glm_nml_list,
                        var_name_list = update_glm_nml_names,
                        working_directory = working_directory,
                        nml = "glm3.nml")
  }
  #GLM3r::run_glm(sim_folder = working_directory, nml_file = "glm3.nml", verbose = FALSE)
#   system2("/workspaces/focal_py39/AED_Tools/binaries/ubuntu/20.04/glm_latest/glm", stdout = paste0("txt/output", i, ".txt"))
  system2("/workspaces/focal_py39/2025/AED_Tools_Private/binaries/ubuntu/20.04/glm_latest/glm", stdout = paste0("txt/output", i, ".txt"))
  GLM_temp_wq_out <- get_glm_nc_var(ncFile = "/output.nc",
                                    working_dir = working_directory,
                                    z_out = NA,
                                    vars_depth = states,
                                    vars_no_depth = NULL,
                                    diagnostic_vars = NULL)
  file.copy("output.nc", file.path("nc", paste0("output_", i, ".nc")), overwrite = TRUE)
  system2("rm", args = "output.nc")
  mixing_vars[i, ] <- GLM_temp_wq_out$mixing_vars
  mixer_count[i] <- GLM_temp_wq_out$mixer_count
  avg_surf_temp[i] <- GLM_temp_wq_out$avg_surf_temp
  white_ice_thickness[i] <- GLM_temp_wq_out$snow_wice_bice[2]
  blue_ice_thickness[i] <- GLM_temp_wq_out$snow_wice_bice[3]
  lake_depth[i] <- GLM_temp_wq_out$lake_depth
  num_glm_heights <- length(GLM_temp_wq_out$heights)
  z[i, 1:num_glm_heights] <- GLM_temp_wq_out$heights
  if("PHY_cyano" %in% states){
    chla[i, 1:num_glm_heights]<- GLM_temp_wq_out$chla
  }
  for(s in 1:length(states)){
    one_step_output[i,1:num_glm_heights , s]<- GLM_temp_wq_out$output[ ,s]
  }
}
##### DO CONTINUIOUS RUN ####
# system2("/workspaces/focal_py39/AED_Tools/binaries/ubuntu/20.04/glm_latest/glm", args = "--nml glm3_full.nml", stdout = "output.txt")
system2("/workspaces/focal_py39/2025/AED_Tools_Private/binaries/ubuntu/20.04/glm_latest/glm", args = "--nml glm3_full.nml", stdout = "output.txt")
#GLM3r::run_glm(sim_folder = working_directory, nml_file = "glm3_full.nml")
full_output <- get_glm_nc_cont(ncFile = "/output.nc",
                               working_dir = working_directory,
                               z_out = NA,
                               vars_depth = states,
                               vars_no_depth = NULL,
                               diagnostic_vars = NULL)
# Layers
t <- length(full_time)-1
layer_full <- full_output$heights[full_output$heights[,t] < 9.969210e+35, t]
layer_step <- z[t,!is.na(z[t, ])]
layer_full
layer_step
length(full_output$heights[full_output$heights[,t] < 9.969210e+35, t])
length(z[t,!is.na(z[t, ])])
#Mixing vars
#full_output$mixing_vars - mixing_vars[t, ]
#Temp
temp_full <- full_output$output[full_output$heights[,t] < 9.969210e+35, t, 1]
temp_step <- one_step_output[t, !is.na(z[t, ]), 1]
combined <- bind_rows(tibble(height = layer_full,
                             temp = temp_full,
                             type = "full"),
                      tibble(height = layer_step,
                             temp = temp_step,
                             type = "step"))
ggplot(combined, aes(x = temp, y = height, color = type)) +
  geom_point(aes(size = type)) +
  scale_size_manual(values = c("full" = 4, "step" = 2))
oxygen_full <- full_output$output[full_output$heights[,t] < 9.969210e+35, t, 3]


