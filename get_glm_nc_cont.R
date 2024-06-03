get_glm_nc_cont <- function(ncFile, working_dir, z_out, vars_depth, vars_no_depth,
                           diagnostic_vars)
{
  glm_nc <- ncdf4::nc_open(paste0(working_dir, ncFile))
  glm_vars <- names(glm_nc$var)
  tallest_layer <- ncdf4::ncvar_get(glm_nc, "NS")
  heights <- matrix(ncdf4::ncvar_get(glm_nc, "z"), ncol = length(tallest_layer))
  snow <- matrix(ncdf4::ncvar_get(glm_nc, "snow_thickness"), ncol = length(tallest_layer))
  ice_white <- matrix(ncdf4::ncvar_get(glm_nc, "white_ice_thickness"), ncol = length(tallest_layer))
  ice_blue <- matrix(ncdf4::ncvar_get(glm_nc, "blue_ice_thickness"), ncol = length(tallest_layer))
  avg_surf_temp <- matrix(ncdf4::ncvar_get(glm_nc, "avg_surf_temp"), ncol = length(tallest_layer))

  
  hours_since <- as.numeric(ncvar_get(glm_nc, "time"))
  time_units <- ncatt_get(glm_nc, "time", "units")$value
  
  start_time <- str_sub(time_units, start = 13, end = 31)
  times <- as_datetime(start_time) + hours(hours_since) 

  output <- array(NA, dim = c(500, length(times), length(vars_depth)))
  for (v in 1:length(vars_depth)) {
    output[, , v] <- matrix(ncdf4::ncvar_get(glm_nc, vars_depth[v]), ncol = length(tallest_layer))
  }
  
  mixing_vars <- ncdf4::ncvar_get(glm_nc, "restart_variables")
  salt <- matrix(ncdf4::ncvar_get(glm_nc, "salt"), ncol = length(tallest_layer))
  if("PHY_cyano" %in% states){
    chla <- matrix(ncdf4::ncvar_get(glm_nc, "PHY_tchla"), ncol = length(tallest_layer))
  }
  
  ncdf4::nc_close(glm_nc)
  
  return(list(time = times, output = output, heights = heights, snow = snow, ice_white = ice_white, ice_blue = ice_blue,
              avg_surf_temp = avg_surf_temp, mixing_vars = mixing_vars,
              salt = salt, chla = chla))
}