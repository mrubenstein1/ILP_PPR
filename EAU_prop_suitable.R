### Establish Percent of Suitable Habitat for EAUs ####
  # Characterizing the percent of suitable habitat for each EAU
  # is a critical first step to establishing both benefit, cost, and transition estimates


#Load land cover data under RCP 4.5 and 8.5 scenarios
  # note: this uses the PPR BAU scenario
  # https://www.sciencebase.gov/catalog/item/5b683db2e4b006a11f75b06a
lu_4.5 <- "input_data/Land-use and vegetation change 2018/"
lu_8.5 <- 

#mask to EAU resolution
aoi_mask <- rast("input_data/wmd_raster_equal_area.tif")

file_names <- grep("gcam_ref_rcp85_2", list.files(LU_in), value = TRUE)

rcl_table <- matrix(c(19,1,28,1,29,1), 
                    byrow = TRUE, ncol = 2)

scf <- (res(aoi_mask)[1]^2) / (30^2)

lu_prop <- matrix(NA, 
                  nrow = length(values(aoi_mask, na.rm = TRUE)-1), 
                  ncol = length(file_names))

colnames(lu_prop) <- sprintf("rcp85_%s", c(2014, 2020, 2030, 2040, 2050,
                                           2060, 2070, 2080, 2090, 2100))

for(i in seq_along(file_names)){
  
  #load raster
  lu_r <- rast(paste0(LU_in, file_names[i]))
  
  #mask to AOI
  lu_r <- crop(lu_r, aoi_mask)
  
  #keep only chosen values
  lu_r <- classify(lu_r, rcl = rcl_table, others = NA)
  
  #aggregate using proportion
  lu_r <- resample(lu_r, aoi_mask, method = 'sum')
  lu_r <- mask(lu_r/scf, aoi_mask)
  
  #extract values / save?
  lu_prop[,i] <- values(lu_r, na.rm = TRUE)
  
  rm(lu_r)
  gc(ful = TRUE)
  
  write.csv(lu_prop, "Data/LU_modified/Lu_prop.csv", row.names = FALSE)
  
}

