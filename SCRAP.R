#Scrap code for testing

#Step 1----
#Process detections and stations
proc.det <- proc_dets(det = acoustic$detections, sta = acoustic$stations)

#Plot station detection history
plot_sta_history(proc.det)

#Map station detections
map_dets(proc.det)

#Map station detections with options
#1
map_dets(proc.det, base.layers = list(acoustic$land),
         base.cols = list("wheat"), base.borders = list("seagreen"))
#2
map_dets(proc_det = proc.det, sta.crs = 4326,
         base.layers = list(acoustic$study_area, acoustic$land),
         base.cols = list("gray30", "wheat"), base.borders = list(NA, "seagreen"))
#3
map_dets(proc_det = proc.det, sta.crs = 4326,
         base.layers = list(acoustic$study_area, acoustic$land),
         base.cols = list("gray30", "wheat"), base.borders = list(NA, "seagreen"),
         xlim=c(-64.628, -64.612), ylim=c(17.770, 17.795),
         leg.pos = "topleft", sta.col = "blue", sta.bg = "orange")

#Step 2----

#Step 3----
#1 hour COAs (default)
coas.60 <- coa_locs(proc.det)

#20 minute COAs
coas.20 <- coa_locs(proc.det, Delta_t = "20 min")

#30 minute COAS, harmonic mean
coas.30.harm <- coa_locs(proc.det, Delta_t = "30 min", mean_type="harmonic")

#Map COAs
map_coas(proc.det, coas.60)

#Map COAs with options
#1
map_coas(proc.det, coas = coas.60, base.layers = list(acoustic$land),
         base.cols = list("wheat"), base.borders = list("seagreen"))

#2
map_coas(proc.det, coas = coas.60, base.layers = list(acoustic$land),
         base.cols = list("wheat"), base.borders = list("seagreen"),
         det.leg.pos = "topright", coa.leg.pos = "topleft")

#Step 4----
#Initialize raster

#study_area
r.sa <- init_raster(study_area = acoustic$study_area, res = 0.005)

#Rasterize COAs
coa.r <- rasterize_dets(det_locs = coas.60, r = r.sa)

#Plot
plot(coa.r)
plot(acoustic$study_area, add=TRUE)
plot(acoustic$land, add=TRUE, col="wheat")
