
###############################################
##### PART III: Band combination: Indices and thresholding for flood mapping ##############

### NIR experiment with threshold to  map water/flooding

plot(subset(r_before,"NIR"))
plot(subset(r_after,"NIR"))

### Lower NIR often correlates to areas with high water fraction or inundated:
r_rec_NIR_before <- subset(r_before,"NIR") < 0.2
r_rec_NIR_after <- subset(r_after,"NIR") < 0.2

### THis is suggesting flooding!!!
plot(r_rec_NIR_before, main="Before RITA, NIR > 0.2")
plot(r_rec_NIR_after, main="After RITA, NIR > 0.2")

#Compare to actual flooding data
freq_fema_zones <- as.data.frame(freq(r_ref))
xtab_threshold <- crosstab(r_ref,r_rec_NIR_after,long=T)

## % overlap between the flooded area and values below 0.2 in NIR
(xtab_threshold[5,3]/freq_fema_zones[2,2])*100 #agreement with FEMA flooded area in %.
