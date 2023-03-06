## Script to examine & format input layers for spatial prioritisation exercise
## Note: tifs are currently generated from earth engine to gee_exports, and I 
## have manually downloaded these to the outputs directory (at a later date 
## consider running these from R via reticulate).

## housekeeping ----
library(stars); library(sf); library(dplyr); library(data.table)

## read rasters ----
pct_mask <- read_stars("outputs/proportion_unmasked_177m.tif")
fc <- read_stars("outputs/pct_contig_forest_177m_res2.tif")
tc <- read_stars("outputs/pct_treecover_in_pasture_177m_res.tif")
#fc_2021 <- read_stars("outputs/pct_contig_forest_2021_177m_res.tif")
#tc_2021 <- read_stars("outputs/pct_treecover_in_pasture_2021_177m_res.tif")
ele <- read_stars("outputs/elev_ALOS_177m_res.tif")

# format: set pixels that are highly masked (i.e either exist below 1000 m or 
# are dominated by urban or water layers)
lc <- c(fc, pct_mask, tc, ele) %>%
    rename(fc = 1, pct_unmask = 2, tc = 3, ele = 4) %>%
    mutate(id_cell = 1:n(), 
           fc = ifelse(pct_unmask <.8, NA, fc), 
           tc = ifelse(pct_unmask <.8, NA, tc), 
           ele = ifelse(pct_unmask <.8, NA, ele), 
           id_cell = ifelse(is.na(ele), NA, id_cell)) %>%
    select(-pct_unmask)

rm(pct_mask, fc, tc, ele)

## examine hists ----
# forest cover
cuts <- cut(lc[["fc"]], seq(0, 1, by=.05), include.lowest=T)
round(table(cuts)/sum(!is.na(cuts)) * 100, 2)

png("figures/prop_cell_classified_as_forest.png")
hist(lc[["fc"]], 
     xlab = "Proportion of cell classified as forest", 
     main = "")
dev.off()

# proportion woody veg in pasture
cuts <- cut(lc[["tc"]], seq(0, 1, by=.05), include.lowest=T)
round(table(cuts)/sum(!is.na(cuts)) * 100, 2)

png("figures/prop_woody_veg_in_pasture.png")
hist(lc[["tc"]], 
     xlab = "Proportion of cell classified as forest", 
     main = "")
dev.off()

## convert to data.table ----
lc_df <- st_as_sf(lc, as_points = TRUE, na.rm = TRUE)

# convert to DT and format data
lc_dt <- as.data.table(lc_df)
lc_dt[, geometry := NULL]
lc_dt[, `:=`(tc = as.integer(round(tc, 0)),
             fc = as.integer(round(fc*100, 0)), 
             ele = as.integer(ele))]

setcolorder(lc_dt, c("id_cell", "ele", "tc", "fc"))

## save outputs ----
lc_dt <- saveRDS(lc_dt, "outputs/input_layers_dt.rds")
