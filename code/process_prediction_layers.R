library(stars); library(data.table); library(dplyr); library(sf)

## layers created in generate spatial layers script. 
elevation <- read_stars("data/ele_1km.tif")
pct_urban <- read_stars("data/pct_urban.tif")
pct_waterbody <- read_stars("data/pct_water.tif")

# note: these latter two are the percentages of forest in non-urban, 
# non-waterbody area, and the % treecover in non-urban, non-waterbody, non-forest
# areas respectively. 
pct_treecover <- read_stars("data/pct_treecover_in_pasture.tif")
pct_forest <- read_stars("data/pct_contig_forest.tif")
# cells that contain 0 forest are masked out. Set these to 0
# pct_forest[pct_urban[] + pct_waterbody[] == 1] <- 0

# convert all layers into a single indexed data table 
ele_sf <- elevation %>% 
    mutate(index = 1:n()) %>%
    st_as_sf(., as_points = T) %>%
    # filter(!is.na(ele_1km.tif)) %>%
    mutate(ele_ALOS = round(ele_1km.tif, 0)) %>%
    select(-ele_1km.tif) %>%
    filter(!is.na(ele_ALOS))

ele_dt <- as.data.table(ele_sf)
ele_dt[,geometry := NULL]

ele_dt[, `:=`(pct_urban = round(pct_urban[[1]][ele_dt$index], 2), 
              pct_water = round(pct_waterbody[[1]][ele_dt$index], 2), 
              pct_forest = round(pct_forest[[1]][ele_dt$index], 2), 
              avg_treecover = round(pct_treecover[[1]][ele_dt$index]/100, 2)
              )]

ele_dt[, pct_unmasked := 1 - (pct_urban + pct_water)]
ele_dt[is.na(avg_treecover), avg_treecover := 0]
ele_dt[is.na(pct_forest), pct_forest := 0]

saveRDS(ele_dt, "outputs/generic_input_layers.rds")

## Get range information ----
# Range information is encoded in the tdist raster, with values >0.5 out of range
# and <= 0.5 in range. The raster is also in a different projection. We need to
# convert tdist to represent in/out of range in the same projection and dimensions
# as the other rasters. 

# get species list
sp_list <- readRDS("../sparing_sharing_elevation/cu_cov.rds")$species %>%
    unique

# get tdist raster
predict_info <- readRDS("../Colombia/colombiaBeta/outputs/tdist_stars.rds") %>%
    filter(species %in% sp_list)

# convert to in/out of range
predict_info <- predict_info %>%
    mutate(tdist = tdist <= 0.5) %>%
    rename(in_range = tdist)

# crop to same region
bbox_EC <- st_as_sfc(st_bbox(elevation)) 
predict_info <- st_crop(predict_info, st_transform(bbox_EC, st_crs(predict_info)))

# iterate across species, convert to polygon before recasting back to raster 
# with same projection & dimensions and then extracting in/out of range info
grd <- elevation
grd$ele_1km.tif <- NA

catch <- vector("list", length(sp_list))
names(catch) <- sp_list
for(i in seq_along(sp_list)) {
    print(i)
    range_i <- predict_info[,,,i] %>%
        select(in_range) %>%
        st_as_sf(., as_points=FALSE, merge=T) %>%
        st_set_crs(., st_crs(predict_info)) %>%
        st_transform(., st_crs(elevation))

    rast_i <- st_rasterize(range_i, template=grd)
    
    catch[[sp_list[i]]] <- data.table(index = ele_dt$index, 
                                      in_range = rast_i$in_range[ele_dt$index])
}

ranges_all <- rbindlist(catch, use.names = TRUE, idcol = "species") %>%
    filter(in_range == 1)
saveRDS(ranges_all, "outputs/range_info_dt.rds")

## stars lookup ----
# raster for converting back to stars objects at a later date if useful
stars_lookup <- elevation %>% 
    mutate(index = 1:n(), EC_and_above_1000 = !is.na(ele_1km.tif)) %>% 
    select(-ele_1km.tif)
saveRDS(stars_lookup, "outputs/stars_lookup.rds")

## elev to relev ----
cu_cov <- readRDS("../sparing_sharing_elevation/cu_cov.rds") %>%
    select(species, elev_breadth, elev_median) %>%
    unique %>%
    as.data.table

relev = (elev_ALOS - elev_median)/elev_breadth * 1.61

ele_dt

ranges_all[ele_dt, elev_ALOS := ele_ALOS, on="index"]
ranges_all[cu_cov, `:=`(elev_median = elev_median, 
                        elev_breadth = elev_breadth), on="species"]

ranges_all[, relev := (elev_ALOS - elev_median)/elev_breadth * 1.61]
ranges_all[, elev_ALOS := NULL]

saveRDS(ranges_all, "outputs/spatial_input_layers.rds")
