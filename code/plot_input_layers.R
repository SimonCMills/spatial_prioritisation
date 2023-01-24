library(stars); library(data.table); library(dplyr); library(sf)

## layers created in generate spatial layers script. 
elevation <- read_stars("data/ele_1km.tif")
pct_urban <- read_stars("data/pct_urban.tif")
pct_waterbody <- read_stars("data/pct_water.tif")

# note: these latter two are the percentages of forest in non-urban, 
# non-waterbody area, and the % treecover in non-urban, non-waterbody, non-forest
# areas respectively. 
pct_treecover <- read_stars("data/pct_treecover_in_pasture.tif")
pct_forest <- read_stars("data/pct_contig_forest.tif") %>%
    mutate(pct_forest = ifelse(is.na(pct_forest), 0, pct_forest))
# cells that contain 0 forest are masked out. Set these to 0
# pct_forest[pct_urban[] + pct_waterbody[] == 1] <- 0
all <- c(elevation, pct_urban, pct_waterbody, pct_treecover, pct_forest) %>%
    rename(pct_urban = pct_urban.tif, 
           ele = ele_1km.tif, 
           pct_water = pct_water.tif, 
           pct_treecover = pct_treecover_in_pasture.tif, 
           pct_forest = pct_contig_forest.tif) %>%
    mutate(pct_unmasked = 1 - (pct_urban + pct_water), 
           pct_forest = ifelse(is.na(pct_forest) & !is.na(ele), 0, pct_forest), 
           pct_forest2 = pct_forest * pct_unmasked, 
           pct_treecover = pct_treecover/100,
           pct_treecover = ifelse(is.na(pct_treecover) & !is.na(ele), 0, pct_treecover),
           pct_pasture = (1-pct_forest) * pct_unmasked)

bbox <- st_as_sfc(st_bbox(elevation))
dims <- st_bbox(elevation)
# other features ----
library(ggplot2); library(rnaturalearth)
flat <- "epsg:3117"

V_border <- ne_countries(country="Colombia", scale=10) %>% 
    st_as_sf() %>%
    st_transform(., st_crs(elevation)) %>%
    st_crop(., bbox) #%>%
    as_Spatial() %>%
    fortify(.)
    
ele_flat <- st_transform(elevation, flat)

Bogota_df <- data_frame(y=4.7110, x=-74.0721) %>%
    st_as_sf(., coords=c("x", "y"), crs="epsg:4326") %>%
    st_transform(., crs=st_crs(elevation)) %>%
    mutate(label = "Bogotá")

# get mountain ranges
mountains <- readRDS("../range_shifts/data/mountain_polygons.rds") %>% 
    #filter(range == "Cordillera_Oriental_Colombia_Venezuela") %>%
    st_union() %>%
    st_transform(., crs=st_crs(elevation))# %>%
    st_buffer(1) %>%
    st_crop(., country_borders) %>%
    st_transform(., crs=flat) %>%
    st_crop(CO_buffer)

st_bbox(ele_flat)

# elevation
p_elevation <- ggplot() +
    geom_sf(data=mountains, fill="grey", col=NA) +
    geom_sf(data=V_border, fill=NA, col="black") +
    geom_stars(data=elevation) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=17, size=4) +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.text = element_text(colour="black")) +
    labs(x="", y="", fill="", title = "(a) elevation")

# elevation
p_ele <- ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all["ele"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "(a) Elevation")

p_urban <- ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all["pct_urban"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "(b) Proportion urban")


p_forest <- ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all["pct_forest"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "(c) Proportion forest")

p_tc <- ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all["pct_treecover"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "(d) Mean wildlife friendliness")


plot_all <- egg::ggarrange(p_ele, p_urban, p_forest, p_tc, ncol=4)
ggsave("input_layers.png", plot=plot_all, dpi=200, units="mm")


all2 <- all %>%
    mutate(rsparrow = NA, wquail = NA)

rsparrow <- spatial_input_layers %>% filter(species == "Zonotrichia_capensis")
woodquail <- spatial_input_layers %>% filter(species == "Odontophorus_strophium")

all2[["rsparrow"]][rsparrow$index] <- 1
all2[["wquail"]][woodquail$index] <- 1

ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all2["rsparrow"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "") +
    guides(fill="none")

ggsave("Ruf_col_sparrow_range.png")

ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all2["wquail"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "") +
    guides(fill="none")

ggsave("Woodquail_range.png")


sp_i <- spatial_input_layers %>% 
    filter(species == "Zonotrichia_capensis")

all2[["rsparrow"]][sp_i$index] <- sp_i$amt_occ

ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all2["rsparrow"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "") +
    guides(fill="none")

ggsave("Ruf_col_sparrow_preds.png")



sp_i <- spatial_input_layers %>% 
    filter(species == "Odontophorus_strophium")

all2[["wquail"]][sp_i$index] <- sp_i$amt_occ

ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all2["wquail"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "") +
    guides(fill="none")

ggsave("Wquail_preds.png")

## lastly, do species richness
all2[["wquail"]] <- NA
all2[["wquail"]][new$index] <- new$V1

ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_stars(data=all2["wquail"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black")) +
    labs(x="", y="", fill="", title = "") 
ggsave("Species_richness.png")


### plot points 
df_ptInfo <- readRDS("../Colombia_rangeLims/data/point_ele_Eastern Cordillera.rds") %>%
    filter(!grepl("^CC", point), ele_jaxa >= 875) %>%
    st_transform(crs = st_crs(elevation))

ggplot() +
    geom_sf(data=mountains, fill="grey90", col=NA) +
    #geom_stars(data=all2["rsparrow"]) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=df_ptInfo) +
    geom_sf(data=Bogota_df, pch=24, size=3, fill="white",col="black") +
    geom_sf(data=V_border, fill=NA, col="black") +
    coord_sf(xlim=dims[c("xmin", "xmax")], 
             ylim = dims[c("ymin", "ymax")]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.background = element_rect(fill="grey96"), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position = "bottom", 
          legend.text = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x="", y="", fill="", title = "") +
    guides(fill="none")
ggsave("Survey_locations.png")


pts <- st_coordinates(df_ptInfo) %>% as_tibble



    
