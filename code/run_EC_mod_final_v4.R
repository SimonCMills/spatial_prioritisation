# version of model with elevationally varying woody veg effects. 

# housekeeping
library(dplyr); library(flocker); library(brms); library(ape)

# pass job id from bash
cpu_info <- as.numeric(commandArgs(T))
job_id <- cpu_info[1]

# read Jacob's canonical bird dataset
bird_data <- readRDS("data/Jacob_dropbox/birds.RDS")

# points in EC
c_phylo <- readRDS("output/consensus_phylogeny.rds")
taxonomy_lookup <- read.csv("data/Jacob_dropbox/taxonomy_lookup_initialVer.csv") %>%
    select(species = HBW, species_eltontraits = eltontraits) %>%
    mutate_all(function(x) gsub(" ", "_", x))

analysis_df <- readRDS("data/analysis_df_EC.rds") %>%
    left_join(., taxonomy_lookup) %>% 
    filter(species != "Columba_livia")

felicity_df <- read.csv("data/Felicity_data/ECpts_landscapemetricsMASTER.csv") %>%
    as_tibble %>%
    rename(dist_edge = dist_forest_edge, 
           woody_veg = woodyveg_index) %>%
    select(point, dist_edge, woody_veg)

adf2 <- bird_data %>% 
    filter(point %in% analysis_df$point) %>%
    filter(species != "Columba_livia") %>%
    mutate(species = ifelse(species == "Grallaria_saturata", "Grallaria_rufula", species),
          forestDep_birdlife = ifelse(species == "Grallaria_rufula", "Medium", forestDep_birdlife)) %>%
    left_join(., taxonomy_lookup) 

adf2 %>%
    filter(forestDep_birdlife == "unset") %>%
    pull(species) %>% unique


# all points in EC in felicity's df? TRUE
#all(unique(adf2$point) %in% felicity_df$point)

# all species fail to get matched with taxonomy? FALSE 
#any(is.na((adf2$species_eltontraits)))

# join in felicity's covariates
adf3 <- left_join(adf2, felicity_df) %>%
    filter(forestDep_birdlife != "unset") %>% # 4 species
    filter(!grepl("^CC", point)) %>%
    filter(paramo != 1)

sp_list <- unique(analysis_df$species_eltontraits)

drop_index <- which(!(c_phylo$tip.label %in% adf3$species_eltontraits))
phylo_subset <- drop.tip(c_phylo, drop_index)
A <- vcv(phylo_subset)
A <- A/max(A)

obsvr_ptLvl <- adf3 %>%
    select(point, obs1:obs4) %>%
    unique %>%
    mutate(observer_pt = apply(select(., -point), 1, function(x) names(which.max(table(x))))) %>%
    select(-(obs1:obs4))

# pt_subset <- adf3 %>% 
#     select(point, cluster) %>% 
#     unique %>% 
#     group_by(cluster) %>%
#     slice(1:2)

in_range_or_det <- adf3 %>%
    group_by(species) %>%
    summarise(in_range = any(distance_from_range <= 0), 
              detected = any(Q == 1)) %>%
    filter(in_range | detected)

adf4 <- adf3 %>%
    # filter(point %in% pt_subset$point)
    filter(species %in% in_range_or_det$species)

# extract detection matrix
det <- adf4 %>% select(v1:v4) %>% 
    as.matrix

cu_cov <- adf4 %>%
    left_join(., obsvr_ptLvl) %>%
    # note: have removed paramo, so all pasture or forest
    mutate(forest = ifelse(natural == 1 & paramo == 0, 1, -1), 
           pasture = ifelse(pasture==1, 1, -1),
           #paramo = ifelse(paramo == 1, 1, -1),
           obs_sp = interaction(observer_pt, species), 
           dist_edge_sc = scale(dist_edge), 
           woody_veg_sc = scale(woody_veg), 
           mass = BodyMass.Value,
           log_mass_sc = scale(log(mass)), 
           elev_breadth_sc = scale(elev_breadth),
           elev_median_sc = scale(elev_median), 
           diet_inv = ifelse(Diet.5Cat == "Invertebrate", 1, -1),
           diet_carn = ifelse(Diet.5Cat == "VertFishScav", 1, -1),
           diet_fruitnect = ifelse(Diet.5Cat == "FruiNect", 1, -1),
           diet_gran = ifelse(Diet.5Cat == "PlantSeed", 1, -1), 
           strata_index_tot = ForStrat.watbelowsurf + ForStrat.wataroundsurf + 
               ForStrat.ground + ForStrat.understory + ForStrat.midhigh + ForStrat.canopy,
           strata_index = (ForStrat.understory*(1/3) + ForStrat.midhigh*(2/3) + 
                               ForStrat.canopy*1)/strata_index_tot, 
           strata_index = ifelse(strata_index_tot == 0, 0, strata_index),
           strata_aerial = ForStrat.aerial/100, 
           relev = (elev_ALOS - elev_median)/(elev_breadth/2),
           relev2 = relev^2) %>%
    select(point, 
           subregion,
           species, 
           phylo = species_eltontraits,
           sp_cl, 
           sp_sr,
           sr = subregion,
           obs_sp,
           pasture, 
           paramo, 
           dist_edge, 
           dist_edge_sc,
           woody_veg,
           woody_veg_sc,
           forest_dep = forestDep_birdlife,
           elev_breadth, 
           elev_breadth_sc,
           elev_median, 
           elev_median_sc,
           mass,
           log_mass_sc,
           diet_inv,
           diet_carn,
           diet_fruitnect,
           diet_gran,
           strata_index, 
           strata_aerial,
           relev, 
           relev2,
           observer_pt)

visit_cov <- adf4 %>% select(hps1:hps4) %>%
    as.matrix %>%
    as.numeric %>%
    scale %>%
    matrix(., ncol=4) %>%
    list(time = .)

fd <- make_flocker_data(det, cu_cov, event_covs = visit_cov)
fit_priors <- c(brms::set_prior("normal(0, 2)", class = "b"), 
                brms::set_prior("logistic(0, 1)", class="Intercept"))

## remove redundant objects (frees up ~1.9GB) 
rm(list = c("det", "visit_cov", "cu_cov", "bird_data", "analysis_df", "adf2", 
            "adf3", "adf4", "taxonomy_lookup"))
            
## fit model 
dirname <- paste0(getwd(), "/stan_out/EC_model_run5_woody_veg_effs")
if(!dir.exists(dirname)) dir.create(dirname)
print(paste0("directory: ", dirname))

writename <- paste0(format(Sys.time(), "warmup_%d-%m-%Y_%H%M_j"), job_id)

saveRDS(A, paste0(dirname, "/A.rds"))
saveRDS(fd, paste0(dirname, "/fd.rds"))

fit <- flock(
    # occupancy
    f_occ = ~ 1 + forest_dep + pasture +
        pasture:forest_dep +
        woody_veg_sc + 
        woody_veg_sc:forest_dep +
        # dietary traits
        diet_inv + diet_inv:pasture + diet_inv:woody_veg_sc +
        diet_carn + diet_carn:pasture + diet_carn:woody_veg_sc +
        diet_fruitnect + diet_fruitnect:pasture + diet_fruitnect:woody_veg_sc +
        diet_gran + diet_gran:pasture + diet_gran:woody_veg_sc +
        # mass and foraging traits
        log_mass_sc + log_mass_sc:pasture + log_mass_sc:woody_veg_sc +
        strata_index + strata_index:pasture + strata_index:woody_veg_sc +
        strata_aerial + strata_aerial:pasture + strata_aerial:woody_veg_sc +
        # elevational breadth
        elev_breadth_sc + elev_breadth_sc:pasture + elev_breadth_sc:woody_veg_sc +
        # elevational position
        elev_median_sc + 
        elev_median_sc:forest_dep +
        elev_median_sc:pasture + 
        elev_median_sc:woody_veg_sc + 
        elev_median_sc:pasture:forest_dep +
        elev_median_sc:woody_veg_sc:forest_dep + 
        # elevational range
        relev + relev2 + 
        # ranefs
        (1 + pasture + woody_veg_sc + relev + relev2||species) + 
        (1|sr) + (1|sp_sr) + (1|sp_cl) +
        (1 + pasture|gr(phylo, cov = A)), # try to fit phylo woody veg too?
    # detection            
    f_det = ~ 1 + time + observer_pt + strata_index + strata_aerial + pasture + 
        woody_veg_sc + (1 + pasture|species) + (1|obs_sp), 
    # flocker info            
    data2 = list(A = A),
    flocker_data = fd, 
    rep_constant = F, 
    # run & save info 
    backend = "cmdstanr", 
    chains = 4, cores = 4, 
    max_treedepth = 9, 
    step_size = 0.02,
    iter=1500, warmup=1000,
    file = paste0(dirname, "/", writename),
    save_warmup = FALSE,
    output_dir = dirname, 
    output_basename = writename, 
    threads_per_chain = 1)

# write model name to file
write(fit$fit@model_name, file = paste0(dirname, "/completed_warmup.txt"), append=T)