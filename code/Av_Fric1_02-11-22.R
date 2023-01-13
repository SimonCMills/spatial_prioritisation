#####################################
#### Targeted Functional Indices ####
#####################################



#### Pacakges ####
library(tidyverse)
library(mFD)


### Read in data ####

## Occupancy data where last 50 columns are posterior samples, rows are species
Occ_dat <- read_rds("data/Birds/temp_output_occupancy_elev.rds")

## Trait data
HWI_dat <- read.csv("data/Birds/Dataset_HWI.csv") %>% 
  mutate(Species.name = str_replace(Species.name, " ", "_"), IUCN.name = str_replace(IUCN.name, " ", "_"))

Elton_dat <- read.delim("data/Birds/BirdFuncDat.txt")%>% 
  mutate(Scientific = str_replace(Scientific, " ", "_")) %>% 
  select(Scientific, Nocturnal)

#### Bind all traits ####

Trait_occ_HWI <- left_join(Occ_dat, HWI_dat, by = c("species_eltontraits" = "Tree.name"))
Trait_occ <- left_join(Trait_occ_HWI, Elton_dat, by = c("species_eltontraits" = "Scientific"))

## Format traits
Trait_db <- Trait_occ %>% ungroup() %>%
  ## Correct the naming oddity needs following up is this duplicate/ssp/regional species etc check with SCM
  mutate(species = ifelse(sp_obs1 %in% c("Grallaria_saturata__JBS", "Grallaria_saturata__SCM"), 
                          "Grallaria_saturata", species)) %>% 
  ## Select only species that are in the filter dataset
  filter(species %in% rownames(Bird_comm_filter)) %>%
  ## Select traits
  group_by(species,
           ## Diet
           Diet.Fruit, Diet.Vect, Diet.Inv, Diet.Vfish, Diet.Vunk,
           Diet.Scav, Diet.Nect, Diet.PlantO, Diet.Seed, Diet.Vend,
           ## Foraging
           ForStrat.aerial, ForStrat.ground, ForStrat.midhigh, ForStrat.canopy,
           ForStrat.watbelowsurf, ForStrat.wataroundsurf, ForStrat.understory, 
           ## Bodymass
           BodyMass.Value,
           ## Activity?
           # Nocturnal
           ## HWI
           HWI) %>% tally() %>% ungroup() %>% select(-n) %>%
  ## take the log of bodymass
  mutate(BodyMass = log(BodyMass.Value)) %>%
  ## prep data by setting species as rownames
  select(-BodyMass.Value) %>% column_to_rownames("species")

## Create a dataframe describing the structure of the trait database.
## What type are traits and which fuzzy traits correspond to what overall trait
Trait_Descr <- data.frame(trait_name = c("Diet.Fruit", "Diet.Vect", "Diet.Inv", "Diet.Vfish", "Diet.Vunk", 
                                         "Diet.Scav", "Diet.Nect", "Diet.PlantO", "Diet.Seed", "Diet.Vend",
                                         "ForStrat.aerial", "ForStrat.ground", "ForStrat.midhigh", "ForStrat.canopy",
                                         "ForStrat.watbelowsurf", "ForStrat.wataroundsurf", "ForStrat.understory",
                                         "BodyMass", "HWI"), 
                          trait_type = c("F", "F", "F", "F", "F", "F", "F", "F", "F", "F",
                                         "F", "F", "F", "F", "F", "F", "F",
                                         "Q", "Q"), 
                          fuzzy_name = c("Diet", "Diet", "Diet", "Diet", "Diet", 
                                         "Diet", "Diet", "Diet", "Diet", "Diet", 
                                         "Foraging", "Foraging", "Foraging", "Foraging",
                                         "Foraging", "Foraging", "Foraging",
                                         NA, NA))


## Computing distances between species based on functional traits
sp_dist <- funct.dist(Trait_db, Trait_Descr, metric = "gower")

## Compute multidimensional functional spaces and assess their quality for use in the loop to check whether limited by species number.
dist_qual <- quality.fspaces(sp_dist, maxdim_pcoa = 10,
                             deviation_weighting = "absolute", fdist_scaling = FALSE, fdendro = "average")

#### Binomial sampling ####

## Pivot the data
Occ_long <- Occ_dat  %>%
  filter(elev_ALOS > 880) %>%
  ## "Correct" this naming so all species are unique.
  ## Needs confirming with SCM whether distinct sp or ssp etc.
  mutate(species = ifelse(sp_obs1 %in% c("Grallaria_saturata__JBS", "Grallaria_saturata__SCM"), 
                          "Grallaria_saturata", species)) %>% 
  select(species, point, elev_ALOS, 158:207) %>% pivot_longer(cols = 4:53) 


## corrected number of 1173 sp
length(unique(Occ_long$species))

## Sample 10 communities from each iteration
## 50 iterations, 10 samples, 500 communities
Occ_sim <- Occ_long %>%
  group_by(species, point, name) %>% 
  summarize(Presence = rbinom(10, 1, value), sample = seq(from = 1, to = 10, length.out = 10)) %>%
  unite("iter_sample", name, sample, sep= "_", remove = FALSE) %>%
  mutate(group_id = cur_group_id())

Occ_sim <- Occ_sim %>% group_by(iter_sample) %>% mutate(group_id = cur_group_id())

Storage <- data.frame(matrix(NA, nrow = 0, ncol = 0))

## Move the pcoa out of the loop

for(i in 1:500)
{
  print(i)
  ## add indexing here when doing loops
  Occ_sample <- Occ_sim %>% filter(group_id == i) %>% ungroup() %>% select(species, point, Presence)
  Bird_comm <- pivot_wider(Occ_sample, names_from = point, values_from = Presence, values_fill = 0)
  
  ## Filter out species that never occur and sites that have no species
  ## NOTE 02/11/22: not 100% sure why I had this here (the site part is neccessary or the alpha.dim FD function will fail I remember,
  ## but I feel like including species that don't occur shouldnt actually change anything at all).
  Bird_comm_filter <- Bird_comm %>% column_to_rownames(var = "species") %>% 
    select(which(colSums(select(Bird_comm, - species)) > 0)) %>%
    filter_all(any_vars(. != 0)) %>% as.matrix()
  
  ## Transpose the bird data ready
  Bird_comm_t <- Bird_comm_filter %>% t()
  
  
  ## Extract the smallest number of species
  Min_sp_num <- rowSums(Bird_comm_t) %>% as.data.frame() %>% summarise(min(.))
  
  ## Select the number of axes based on the minimum number of species present at a site
  Optimum_pcoa <- round(dist_qual$"quality_fspaces", 3) %>% as.data.frame() %>% 
    ## add col for the number of pcoa axes
    mutate(pcoa = seq(from = 1, to = 11)) %>% 
    ## remove the tree average value to keep only axes 1 - 10
    filter(pcoa < 11) %>%
    ## filter to keep only axes less than the number of species
    filter(pcoa < Min_sp_num[,1]) %>%
    ## select the minimum number of axes by the MAD
    slice(which.min(mad))
  
  ## Get all axes coords
  sp_faxes_coord <- dist_qual$"details_fspaces"$"sp_pc_coord"
  
  ## Select the number of axes from coord object
  ## Coerce to matrix in the case of only axis defaults to vector
  selected_sp_axes_coord <- as.matrix(select(as.data.frame(sp_faxes_coord), c(1:Optimum_pcoa[,2])))
  
  
  ## Functional alpha diversity indices in a multidimensional space
  test <- alpha.fd.multidim(sp_faxes_coord = selected_sp_axes_coord,
                            asb_sp_w = Bird_comm_t,
                            ind_vect = c("fdis", "fmpd", "fnnd", "fric", "fori", "fspe"),
                            scaling = TRUE, check_input = TRUE, details_returned = TRUE)
  
  Metrics_iter <- test$functional_diversity_indices %>% rownames_to_column() %>% select(1:10) %>% mutate(iter = i, min_sp = Min_sp_num)
  Storage <- rbind(Storage, Metrics_iter)
  
  write_rds(Storage,
            paste0("Output/Func_metrics/Data/Out_",
                   paste(i), ".rds"))
}
