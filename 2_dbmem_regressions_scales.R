#### 2. dbMEM + regressions at different scales ####

library(tidyverse)
library(adespatial) # v0.3.8
library(vegan) # v2.5.6
library(brms) # v2.13.0
library(tidybayes) # v2.1.1
library(performance) # v0.4.7
library(viridis)
library(extrafont)
library(patchwork)

#### Data ####
amph <- read.csv("amph_data.csv", head = TRUE)

amph$urban <- as.factor(amph$urban)
amph$species <- as.factor(amph$species)

## Subsets ##
# Gene diversity and allelic richness subsets
amphxy <- dplyr::select(amph, lon, lat)

# Ne subsets
amph.ne <- amph[!(is.na(amph$Ne)),]
amphxy.ne <- dplyr::select(amph.ne, lon, lat)

# Fst subsets
amph.fst <- amph[!is.na(amph$global_fst),]
amphxy.fst <- dplyr::select(amph.fst, lon, lat)

#### dbMEM analysis ####
# Check for linear gradient
anova(lm(amph$gene_diversity ~ ., data=amphxy)) # yes
anova(lm(amph$allelic_richness ~ ., data=amphxy)) # yes
anova(lm(amph.ne$Ne ~ ., data=amphxy.ne)) # no
anova(lm(amph.fst$global_fst ~ ., data=amphxy.fst)) # yes

# Detrend data
amph.det.gd  <-  resid(lm(amph$gene_diversity ~ .,   data=amphxy)) # gene diversity
amph.det.ar  <-  resid(lm(amph$allelic_richness ~ ., data=amphxy)) # allelic richness
amph.det.fst <-  resid(lm(amph.fst$global_fst ~ ., data=amphxy.fst)) # Fst

# Construct the matrix of dbMEM variables
amph.dbmem <-    as.data.frame(dbmem(amphxy, silent=FALSE)) # same sites for GD and AR
amph.dbmem.ne <- as.data.frame(dbmem(amphxy.ne, silent=FALSE))
amph.dbmem.fst <-as.data.frame(dbmem(amphxy.fst, silent=FALSE))

# Select significant MEMs: Global significance test
### gene diversity
GDlm <- lm(amph.det.gd ~., amph.dbmem)
summary(GDlm)

(GDr2da <- RsquareAdj(GDlm)$adj.r.squared)
# forward selection
GDmemfwd <- forward.sel(amph.det.gd, as.matrix(amph.dbmem), 
                        adjR2thresh = GDr2da)
# sort & extract selected MEMs
(GDmems <- sort(GDmemfwd[,2]))

### allelic richness
ARlm <- lm(amph.det.ar ~., amph.dbmem)
summary(ARlm)

(ARr2da <- RsquareAdj(ARlm)$adj.r.squared)
# forward selection
ARmemfwd <- forward.sel(amph.det.ar, as.matrix(amph.dbmem), 
                        adjR2thresh = ARr2da)
# sort & extract selected MEMs
(ARmems <- sort(ARmemfwd[,2]))

### effective population size
NElm <- lm(amph.ne$Ne ~., amph.dbmem.ne)
summary(NElm) # none

### FST
FSTlm <- lm(amph.det.fst ~., amph.dbmem.fst)
summary(FSTlm)

(FSTr2da <- RsquareAdj(FSTlm)$adj.r.squared)
# forward selection
FSTmemfwd <- forward.sel(amph.det.fst, as.matrix(amph.dbmem.fst), 
                         adjR2thresh = FSTr2da)
# sort & extract selected MEMs
(FSTmems <- sort(FSTmemfwd[,2]))

#### brms models ####
### Data
# combine dbMEMs with gdata. Gene diversity & allelic richness:
amph.db <- cbind(amph, amph.dbmem)

# Ne
amph.db.ne <- cbind(amph.ne, amph.dbmem.ne)

# FST
amph.db.fst <- cbind(amph.fst, amph.dbmem.fst)

## BRM 
brmfun <- function(model_spec, data_group){
  brm(model_spec, cores = 4, chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.999, max_treedepth = 15),
      data = data_group)
}

#### gene diversity ####
## Urban/rural
UR_gd_1 <- bf(scale(gene_diversity) ~ urban + 
                MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (urban|species))
## population density
PD_gd_1km <- bf(scale(gene_diversity) ~ scale(log(pd_1_km+1)) + 
                  MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(log(pd_1_km+1))|species))
PD_gd_5km <- bf(scale(gene_diversity) ~ scale(log(pd_5_km+1)) + 
                  MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(log(pd_5_km+1))|species))
PD_gd_10km <- bf(scale(gene_diversity) ~ scale(log(pd_10_km+1)) + 
                 MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(log(pd_10_km+1))|species))
PD_gd_15km <- bf(scale(gene_diversity) ~ scale(log(pd_15_km+1)) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(log(pd_15_km+1))|species))
## HFI
HFI_gd_1km <- bf(scale(gene_diversity) ~ scale(hfi_1_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(hfi_1_km)|species))
HFI_gd_5km <- bf(scale(gene_diversity) ~ scale(hfi_5_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(hfi_5_km)|species))
HFI_gd_10km <- bf(scale(gene_diversity) ~ scale(hfi_10_km) + 
              MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(hfi_10_km)|species))
HFI_gd_15km <- bf(scale(gene_diversity) ~ scale(hfi_15_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(hfi_15_km)|species))
## Roads
road_gd_1km <- bf(scale(gene_diversity) ~ scale(road_1_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(road_1_km)|species))
road_gd_5km <- bf(scale(gene_diversity) ~ scale(road_5_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(road_5_km)|species))
road_gd_10km <- bf(scale(gene_diversity) ~ scale(road_10_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(road_10_km)|species))
road_gd_15km <- bf(scale(gene_diversity) ~ scale(road_15_km) + 
                     MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(road_15_km)|species))

## % Impervious surface
IS_gd_1km <- bf(scale(gene_diversity) ~ scale(imp_1_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(imp_1_km)|species))
IS_gd_5km <- bf(scale(gene_diversity) ~ scale(imp_5_km) + 
                  MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(imp_5_km)|species))
IS_gd_10km <- bf(scale(gene_diversity) ~ scale(imp_10_km) + 
               MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(imp_10_km)|species))
IS_gd_15km <- bf(scale(gene_diversity) ~ scale(imp_15_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(imp_15_km)|species))
IS_gd_30m <- bf(scale(gene_diversity) ~ scale(imp_30m) +
               MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(imp_30m)|species))

## None ##
NU_gd <- bf(scale(gene_diversity) ~ MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (1|species))

uGD1 <-brmfun(UR_gd_1, amph.db)
pGD1 <-brmfun(PD_gd_1km, amph.db)
pGD2 <-brmfun(PD_gd_5km, amph.db)
pGD3 <-brmfun(PD_gd_10km, amph.db)
pGD4 <-brmfun(PD_gd_15km, amph.db)
hGD1 <-brmfun(HFI_gd_1km, amph.db)
hGD2 <-brmfun(HFI_gd_5km, amph.db)
hGD3 <-brmfun(HFI_gd_10km, amph.db)
hGD4 <-brmfun(HFI_gd_15km, amph.db)
rGD1 <-brmfun(road_gd_1km, amph.db)
rGD2 <-brmfun(road_gd_5km, amph.db)
rGD3 <-brmfun(road_gd_10km, amph.db)
rGD4 <-brmfun(road_gd_15km, amph.db)
iGD1 <-brmfun(IS_gd_1km, amph.db)
iGD2 <-brmfun(IS_gd_5km, amph.db)
iGD3 <-brmfun(IS_gd_10km, amph.db)
#iGD4 <-brmfun(IS_gd_15km, amph.db) # didn't run; increase iterations:
iGD4 <-brm(IS_gd_15km, cores = 4, chains = 4, iter = 10000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
    data = amph.db)
iGD5 <-brmfun(IS_gd_30m, amph.db)
nGD1 <-brmfun(NU_gd, amph.db)

#### allelic richness ####
## urban/rural
UR_ar_1 <- bf(scale(allelic_richness) ~ urban + 
               MEM2 + MEM3 + MEM4 + MEM7 + (urban|species))

## population density
PD_ar_1km <- bf(scale(allelic_richness) ~ scale(log(pd_1_km+1)) + 
                  MEM2 + MEM3 + MEM4 + MEM7 + (scale(log(pd_1_km+1))|species))
PD_ar_5km <- bf(scale(allelic_richness) ~ scale(log(pd_5_km+1)) + 
                  MEM2 + MEM3 + MEM4 + MEM7 + (scale(log(pd_5_km+1))|species))
PD_ar_10km <- bf(scale(allelic_richness) ~ scale(log(pd_10_km+1)) + 
             MEM2 + MEM3 + MEM4 + MEM7 + (scale(log(pd_10_km+1))|species))
PD_ar_15km <- bf(scale(allelic_richness) ~ scale(log(pd_15_km+1)) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + (scale(log(pd_15_km+1))|species))

## HFI
HFI_ar_1km <- bf(scale(allelic_richness) ~ scale(hfi_1_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + (scale(hfi_1_km)|species))
HFI_ar_5km <- bf(scale(allelic_richness) ~ scale(hfi_5_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + (scale(hfi_5_km)|species))
HFI_ar_10km <- bf(scale(allelic_richness) ~ scale(hfi_10_km) + 
              MEM2 + MEM3 + MEM4 + MEM7 + (scale(hfi_10_km)|species))
HFI_ar_15km <- bf(scale(allelic_richness) ~ scale(hfi_15_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + (scale(hfi_15_km)|species))

## roads
road_ar_1km <- bf(scale(allelic_richness) ~ scale(road_1_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + (scale(road_1_km)|species))
road_ar_5km <- bf(scale(allelic_richness) ~ scale(road_5_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + (scale(road_5_km)|species))
road_ar_10km <- bf(scale(allelic_richness) ~ scale(road_10_km) + 
                    MEM2 + MEM3 + MEM4 + MEM7 + (scale(road_10_km)|species))
road_ar_15km <- bf(scale(allelic_richness) ~ scale(road_15_km) + 
                     MEM2 + MEM3 + MEM4 + MEM7 + (scale(road_15_km)|species))

## % impervious surface
IS_ar_1km <- bf(scale(allelic_richness) ~ scale(imp_1_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + (scale(imp_1_km)|species))
IS_ar_5km <- bf(scale(allelic_richness) ~ scale(imp_5_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + (scale(imp_5_km)|species))
IS_ar_10km <- bf(scale(allelic_richness) ~ scale(imp_10_km) + 
               MEM2 + MEM3 + MEM4 + MEM7 + (scale(imp_10_km)|species))
IS_ar_15km <- bf(scale(allelic_richness) ~ scale(imp_15_km) + 
                   MEM2 + MEM3 + MEM4 + MEM7 + (scale(imp_15_km)|species))
IS_ar_30m <- bf(scale(allelic_richness) ~ scale(imp_30m) + 
               MEM2 + MEM3 + MEM4 + MEM7 + (scale(imp_30m)|species))

## None
NU_ar <- bf(scale(allelic_richness) ~ MEM2 + MEM3 + MEM4 + MEM7 + (1|species))

uAR1 <- brmfun(UR_ar_1, amph.db)
pAR1 <- brmfun(PD_ar_1km, amph.db)
pAR2 <- brmfun(PD_ar_5km, amph.db)
pAR3 <- brmfun(PD_ar_10km, amph.db)
pAR4 <- brmfun(PD_ar_15km, amph.db)
hAR1 <-brmfun(HFI_ar_1km, amph.db)
hAR2 <-brmfun(HFI_ar_5km, amph.db)
hAR3 <-brmfun(HFI_ar_10km, amph.db)
hAR4 <-brmfun(HFI_ar_15km, amph.db)
rAR1 <-brmfun(road_ar_1km, amph.db)
rAR2 <-brmfun(road_ar_5km, amph.db)
rAR3 <-brmfun(road_ar_10km, amph.db)
rAR4 <-brmfun(road_ar_15km, amph.db)
iAR1 <-brmfun(IS_ar_1km, amph.db)
iAR2 <-brmfun(IS_ar_5km, amph.db)
iAR3 <-brmfun(IS_ar_10km, amph.db)
iAR4 <-brmfun(IS_ar_15km, amph.db)
iAR5 <-brmfun(IS_ar_30m, amph.db)
nAR1 <-brmfun(NU_ar, amph.db)

#### effective population size ####
## Urban/rural
UR_ne_1 <- bf(scale(log(Ne)) ~ urban +(urban|species))

## population density
PD_ne_1km <- bf(scale(log(Ne)) ~ scale(log(pd_1_km+1)) + (scale(log(pd_1_km+1))|species))
PD_ne_5km <- bf(scale(log(Ne)) ~ scale(log(pd_5_km+1)) + (scale(log(pd_5_km+1))|species))
PD_ne_10km <- bf(scale(log(Ne)) ~ scale(log(pd_10_km+1)) + (scale(log(pd_10_km+1))|species))
PD_ne_15km <- bf(scale(log(Ne)) ~ scale(log(pd_15_km+1)) + (scale(log(pd_15_km+1))|species))

# HFI
HFI_ne_1km <- bf(scale(log(Ne)) ~ scale(hfi_1_km) + (scale(hfi_1_km)|species))
HFI_ne_5km <- bf(scale(log(Ne)) ~ scale(hfi_5_km) + (scale(hfi_5_km)|species))
HFI_ne_10km <- bf(scale(log(Ne)) ~ scale(hfi_10_km) + (scale(hfi_10_km)|species))
HFI_ne_15km <- bf(scale(log(Ne)) ~ scale(hfi_15_km) + (scale(hfi_15_km)|species))

# Roads
road_ne_1km  <- bf(scale(log(Ne)) ~ scale(road_1_km) + (scale(road_1_km)|species))
road_ne_5km  <- bf(scale(log(Ne)) ~ scale(road_5_km) + (scale(road_5_km)|species))
road_ne_10km <- bf(scale(log(Ne)) ~ scale(road_10_km) + (scale(road_10_km)|species))
road_ne_15km <- bf(scale(log(Ne)) ~ scale(road_15_km) + (scale(road_15_km)|species))

# % impervious surface
## Note some models did not convergence with random slopes
IS_ne_1km <- bf(scale(log(Ne)) ~ scale(imp_1_km) + (1|species)) #
IS_ne_5km <- bf(scale(log(Ne)) ~ scale(imp_5_km) + (scale(imp_5_km)|species))
IS_ne_10km <- bf(scale(log(Ne)) ~ scale(imp_10_km) + (scale(imp_10_km)|species))
IS_ne_15km <- bf(scale(log(Ne)) ~ scale(imp_15_km) + (scale(imp_15_km)|species))
IS_ne_30m <- bf(scale(log(Ne)) ~ scale(imp_30m) + (1|species)) #

# None
NU_ne <- bf(scale(log(Ne)) ~ (1|species))

uNE1 <- brmfun(UR_ne_1,  amph.db.ne)
pNE1 <- brmfun(PD_ne_1km,  amph.db.ne)
pNE2 <- brmfun(PD_ne_5km,  amph.db.ne)
pNE3 <- brmfun(PD_ne_10km, amph.db.ne)
pNE4 <- brmfun(PD_ne_15km, amph.db.ne)
hNE1 <-brmfun(HFI_ne_1km,  amph.db.ne)
hNE2 <-brmfun(HFI_ne_5km,  amph.db.ne)
hNE3 <-brmfun(HFI_ne_10km, amph.db.ne)
hNE4 <-brmfun(HFI_ne_15km, amph.db.ne)
rNE1 <-brmfun(road_ne_1km,  amph.db.ne)
rNE2 <-brmfun(road_ne_5km,  amph.db.ne)
rNE3 <-brmfun(road_ne_10km, amph.db.ne)
rNE4 <-brmfun(road_ne_15km, amph.db.ne)
iNE1 <-brmfun(IS_ne_1km, amph.db.ne)
iNE2 <-brmfun(IS_ne_5km, amph.db.ne)
iNE3 <-brmfun(IS_ne_10km, amph.db.ne)
iNE4 <-brmfun(IS_ne_15km, amph.db.ne)
iNE5 <-brmfun(IS_ne_30km, amph.db.ne)
nNE1 <-brmfun(NU_ne, amph.db.ne)

#### FST ####
# urban/rural
UR_fst_1 <- bf(scale(global_fst) ~ urban + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (urban|species))

# Population density
PD_fst_1km <- bf(scale(global_fst) ~ scale(log(pd_1_km+1)) + 
                   MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(log(pd_1_km+1))|species))
PD_fst_5km  <- bf(scale(global_fst) ~ scale(log(pd_5_km+1)) + 
                    MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(log(pd_5_km+1))|species))
PD_fst_10km <- bf(scale(global_fst) ~ scale(log(pd_10_km+1)) + 
               MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(log(pd_10_km+1))|species))
PD_fst_15km <- bf(scale(global_fst) ~ scale(log(pd_15_km+1)) + 
                    MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(log(pd_15_km+1))|species))


# HFI
HFI_fst_1km <- bf(scale(global_fst) ~ scale(hfi_1_km) +
                    MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(hfi_1_km)|species))
HFI_fst_5km <- bf(scale(global_fst) ~ scale(hfi_5_km) +
                    MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(hfi_5_km)|species))
HFI_fst_10km <- bf(scale(global_fst) ~ scale(hfi_10_km) +
               MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(hfi_10_km)|species))
HFI_fst_15km <- bf(scale(global_fst) ~ scale(hfi_15_km) +
                     MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(hfi_15_km)|species))

# roads
roads_fst_1km <- bf(scale(global_fst) ~ scale(road_1_km) +
                      MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(road_1_km)|species))
roads_fst_5km <- bf(scale(global_fst) ~ scale(road_5_km) +
                      MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(road_5_km)|species))
roads_fst_10km <- bf(scale(global_fst) ~ scale(road_10_km) +
                      MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(road_10_km)|species))
roads_fst_15km <- bf(scale(global_fst) ~ scale(road_15_km) +
                       MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + MEM24 + (scale(road_15_km)|species))

# % impervious surface
IS_fst_1km <- bf(scale(global_fst) ~ scale(imp_1_km) + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(imp_1_km)|species))
IS_fst_5km <- bf(scale(global_fst) ~ scale(imp_5_km) + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(imp_5_km)|species))
IS_fst_10km <- bf(scale(global_fst) ~ scale(imp_10_km) + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(imp_10_km)|species))
IS_fst_15km <- bf(scale(global_fst) ~ scale(imp_15_km) + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(imp_15_km)|species))
IS_fst_30m <- bf(scale(global_fst) ~ scale(imp_30m) + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(imp_30m)|species))

# None
NU_fst <- bf(scale(global_fst) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (1|species))

uFST1 <- brmfun(UR_fst_1,  amph.db.fst)
pFST1 <- brmfun(PD_fst_1km,  amph.db.fst)
pFST2 <- brmfun(PD_fst_5km,  amph.db.fst)
pFST3 <- brmfun(PD_fst_10km, amph.db.fst)
pFST4 <- brmfun(PD_fst_15km, amph.db.fst)
hFST1 <-brmfun(HFI_fst_1km,  amph.db.fst)
hFST2 <-brmfun(HFI_fst_5km,  amph.db.fst)
hFST3 <-brmfun(HFI_fst_10km, amph.db.fst)
hFST4 <-brmfun(HFI_fst_15km, amph.db.fst)
rFST1 <-brmfun(roads_fst_1km,  amph.db.fst)
rFST2 <-brmfun(roads_fst_5km,  amph.db.fst)
rFST3 <-brmfun(roads_fst_10km, amph.db.fst)
rFST4 <-brmfun(roads_fst_15km, amph.db.fst)
iFST1 <-brmfun(IS_fst_1km, amph.db.fst)
iFST2 <-brmfun(IS_fst_5km, amph.db.fst)
iFST3 <-brmfun(IS_fst_10km, amph.db.fst)
#iFST4 <-brmfun(IS_fst_15km, amph.db.fst) # didn't run; increase iterations:
iFST4 <-brm(IS_fst_15km, cores = 4, chains = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
           data = amph.db.fst)
iFST5 <-brmfun(IS_fst_30m, amph.db.fst)
nFST1 <-brmfun(NU_fst, amph.db.fst)

#### Bayes R2 (performance) ####
lapply(list(uGD1, pGD3, hGD3, iGD5, AGD5, nGD1,
            uAR1, pAR3, hAR3, iAR5, AAR5, nAR1,
            uNE1, pNE3, hNE3, iNE5, ANE5, nNE1,
            uFST1, pFST3, hFST3, iFST5, AFST5, nFST1), r2_bayes)

#### Plots ####
# plot overall effects:
plotdf <- function(x, response){data.frame(Variable = rownames(summary(x)$fixed)[2],
                                           Coefficient = summary(x)$fixed[2,1],
                                           CIlo95 = summary(x)$fixed[2,3],
                                           CIup95 = summary(x)$fixed[2,4],
                                           CIlo90 = posterior_interval(x, pars = get_variables(x)[2], prob = 0.90)[,1],
                                           CIup90 = posterior_interval(x, pars = get_variables(x)[2], prob = 0.90)[,2],
                                           Response_var = response)
}

# species-specific effect plots:
spp_plotdat <- function(data){
  data %>% 
    mutate(Order = ifelse(species %in% c("Ambystoma_barbouri", "Ambystoma_maculatum",
                                         "Desmognathus_fuscus", "Dicamptodon_aterrimus",
                                         "Dicamptodon_copei", "Ensatina_eschscholtzii",
                                         "Hydromantes_brunus", "Hydromantes_platycephalus",
                                         "Plethodon_albagula", "Plethodon_cinereus", 
                                         "Taricha_granulosa"), "Caudata", "Anura")) %>% 
    mutate(species = gsub("_", " ", species)) %>% 
    group_by(Order) %>% 
    arrange(species, .by_group = TRUE) %>% 
    mutate(species = factor(species, levels = unique(.$species)))
}

spp_plot <- function(data, xlab){
  ggplot(data, aes(y = species, x = species_mean, color=Order)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    stat_pointinterval(.width = c(0.9, 0.95)) +
    theme_minimal(base_size = 11) +
    labs(title = "", 
         x = xlab,
         y = "") +
    scale_color_viridis(option="D", discrete=TRUE, begin = 0.3, end=0.8) + 
    scale_y_discrete(limits = rev(levels(data$species))) +
    theme(text=element_text(family="Lato Black"),
          axis.text.y = element_text(face = "italic"))
}

#### Main text coefficient plot (Fig. 2) ####
m1ar <- plotdf(AAR1, "allelic richness")
m2ar <- plotdf(AAR2, "allelic richness")
m3ar <- plotdf(AAR3, "allelic richness")
m4ar <- plotdf(AAR7, "allelic richness")
m1gd <- plotdf(AGD1, "gene diversity")
m2gd <- plotdf(AGD2, "gene diversity")
m3gd <- plotdf(AGD3, "gene diversity")
m4gd <- plotdf(AGD7, "gene diversity")
m1ne <- plotdf(ANE1, "effective population size")
m2ne <- plotdf(ANE2, "effective population size")
m3ne <- plotdf(ANE3, "effective population size")
m4ne <- plotdf(ANE7, "effective population size")
m1fst <- plotdf(AFST1, "FST")
m2fst <- plotdf(AFST2, "FST")
m3fst <- plotdf(AFST3, "FST")
m4fst <- plotdf(AFST7, "FST")

allModelFrame.amph <- data.frame(rbind(m1ar, m2ar, m3ar, m4ar,
                                       m1gd, m2gd, m3gd, m4gd,
                                       m1ne, m2ne, m3ne, m4ne,
                                       m1fst, m2fst, m3fst, m4fst))

allModelFrame.amph$Variablecol <- rep(c("urban/rural", "human population density", 
                                    "Human Footprint Index", "% impervious surface*"), 4)

allModelFrame.amph$Response_var <- factor(allModelFrame.amph$Response_var, 
                                          levels = c("allelic richness", "gene diversity", 
                                                     "effective population size", "FST")) # reorder factor levels for plot
# plot
mp2 <- ggplot(allModelFrame.amph, aes(colour = Response_var))

mp2 + geom_hline(yintercept=seq(-1.5, 1.5, 0.5),  # x axis lines
                 lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = Variablecol, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = Variablecol, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame.amph$Variablecol))-0.5, 1), # y axis lines appear btwn groups
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE),
                      labels=c("allelic richness", "gene diversity", "effective population size", bquote(F[ST]))) +
  labs(x= "", title = "Model coefficients") +
  theme(text=element_text(family="Lato Black"))

#### Species-specific HFI plot (Fig. 3) ####
gdphfi <- hGD3 %>%
  spread_draws(b_scalehfi_10_km, r_species[species,]) %>%
  mutate(species_mean = b_scalehfi_10_km + r_species) 
gdphfi <- spp_plotdat(gdphfi) 

plot1 <- spp_plot(gdphfi, "\u03b2 (Human Footprint Index)") + 
  labs(title = 'Species-specific effects of HFI on gene diversity')


plot2 <- amph %>%
  mutate(species = gsub("_", " ", species)) %>% 
  group_by(order) %>% 
  arrange(species, .by_group = TRUE) %>% 
  mutate(species = factor(species, levels = unique(.$species))) %>% 
  ggplot(aes(x=hfi, y=gene_diversity, color=order)) +
  geom_point()+
  geom_smooth(method="lm")+
  scale_color_viridis(discrete=T, begin = 0.3, end=0.8,
                      guide = FALSE) +
  labs(x="Human Footprint Index", y="gene diversity") +
  facet_wrap(~species) +
  theme_minimal() +
  theme(text=element_text(family="Lato Black")) +
  theme(strip.text = element_text(face = "italic"))

layout <- "
AABBBB"
plot1 + plot2 + plot_layout(guides = "collect", design=layout)

#### Scales plot (Fig. S1) ####
gdpd1km <-   plotdf(pGD1,   "gene diversity")
gdpd5km <-   plotdf(pGD2,   "gene diversity")
gdpd10km <-  plotdf(pGD3,  "gene diversity")
gdpd15km <-  plotdf(pGD4,  "gene diversity")
gdhfi1km <-  plotdf(hGD1,  "gene diversity")
gdhfi5km <-  plotdf(hGD2,  "gene diversity")
gdhfi10km <- plotdf(hGD3, "gene diversity")
gdhfi15km <- plotdf(hGD4, "gene diversity")
gdro1km <-   plotdf(rGD1,  "gene diversity")
gdro5km <-   plotdf(rGD2,  "gene diversity")
gdro10km <-  plotdf(rGD3, "gene diversity")
gdro15km <-  plotdf(rGD4, "gene diversity")
gdis1km <-   plotdf(iGD1,  "gene diversity")
gdis5km <-   plotdf(iGD2,  "gene diversity")
gdis10km <-  plotdf(iGD3, "gene diversity")
gdis15km <-  plotdf(iGD4, "gene diversity")

arpd1km <-   plotdf(pAR1,   "allelic richness")
arpd5km <-   plotdf(pAR2,   "allelic richness")
arpd10km <-  plotdf(pAR3,  "allelic richness")
arpd15km <-  plotdf(pAR4,  "allelic richness")
arhfi1km <-  plotdf(hAR1,  "allelic richness")
arhfi5km <-  plotdf(hAR2,  "allelic richness")
arhfi10km <- plotdf(hAR3, "allelic richness")
arhfi15km <- plotdf(hAR4, "allelic richness")
arro1km <-   plotdf(rAR1,  "allelic richness")
arro5km <-   plotdf(rAR2,  "allelic richness")
arro10km <-  plotdf(rAR3, "allelic richness")
arro15km <-  plotdf(rAR4, "allelic richness")
aris1km <-   plotdf(iAR1,  "allelic richness")
aris5km <-   plotdf(iAR2,  "allelic richness")
aris10km <-  plotdf(iAR3, "allelic richness")
aris15km <-  plotdf(iAR4, "allelic richness")

nepd1km <-   plotdf(pNE1,   "effective population size")
nepd5km <-   plotdf(pNE2,   "effective population size")
nepd10km <-  plotdf(pNE3,  "effective population size")
nepd15km <-  plotdf(pNE4,  "effective population size")
nehfi1km <-  plotdf(hNE1,  "effective population size")
nehfi5km <-  plotdf(hNE2,  "effective population size")
nehfi10km <- plotdf(hNE3, "effective population size")
nehfi15km <- plotdf(hNE4, "effective population size")
nero1km <-  plotdf(rNE1,  "effective population size")
nero5km <-  plotdf(rNE2,  "effective population size")
nero10km <- plotdf(rNE3, "effective population size")
nero15km <- plotdf(rNE4, "effective population size")
neis1km <-  plotdf(iNE1,  "effective population size")
neis5km <-  plotdf(iNE2,  "effective population size")
neis10km <- plotdf(iNE3, "effective population size")
neis15km <- plotdf(iNE4, "effective population size")

fstpd1km <-   plotdf(pFST1,   "FST")
fstpd5km <-   plotdf(pFST2,   "FST")
fstpd10km <-  plotdf(pFST3,  "FST")
fstpd15km <-  plotdf(pFST4,  "FST")
fsthfi1km <-  plotdf(hFST1,  "FST")
fsthfi5km <-  plotdf(hFST2,  "FST")
fsthfi10km <- plotdf(hFST3, "FST")
fsthfi15km <- plotdf(hFST4, "FST")
fstro1km <-  plotdf(rFST1,  "FST")
fstro5km <-  plotdf(rFST2,  "FST")
fstro10km <- plotdf(rFST3, "FST")
fstro15km <- plotdf(rFST4, "FST")
fstis1km <-  plotdf(iFST1,  "FST")
fstis5km <-  plotdf(iFST2,  "FST")
fstis10km <- plotdf(iFST3, "FST")
fstis15km <- plotdf(iFST4, "FST")

allModelFramepd <- data.frame(rbind(gdpd1km, gdpd5km, gdpd10km, gdpd15km,
                                    arpd1km, arpd5km, arpd10km, arpd15km,
                                    nepd1km, nepd5km, nepd10km, nepd15km,
                                    fstpd1km, fstpd5km, fstpd10km, fstpd15km))

allModelFramehfi <- data.frame(rbind(gdhfi1km, gdhfi5km, gdhfi10km, gdhfi15km,
                                     arhfi1km, arhfi5km, arhfi10km, arhfi15km,
                                     nehfi1km, nehfi5km, nehfi10km, nehfi15km,
                                     fsthfi1km, fsthfi5km, fsthfi10km, fsthfi15km))

allModelFrameroad <- data.frame(rbind(gdro1km, gdro5km, gdro10km, gdro15km,
                                      arro1km, arro5km, arro10km, arro15km,
                                      nero1km, nero5km, nero10km, nero15km,
                                     fstro1km,fstro5km,fstro10km,fstro15km))

allModelFrameimp <-  data.frame(rbind(gdis1km, gdis5km, gdis10km, gdis15km,
                                      aris1km, aris5km, aris10km, aris15km,
                                      neis1km, neis5km, neis10km, neis15km,
                                      fstis1km,fstis5km,fstis10km,fstis15km))

allModelFramepd$scale <- rep(c("1 km", "5 km", "10 km", "15 km"), 4)
allModelFramepd$scale <- factor(allModelFramepd$scale, levels = c("15 km", "10 km","5 km", "1 km")) # reorder factor levels for plot
allModelFramepd$Response_var <- factor(allModelFramepd$Response_var, levels = c("FST", "effective population size", "gene diversity", "allelic richness")) # reorder factor levels for plot

allModelFramehfi$scale <- rep(c("1 km", "5 km", "10 km", "15 km"), 4)
allModelFramehfi$scale <- factor(allModelFramehfi$scale, levels = c("15 km", "10 km","5 km", "1 km"))
allModelFramehfi$Response_var <- factor(allModelFramehfi$Response_var, levels = c("FST", "effective population size", "gene diversity", "allelic richness")) # reorder factor levels for plot

allModelFrameroad$scale <- rep(c("1 km", "5 km", "10 km", "15 km"), 4)
allModelFrameroad$scale <- factor(allModelFrameroad$scale, levels = c("15 km", "10 km","5 km", "1 km")) 
allModelFrameroad$Response_var <- factor(allModelFrameroad$Response_var, levels = c("FST", "effective population size", "gene diversity", "allelic richness")) # reorder factor levels for plot

allModelFrameimp$scale <- rep(c("1 km", "5 km", "10 km", "15 km"), 4)
allModelFrameimp$scale <- factor(allModelFrameimp$scale, levels = c("15 km", "10 km","5 km", "1 km")) 
allModelFrameimp$Response_var <- factor(allModelFrameimp$Response_var, levels = c("FST", "effective population size", "gene diversity", "allelic richness")) # reorder factor levels for plot


#### S1: Population Density plot ####
zp2 <- ggplot(allModelFramepd, aes(colour = scale)) +
  geom_hline(yintercept=seq(-1, 1.0, 0.25),  
             lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = Response_var, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = Response_var, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFramepd$Response_var))-0.5, 1), 
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", y = "model coefficient (population density)", title = "") +
  theme(text=element_text(family="Lato Black"))

#### S1: Human Footprint Index plot ####
zp3 <- ggplot(allModelFramehfi, aes(colour = scale)) +
  geom_hline(yintercept=seq(-1.5, 1.0, 0.25),  
             lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = Response_var, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = Response_var, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFramehfi$Response_var))-0.5, 1), 
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", y = "model coefficient (HFI)", title = "") +
  theme(text=element_text(family="Lato Black"))

#### S1: Road plot ####
zp4 <- ggplot(allModelFrameroad, aes(colour = scale)) +
  geom_hline(yintercept=seq(-1.25, 1.0, 0.25),  
             lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = Response_var, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = Response_var, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrameroad$Response_var))-0.5, 1), 
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", y = "model coefficient (roads)", title = "") +
  theme(text=element_text(family="Lato Black"))

#### S1: Impervious surface plot ####
zp5 <- ggplot(allModelFrameimp, aes(colour = scale)) +
  geom_hline(yintercept=seq(-2, 4.0, 1),  
             lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = Response_var, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = Response_var, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrameimp$Response_var))-0.5, 1), 
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
    scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", y = "model coefficient (impervious surface cover)", title = "") +
  theme(text=element_text(family="Lato Black"))

#### S1: Combined plot ####
(zp2 + zp3)/(zp4 + zp5) + plot_layout(guides = "collect")

#### Species specific effects (Fig. S2-S5) ####
## Gene diversity models
gdpu <- uGD1 %>%
  spread_draws(b_urban1, r_species[species,]) %>%
  mutate(species_mean = b_urban1 + r_species) 
gdpu <- spp_plotdat(gdpu)  

gdphpd <- pGD3 %>%
  spread_draws(b_scalelogpd_10_kmP1, r_species[species,]) %>%
  mutate(species_mean = b_scalelogpd_10_kmP1 + r_species) 
gdphpd <- spp_plotdat(gdphpd) 

gdphfi <- hGD3 %>%
  spread_draws(b_scalehfi_10_km, r_species[species,]) %>%
  mutate(species_mean = b_scalehfi_10_km + r_species) 
gdphfi <- spp_plotdat(gdphfi) 

gdpimp <- iGD5 %>%
  spread_draws(b_scaleimp_30m, r_species[species,]) %>%
  mutate(species_mean = b_scaleimp_30m + r_species) 
gdpimp <- spp_plotdat(gdpimp)

plotGD1 <- spp_plot(gdpu,   "\u03b2 (urban/rural category)")
plotGD2 <- spp_plot(gdphpd, "\u03b2 (human population density)")
plotGD3 <- spp_plot(gdphfi, "\u03b2 (Human Footprint Index)")
plotGD4 <- spp_plot(gdpimp, "\u03b2 (% impervious surface cover)")

## Allelic richness models
arpu <- uAR1 %>%
  spread_draws(b_urban1, r_species[species,]) %>%
  mutate(species_mean = b_urban1 + r_species) 
arpu <- spp_plotdat(arpu)  

arphpd <- pAR3 %>%
  spread_draws(b_scalelogpd_10_kmP1, r_species[species,]) %>%
  mutate(species_mean = b_scalelogpd_10_kmP1 + r_species) 
arphpd <- spp_plotdat(arphpd) 

arphfi <- hAR3 %>%
  spread_draws(b_scalehfi_10_km, r_species[species,]) %>%
  mutate(species_mean = b_scalehfi_10_km + r_species) 
arphfi <- spp_plotdat(arphfi)

arpimp <- iAR5 %>%
  spread_draws(b_scaleimp_30m, r_species[species,]) %>%
  mutate(species_mean = b_scaleimp_30m + r_species) 
arpimp <- spp_plotdat(arpimp) 

plotAR1 <- spp_plot(arpu,   "\u03b2 (urban/rural category)")
plotAR2 <- spp_plot(arphpd, "\u03b2 (human population density)")
plotAR3 <- spp_plot(arphfi, "\u03b2 (Human Footprint Index)")
plotAR4 <- spp_plot(arpimp, "\u03b2 (% impervious surface cover)")

## Ne models
nepu <- uNE1 %>%
  spread_draws(b_urban1, r_species[species,]) %>%
  mutate(species_mean = b_urban1 + r_species) 
nepu <- spp_plotdat(nepu)  

nephpd <- pNE3 %>%
  spread_draws(b_scalelogpd_10_kmP1, r_species[species,]) %>%
  mutate(species_mean = b_scalelogpd_10_kmP1 + r_species) 
nephpd <- spp_plotdat(nephpd) 

nephfi <- hNE3 %>%
  spread_draws(b_scalehfi_10_km, r_species[species,]) %>%
  mutate(species_mean = b_scalehfi_10_km + r_species) 
nephfi <- spp_plotdat(nephfi)

plotNE1 <- spp_plot(nepu,   "\u03b2 (urban/rural category)")
plotNE2 <- spp_plot(nephpd, "\u03b2 (human population density)")
plotNE3 <- spp_plot(nephfi, "\u03b2 (Human Footprint Index)")
# no random slope for Ne 30m impervious surface

## FST
fstpu <- uFST1 %>%
  spread_draws(b_urban1, r_species[species,]) %>%
  mutate(species_mean = b_urban1 + r_species) 
fstpu <- spp_plotdat(fstpu)  

fstphpd <- pFST3 %>%
  spread_draws(b_scalelogpd_10_kmP1, r_species[species,]) %>%
  mutate(species_mean = b_scalelogpd_10_kmP1 + r_species) 
fstphpd <- spp_plotdat(fstphpd) 

fstphfi <- hFST3 %>%
  spread_draws(b_scalehfi_10_km, r_species[species,]) %>%
  mutate(species_mean = b_scalehfi_10_km + r_species) 
fstphfi <- spp_plotdat(fstphfi)

fstpimp <- iFST5 %>%
  spread_draws(b_scaleimp_30m, r_species[species,]) %>%
  mutate(species_mean = b_scaleimp_30m + r_species) 
fstpimp <- spp_plotdat(fstpimp)

plotFST1 <- spp_plot(fstpu,   "\u03b2 (urban/rural category)")
plotFST2 <- spp_plot(fstphpd, "\u03b2 (human population density)")
plotFST3 <- spp_plot(fstphfi, "\u03b2 (Human Footprint Index)")
plotFST4 <- spp_plot(fstpimp, "\u03b2 (% impervious surface cover)")

# Gene diversity plot
(plotGD1 + plotGD2) /(plotGD3 + plotGD4) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = 'Species specific effects',
                  subtitle = 'Gene diversity') & theme(text=element_text("Lato Black"))

# Allelic richness plot
(plotAR1 + plotAR2)/(plotAR3 + plotAR4) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = 'Species specific effects',
                  subtitle = 'Allelic richness') & theme(text=element_text("Lato Black"))
# FST plot
(plotFST1 + plotFST2)/(plotFST3 + plotFST4) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = 'Species specific effects',
                  subtitle = expression(F[ST])) & theme(text=element_text("Lato Black"))
# Ne plot
ly <- '
AB
C#
'
plotNE1 + plotNE2 + plotNE3 +
  plot_layout(guides = "collect", design = ly) +
  plot_annotation(title = 'Species specific effects',
                  subtitle = "Effective population size") & theme(text=element_text("Lato Black"))