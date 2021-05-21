## 3. Phylogenetic regressions ##
library(taxize) # v0.9.97
library(ape) # v5.4
library(brms) # v2.13.0


#### Load data ####
amph <- read.csv("amph_data.csv", head = TRUE)

amph$urban <- as.factor(amph$urban)
amph$species <- as.factor(amph$species)

# Ne subsets
amph.ne <- amph[!(is.na(amph$Ne)),]

# Fst subsets
amph.fst <- amph[!is.na(amph$global_fst),]

### Phylogenetic tree ####

# 1) make a tree using species names
amph$species <- gsub("_", " ", amph$species) # remove underscore to match names
aspecies <- unique(amph$species)
species <- as.list(aspecies) # list of species names

amph.ne$species <- gsub("_", " ", amph.ne$species) 
aspeciesne <- unique(amph.ne$species)
nspecies <- as.list(aspeciesne)

amph.fst$species <- gsub("_", " ", amph.fst$species) 
aspeciesf <- unique(amph.fst$species)
fspecies <- as.list(aspeciesf)

# get NCBI unique identifier (UID) for each species:
uids <- lapply(species, function(x) get_uid(x, messages = FALSE)[1])

# create a tree
taxize_class <- classification(uids, db = "ncbi")
taxize_tree <- class2tree(taxize_class, check = TRUE)
plot(taxize_tree)
amphtree <- taxize_tree$phylo # class "phylo" object
amphtree.nef <- drop.tip(amphtree, "Plethodon cinereus") ## Drop species for Ne & FST trees

#### phylogenetic regression with brms (main text models only) ####

A <- ape::vcv.phylo(amphtree) # create phylo variance covariance matrix
A <- A/det(A)^(1/length(aspecies)) # standardize the variance covariance matrix

amph.db$species <- gsub("_", " ", amph.db$species) # amph.db is data + MEMs from urban_analysis script
amph.db$phylo <- amph.db$species # create phylo column separate from species column

amph.db.ne$species <- gsub("_", " ", amph.db.ne$species)
amph.db.ne$phylo <- amph.db.ne$species

amph.db.fst$species <- gsub("_", " ", amph.db.fst$species)
amph.db.fst$phylo <- amph.db.fst$species

brmfun_phy <- function(model_spec, data_group){
  brm(model_spec, cores = 4, chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.999, max_treedepth = 15),
      data = data_group, data2 = list(A = A))
}

# gene diversity
a_gd_1p <- bf(scale(gene_diversity) ~ urban + 
               MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (urban|species) + (1|gr(phylo, cov = A)))

a_gd_2p <- bf(scale(gene_diversity) ~ scale(log(pd_10_km+1)) + 
               MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(log(pd_10_km+1))|species) + (1|gr(phylo, cov = A)))

a_gd_3p <- bf(scale(gene_diversity) ~ scale(hfi_10_km) + 
                MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(hfi_10_km)|species) + (1|gr(phylo, cov = A)))

a_gd_4p <- bf(scale(gene_diversity) ~ MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (1|species) + (1|gr(phylo, cov = A)))

IS_gd_30mp <- bf(scale(gene_diversity) ~ scale(imp_30m) +
                  MEM2 + MEM3 + MEM4 + MEM7 + MEM10 + (scale(imp_30m)|species) + (1|gr(phylo, cov = A)))

AGD1p <-brmfun_phy(a_gd_1p, amph.db)
AGD2p <-brmfun_phy(a_gd_2p, amph.db)
AGD3p <-brmfun_phy(a_gd_3p, amph.db)
AGD4p <-brmfun_phy(a_gd_4p, amph.db)
AGD7p <-brmfun_phy(IS_gd_30mp, amph.db)

# allelic richness
a_ar_1p <- bf(scale(allelic_richness) ~ urban + 
               MEM2 + MEM3 + MEM4 + MEM7 + (urban|species) + (1|gr(phylo, cov = A)))

a_ar_2p <- bf(scale(allelic_richness) ~ scale(log(pd_10_km+1)) + 
               MEM2 + MEM3 + MEM4 + MEM7 + (scale(log(pd_10_km+1))|species) + (1|gr(phylo, cov = A)))

a_ar_3p <- bf(scale(allelic_richness) ~ scale(hfi_10_km) + 
               MEM2 + MEM3 + MEM4 + MEM7 + (scale(hfi_10_km)|species) + (1|gr(phylo, cov = A)))

a_ar_4p <- bf(scale(allelic_richness) ~ MEM2 + MEM3 + MEM4 + MEM7 + (1|species) + (1|gr(phylo, cov = A)))

IS_ar_30mp <- bf(scale(allelic_richness) ~ scale(imp_30m) + MEM2 + MEM3 + MEM4 + MEM7 + (scale(imp_30m)|species) + (1|gr(phylo, cov = A)))

AAR1p <-brmfun_phy(a_ar_1p, amph.db)
AAR2p <-brmfun_phy(a_ar_2p, amph.db) 
AAR3p <-brmfun_phy(a_ar_3p, amph.db)
AAR4p <-brmfun_phy(a_ar_4p, amph.db)
AAR7p <-brmfun_phy(IS_ar_30mp, amph.db)

## Reset matrix A for Ne & FST (different tree)
A <- ape::vcv.phylo(amphtree.nef)
A <- A/det(A)^(1/length(aspeciesne))

# effective population size
a_ne_1p <- bf(scale(log(Ne)) ~ urban + (urban|species) + (1|gr(phylo, cov = A)))

a_ne_2p <- bf(scale(log(Ne)) ~ scale(log(pd_10_km+1)) + (scale(log(pd_10_km+1))|species) + (1|gr(phylo, cov = A)))

a_ne_3p <- bf(scale(log(Ne)) ~ scale(hfi_10_km) + (scale(hfi_10_km)|species) + (1|gr(phylo, cov = A)))

a_ne_4p <- bf(scale(log(Ne)) ~ (1|species) + (1|gr(phylo, cov = A)))

IS_ne_30mp <- bf(scale(log(Ne)) ~ scale(imp_30m) + (1|species) + (1|gr(phylo, cov = A)))

# Note models that didn't converge were run with longer chains
ANE1p <- brm(a_ne_1p, cores = 4, chains = 4, iter = 6000, warmup = 1000, 
             control = list(adapt_delta = 0.9999, max_treedepth = 15),
             data = amph.db.ne, data2 = list(A = A))
ANE2p <-brmfun_phy(a_ne_2p, amph.db.ne)
ANE3p <- brm(a_ne_3p, cores = 4, chains = 4, iter = 6000, warmup = 1000, 
             control = list(adapt_delta = 0.9999, max_treedepth = 15),
             data = amph.db.ne, data2 = list(A = A))
ANE4p <-brmfun_phy(a_ne_4p, amph.db.ne)
ANE7p <-brmfun_phy(IS_ne_30mp, amph.db.ne)


# FST
a_fst_1p <- bf(scale(global_fst) ~ urban + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (urban|species) + (1|gr(phylo, cov = A)))

a_fst_2p <- bf(scale(global_fst) ~ scale(log(pd_10_km+1)) + 
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(log(pd_10_km+1))|species) + (1|gr(phylo, cov = A)))

a_fst_3p <- bf(scale(global_fst) ~ scale(hfi_10_km) +
                MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(hfi_10_km)|species) + (1|gr(phylo, cov = A)))

a_fst_4p <- bf(scale(global_fst) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (1|species) + (1|gr(phylo, cov = A)))

IS_fst_30m <- bf(scale(global_fst) ~ scale(imp_30m) + 
                 MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11 + (scale(imp_30m)|species) + (1|gr(phylo, cov = A)))

AFST1p <-brmfun_phy(a_fst_1p, amph.db.fst)
AFST2p <-brmfun_phy(a_fst_2p, amph.db.fst)
AFST3p <-brmfun_phy(a_fst_3p, amph.db.fst)
AFST4p <-brmfun_phy(a_fst_4p, amph.db.fst)
AFST7p <-brmfun_phy(IS_fst_30m, amph.db.fst)

#### Bayes R2 ####
library(performance)
lapply(list(AGD1p, AGD2p, AGD3p, AGD4p, AGD7p,
            AAR1p, AAR2p, AAR3p, AAR4p, AAR7p,
            ANE1p, ANE2p, ANE3p, ANE4p, ANE7p,
            AFST1p, AFST2p, AFST3p, AFST4p, AFST7p), r2_bayes)

#### Plots ####
plotdf <- function(x, response){data.frame(Variable = rownames(summary(x)$fixed)[2],
                                           Coefficient = summary(x)$fixed[2,1],
                                           CIlo95 = summary(x)$fixed[2,3],
                                           CIup95 = summary(x)$fixed[2,4],
                                           CIlo90 = posterior_interval(x, pars = get_variables(x)[2], prob = 0.90)[,1],
                                           CIup90 = posterior_interval(x, pars = get_variables(x)[2], prob = 0.90)[,2],
                                           Response_var = response)
}

#### Coefficient plot ####
m1arp <- plotdf(AAR1p, "allelic richness")
m2arp <- plotdf(AAR2p, "allelic richness")
m3arp <- plotdf(AAR3p, "allelic richness")
m4arp <- plotdf(AAR7p, "allelic richness")
m1gdp <- plotdf(AGD1p, "gene diversity")
m2gdp <- plotdf(AGD2p, "gene diversity")
m3gdp <- plotdf(AGD3p, "gene diversity")
m4gdp <- plotdf(AGD7p, "gene diversity")
m1nep <- plotdf(ANE1p, "effective population size")
m2nep <- plotdf(ANE2p, "effective population size")
m3nep <- plotdf(ANE3p, "effective population size")
m4nep <- plotdf(ANE7p, "effective population size")
m1fstp <- plotdf(AFST1p, "FST")
m2fstp <- plotdf(AFST2p, "FST")
m3fstp <- plotdf(AFST3p, "FST")
m4fstp <- plotdf(AFST7p, "FST")

allModelFrame.amphp <- data.frame(rbind(m1arp, m2arp, m3arp, m4arp,
                                       m1gdp, m2gdp, m3gdp, m4gdp,
                                       m1nep, m2nep, m3nep, m4nep,
                                       m1fstp, m2fstp, m3fstp, m4fstp))

allModelFrame.amphp$Variablecol <- rep(c("urban/rural", "human population density", 
                                        "Human Footprint Index", "% impervious surface*"), 4)

allModelFrame.amphp$Response_var <- factor(allModelFrame.amphp$Response_var, 
                                          levels = c("allelic richness", "gene diversity", 
                                                     "effective population size", "FST"))

zp2 <- ggplot(allModelFrame.amphp, aes(colour = Response_var))

zp2 + geom_hline(yintercept=seq(-1.5, 1.5, 0.5),  # x axis lines
                 lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = Variablecol, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = Variablecol, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame.amphp$Variablecol))-0.5, 1), # y axis btwn groups
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", title = "Model coefficients (phylogenetic correction)") +
  theme(text=element_text(family="Lato Black"))
