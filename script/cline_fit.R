rm(list = ls())
.packages = c("ggplot2", "tibble", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo", "extrafont", "ggthemes", "bbmle", "svglite", "gdata", "RColorBrewer", "segmented", "dplyr")
.packagesdev = "thomasp85/patchwork"
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

vers = 9

spp = 14
# setwd("Documents/Baltic/")
carl_sa = read.xls(list.files(path = getwd(), pattern = "edit", full.names = TRUE), sheet = spp+1, header = TRUE)

carl = rbindlist(lapply(1:spp, function(x) {
  read.xls(list.files(path = getwd(), pattern = "edit", full.names = TRUE), sheet = x, header = TRUE)[, 1:6]
}))
data.frame(table(as.character(carl$species)))
# carl[30, ]
carl = carl[rowSums(is.na(carl)) == 0, ]
head(carl)

carl <- carl[!grepl(pattern = "wrasse", x = carl$species), ]
spp <- nrow(data.frame(table(as.character(carl$species))))

colnames(carl) = c("loc", "wrong_km", "Fst", "km", "rel_fst", "species")

(tar_sp = as.character(unique(carl$species)))
# carl[carl$species=='Ballan_wrasse', 'km'] <- carl[carl$species=='Ballan_wrasse', 'km'] - 2000
lapply(1:spp, function(x) {
  ggplot(data = carl[carl$species==tar_sp[x], ]) +
    geom_point(aes(x = km, y = rel_fst, col=species), size=3) +
    # scale_colour_brewer(palette = "Set3") +
    # scale_color_viridis_d(option = "D") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank())
})

# carl = read.csv("data/cline data template.csv", sep = ";")
# carl = read.csv("baltic_sp_div/data/Baltic_clines.csv")
# carl_sa = read.csv("baltic_sp_div/data/Baltic_clines_sal.csv")
# carl = carl[!carl$species == "Flounder",]
# carl = carl[!carl$species == "Macoma",]
# carl = carl[!carl$species == "Cod",]

sample_n(carl, size = 10)
# brewer.pal.info
ggplot(data = carl) +
  geom_point(aes(x = km, y = rel_fst, col=species), size=3) +
  # scale_color_manual(values = c("#a6cee3", "#c51b7d", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  #                               "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#b15928", "#000000", "#969696")) +
  scale_color_manual(values = c("#c51b7d", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                                "#ff7f00", "#cab2d6", "#6a3d9a", "#b15928", "#000000", "#969696")) +
  # scale_colour_brewer(palette = "Set3") +
  # scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(col = guide_legend(nrow = 2))

skel_in = ggplot(data = carl[carl$species=="Skeletonema",]) +
  geom_point(aes(x = km, y = rel_fst, col=species), size=2) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_skeletonema_data.svg"),
       plot=skel_in, width=12, height=8)
idot_in = ggplot(data = carl[carl$species=="Idotea",]) +
  geom_point(aes(x = km, y = rel_fst, col=species), size=2) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_idotea_data.svg"),
       plot=idot_in, width=12, height=8)
maco_in = ggplot(data = carl[carl$species=="Macoma",]) +
  geom_point(aes(x = km, y = rel_fst, col=species), size=2) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_macoma_data.svg"),
       plot=idot_in, width=12, height=8)
cod_in = ggplot(data = carl[carl$species=="Cod",]) +
  geom_point(aes(x = km, y = rel_fst, col=species), size=2) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_cod_data.svg"),
       plot=idot_in, width=12, height=8)
flou_in = ggplot(data = carl[carl$species=="Flounder",]) +
  geom_point(aes(x = km, y = rel_fst, col=species), size=2) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_flounder_data.svg"),
       plot=idot_in, width=12, height=8)

cline <- function(phen,position,centre,w,left,right,sl,sc,sr){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- left + (right-left)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling
  
  s_x <- sqrt(sl^2 + 4*p_x*(1-p_x)*sc^2 + (p_x^2)*(sr^2-sl^2))
  # unimodal variance model as in HZAR, sx is SD for wave/hybrid/crab (assumes variances are additive, unlike Cfit)
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

# carl = carl[!carl$species == "Plaice",]
# carl = carl[!carl$species == "Dab",]
# carl = carl[!carl$species == "Small_Sandeel",]
(tar_sp = as.character(unique(carl$species)))

fvalc = 313
fvalw = 21
dvalc = -67
dvalw = 330
ssvalc = -180
ssvalw = 157
gsvalc = -92
gsvalw = 276
pvalc = -164
pvalw = 268
skvalc = 235
skvalw = 418

# pvalc = runif(n = 1, min = -200, max = 0)
# pvalw = runif(n = 1, min = 0, max = 300)
# dvalc = runif(n = 1, min = -200, max = 0)
# dvalw = runif(n = 1, min = 0, max = 400)
# ssvalc = runif(n = 1, min = -200, max = 200)
# ssvalw = runif(n = 1, min = 0, max = 500)
# gsvalc = runif(n = 1, min = -300, max = -50)
# gsvalw = runif(n = 1, min = 0, max = 500)
# skvalc = runif(n = 1, min = 200, max = 500)
# skvalw = runif(n = 1, min = 0, max = 500)

theta.init <- list(Corkwing_wrasse=list(centre=-700,w=350,left=0.05,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   Ballan_wrasse=list(centre=-900,w=200,left=0.05,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   Atlantic_cod=list(centre=-90,w=100,left=0.1,right=0.9,sl=0.1,sc=0.1,sr=0.1),
                   Herring=list(centre=400,w=400,left=0.05,right=0.95,sl=0.1,sc=0.1,sr=0.1),
                   Turbot=list(centre=300,w=550,left=0.05,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   Flounder=list(centre=fvalc,w=fvalw,left=0.1,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   Mytilus=list(centre=50,w=150,left=0.05,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   Limecola=list(centre=150,w=100,left=0.2,right=0.85,sl=0.2,sc=0.1,sr=0.2),
                   Idotea=list(centre=210,w=490,left=0.3,right=0.9,sl=0.1,sc=0.1,sr=0.1),
                   Europ_plaice=list(centre=pvalc,w=pvalw,left=0.1,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   Dab=list(centre=dvalc,w=dvalw,left=0.1,right=0.85,sl=0.1,sc=0.1,sr=0.1),
                   Small_sandeel=list(centre=ssvalc,w=ssvalw,left=0.1,right=0.9,sl=0.1,sc=0.1,sr=0.1),
                   Greater_sandeel=list(centre=gsvalc,w=gsvalw,left=0.1,right=0.85,sl=0.1,sc=0.1,sr=0.1),
                   Skeletonema=list(centre=skvalc,w=skvalw,left=0.1,right=0.9,sl=0.1,sc=0.1,sr=0.1))

# theta.init <- list(Atlantic_cod=list(centre=-90,w=100,left=0.1,right=0.9,sl=0.1,sc=0.1,sr=0.1),
#                    Herring=list(centre=400,w=400,left=0.05,right=0.95,sl=0.1,sc=0.1,sr=0.1),
#                    Turbot=list(centre=300,w=550,left=0.05,right=0.8,sl=0.1,sc=0.1,sr=0.1),
#                    Flounder=list(centre=fvalc,w=fvalw,left=0.1,right=0.8,sl=0.1,sc=0.1,sr=0.1),
#                    Mytilus=list(centre=50,w=150,left=0.05,right=0.8,sl=0.1,sc=0.1,sr=0.1),
#                    Limecola=list(centre=150,w=100,left=0.2,right=0.85,sl=0.2,sc=0.1,sr=0.2),
#                    Idotea=list(centre=210,w=490,left=0.3,right=0.9,sl=0.1,sc=0.1,sr=0.1),
#                    Europ_plaice=list(centre=pvalc,w=pvalw,left=0.1,right=0.8,sl=0.1,sc=0.1,sr=0.1),
#                    Dab=list(centre=dvalc,w=dvalw,left=0.1,right=0.85,sl=0.1,sc=0.1,sr=0.1),
#                    Small_sandeel=list(centre=ssvalc,w=ssvalw,left=0.1,right=0.9,sl=0.1,sc=0.1,sr=0.1),
#                    Greater_sandeel=list(centre=gsvalc,w=gsvalw,left=0.1,right=0.85,sl=0.1,sc=0.1,sr=0.1),
#                    Skeletonema=list(centre=skvalc,w=skvalw,left=0.1,right=0.9,sl=0.1,sc=0.1,sr=0.1))

write.csv(x = data.frame(val = unlist(theta.init)),
          file = paste0("baltic_sp_div/tables/baltic_div_", max(as.integer(factor(tar_sp))), "species_init_val_v",
                                                                  vers, ".csv"), row.names = TRUE)
# mle2(cline, theta.init, data=list(position=carl[carl$species=="cod", ]$km,
#                                   phen=carl[carl$species=="cod", ]$rel_fst),
#      control=list(parscale=abs(unlist(theta.init))))

mle.cline.m <- lapply(as.character(unique(carl$species)), function(sp) {
  mle2(cline, theta.init[[sp]], data=list(position=carl[carl$species==sp, ]$km,
                                    phen=carl[carl$species==sp, ]$rel_fst),
       control=list(parscale=abs(unlist(theta.init[[sp]]))))
})
lapply(mle.cline.m, summary)
# lapply(mle.cline.m, AIC)

mle.cline.coef = lapply(mle.cline.m, function(x) round(coef(x), digits = 3))

cline_pl <- function(position,centre,w,left,right,sl,sc,sr){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- left + (right-left)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling
  
  s_x <- sqrt(sl^2 + 4*p_x*(1-p_x)*sc^2 + (p_x^2)*(sr^2-sl^2))
  # unimodal variance model as in HZAR, sx is SD for wave/hybrid/crab (assumes variances are additive, unlike Cfit)
  
  # minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  # return(minusll)
  
  phen_cline = data.frame(phen_cline = z_x, sd_cline = s_x, position = position)
  return(phen_cline)
}
#############
# TOY CLINE #
#############
library(ggplot2)
pd <- cline_pl(position = 1:100, centre = 50, w = 25, left = 0, right = 1, sl = 0.2, sc = 0.2, sr = 0.2)

pd$phen_cline <- pd$phen_cline+round(runif(n = 100, min = 0, max = 0.1), 3)
max(pd$phen_cline)
pd$phen_cline <- ifelse(test = pd$phen_cline > 1, yes = pd$phen_cline - 0.080794, no = pd$phen_cline)

pd$position <- pd$position+round(runif(n = 100, min = 0, max = 10), 3)
pp <- cline_pl(position = 1:100, centre = 53, w = 25, left = 0.03, right = 0.975, sl = 0.2, sc = 0.2, sr = 0.2)

xrow <- sample(x = 1:100, size = 50, replace = FALSE)
clp <- ggplot(data = pd[xrow, ]) +
  geom_vline(xintercept = 48, linetype = "dashed", size = 1.5, col = "#1f77b4") +
  geom_point(aes(x = position, y = phen_cline), size = 2.5) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.title = element_blank())
clp
ggsave(filename = "/Volumes/Seagate Remote Backup/phd/thesis/defence/figures/toy_freq_gradient.pdf",
       plot = clp, width = 6, height = 4, dpi = "screen")

clp <- ggplot(data = pd[xrow, ]) +
  geom_vline(xintercept = 53, linetype = "dashed", size = 1.5, col = "#ff7f0e") +
  geom_vline(xintercept = 48, linetype = "dashed", size = 1.5, col = "#1f77b4") +
  geom_point(aes(x = position, y = phen_cline), size = 2.5) +
  # geom_line(aes(x = position, y = phen_cline2), size = 3, col = "grey") +
  # geom_line(data = pp, aes(x = position, y = phen_cline3), size = 3, col = "grey") +
  geom_line(data = pp, aes(x = position, y = phen_cline), size = 2, col = "#ff7f0e") +
  # geom_line(data = pp, aes(x = position, y = phen_cline4), size = 3, col = "blue") +
  # xlim(c(1, 104)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.title = element_blank())
clp
ggsave(filename = "/Volumes/Seagate Remote Backup/phd/thesis/defence/figures/toy_cline_pl_centre.pdf",
       plot = clp, width = 6, height = 4, dpi = "screen")

pp <- cline_pl(position = 1:100, centre = 50, w = 10, left = 0.1, right = 0.8, sl = 0.2, sc = 0.2, sr = 0.2)
pp$phen_cline2 <- cline_pl(position = 1:100, centre = 53, w = 15, left = 0.15, right = 0.9, sl = 0.2, sc = 0.2, sr = 0.2)[,1]
# pp$phen_cline3 <- cline_pl(position = 1:100, centre = 48, w = 5, left = 0.1, right = 0.3, sl = 0.2, sc = 0.2, sr = 0.2)[,1]
# pp$phen_cline4 <- cline_pl(position = 1:100, centre = 49, w = 15, left = 0.3, right = 0.7, sl = 0.2, sc = 0.2, sr = 0.2)[,1]

clp <- ggplot(data = pp) +
  geom_vline(xintercept = 50, linetype = "dashed", size = 1.5, col = "#31a354") +
  geom_vline(xintercept = 53, linetype = "dashed", size = 1.5, col = "#636363") +
  # geom_point(aes(x = position, y = phen_cline), size = 4) +
  geom_line(data = pp, aes(x = position, y = phen_cline2), size = 2, col = "#636363") +
  # geom_line(data = pp, aes(x = position, y = phen_cline3), size = 3, col = "grey") +
  geom_line(data = pp, aes(x = position, y = phen_cline), size = 2, col = "#31a354") +
  # geom_line(data = pp, aes(x = position, y = phen_cline4), size = 3, col = "blue") +
  # xlim(c(1, 104)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.title = element_blank())
clp
ggsave(filename = "/Volumes/Seagate Remote Backup/phd/thesis/defence/figures/toy_cline_chapter2.pdf",
       plot = clp, width = 6, height = 4, dpi = "screen")
# ggsave(filename = "baltic_sp_div/figures/toy_cline_pl.pdf", plot = clp, scale = 0.7, dpi = "screen")
# ggsave(filename = "Documents/Baltic/baltic_sp_div/figures/toy_cline_pl_centre.pdf", plot = clp, scale = 0.7, dpi = "screen")
#############
#############
# range(carl$km)
# cline_fit = cline_pl(phen = carl$FST, position = carl$km, centre = 24.5, w = 54.9, crab = 0.02, wave = 0.61,
#                      sc = 0.01, sh = 0.04, sw = 0.07)
# cline_fit = cline_pl(position = seq(from = -1440, to = 500, by = 10), centre = 24.5, w = 54.9, crab = 0.02, wave = 0.61,
#                      sc = 0.01, sh = 0.04, sw = 0.07)
# cline_pl(position = 25, centre = 24.5, w = 54.9, crab = 0.02, wave = 0.61,
#          sc = 0.01, sh = 0.04, sw = 0.07)

cline_fit = lapply(seq_along(mle.cline.coef), function(x) {
  cline_pl(position = seq(from = min(carl[carl$species==tar_sp[x], ]$km),
                          to = max(carl[carl$species==tar_sp[x], ]$km), by = 10), centre = mle.cline.coef[[x]]['centre'],
           w = mle.cline.coef[[x]]['w'], left = mle.cline.coef[[x]]['left'], right = mle.cline.coef[[x]]['right'],
           sl = mle.cline.coef[[x]]['sl'], sc = mle.cline.coef[[x]]['sc'], sr = mle.cline.coef[[x]]['sr'])
})

# lines(x = cline_fit[[1]]$position, y = cline_fit[[1]]$phen_cline, col='red', lwd=3)
# lines(x = cline_fit[[2]]$position, y = cline_fit[[2]]$phen_cline, col='black', lwd=3)
# abline(v = mle.cline.coef[[1]]['centre'], col='red', lty=2)
# abline(v = mle.cline.coef[[2]]['centre'], col='black', lty=2)

cline_fit_sp = rbindlist(lapply(seq_along(tar_sp), function(y) {
  mutate(cline_fit[[y]], species = tar_sp[y])
}))
cline_coef_sp = rbindlist(lapply(seq_along(tar_sp), function(y) {
  df_tmp = data.frame(pars = names(mle.cline.coef[[y]]), val = mle.cline.coef[[y]])
  mutate(df_tmp, species = tar_sp[y])
}))

clines_img = ggplot(data = cline_fit_sp, aes(col=species)) +
  geom_vline(data = cline_coef_sp[cline_coef_sp$pars == 'centre',], aes(xintercept = val, col=species),
             size = 1, linetype = 'dashed', alpha = 0.7) +
  # geom_vline(xintercept = cline_coef_sp[[1]]$val[cline_coef_sp[[1]]$pars == 'centre'], col='blue', size = 1,
  #            linetype = 'dashed', alpha = 0.7) +
  # geom_vline(xintercept = cline_coef_sp[[2]]$val[cline_coef_sp[[2]]$pars == 'centre'], col='red', size = 1,
  #            linetype = 'dashed', alpha = 0.7) +
  # scale_color_manual(values = c('red', 'blue')) +
  # scale_color_viridis_d() +
  # scale_color_manual(values = c("#a6cee3", "#c51b7d", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  #                               "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#b15928", "#000000", "#969696")) +
  scale_color_manual(values = c("#c51b7d", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                                "#ff7f00", "#cab2d6", "#6a3d9a", "#b15928", "#000000", "#969696")) +
  # scale_colour_brewer(palette = "Paired") +
  geom_point(data = carl, aes(x = km, y = rel_fst), size = 2) +
  geom_line(aes(x = position, y = phen_cline), size = 1.2, alpha = 0.7) +
  theme_bw() +
  theme(legend.position="top", legend.title = element_blank(), axis.text.x = element_blank(),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  # labs(y='Normalized genetic divergence', x='') +
  labs(y='Genetic divergence', x='') +
  guides(col = guide_legend(nrow = 2))
clines_img

clines_img = ggplot(data = cline_fit_sp, aes(col=species)) +
  # geom_vline(data = cline_coef_sp[cline_coef_sp$pars == 'centre',], aes(xintercept = val, col=species),
  #            size = 1, linetype = 'dashed', alpha = 0.7) +
  # geom_vline(xintercept = cline_coef_sp[[1]]$val[cline_coef_sp[[1]]$pars == 'centre'], col='blue', size = 1,
  #            linetype = 'dashed', alpha = 0.7) +
  # geom_vline(xintercept = cline_coef_sp[[2]]$val[cline_coef_sp[[2]]$pars == 'centre'], col='red', size = 1,
  #            linetype = 'dashed', alpha = 0.7) +
  scale_color_manual(values = c("#a6cee3", "#c51b7d", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                                "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#b15928", "#000000", "#969696")) +
  # scale_color_manual(values = rep("#000000", 14)) +
  # geom_point(data = carl, aes(x = km, y = rel_fst), size = 2) +
  geom_line(aes(x = position, y = phen_cline), size = 1.2, alpha = 0.7) +
  # theme_bw() +
  theme(legend.position="none", legend.title = element_blank(),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.4, linetype = "solid",
                                 colour = "black")) +
  labs(y='', x='')
clines_img

dir.create(path = "baltic_sp_div/figures")
dir.create(path = "baltic_sp_div/tables")
# ggsave(file=paste0("baltic_sp_div/figures/baltic_div_", max(as.integer(factor(tar_sp))), "species_", vers, ".svg"),
#        plot=clines_img, width=10, height=8)
# write.csv(x = cline_fit_sp,
#           file = paste0("baltic_sp_div/tables/baltic_div_", max(as.integer(factor(tar_sp))), "species_", vers, ".csv"),
#           row.names = FALSE)
write.csv(x = cline_coef_sp,
          file = paste0("baltic_sp_div/tables/baltic_div_", max(as.integer(factor(tar_sp))), "species_cline_coef_v",
                        vers, ".csv"),
          row.names = FALSE)

head(carl_sa)
sal_img = ggplot(data = carl_sa) +
  geom_point(aes(x = Dist.fr.Entrance, y = average.salinity)) +
  # geom_text(data = carl_sa[1:2,], aes(x = DistEntrance, y = 30, label = Loc))  +
  # geom_text(data = carl_sa[3:6,], aes(x = DistEntrance, y = 33, label = Loc))  +
  # geom_text(data = carl_sa[7:8,], aes(x = DistEntrance, y = 27, label = Loc))  +
  # geom_text(data = carl_sa[9,], aes(x = DistEntrance, y = 35, label = Loc))  +
  # geom_text(data = carl_sa[10:11,], aes(x = DistEntrance, y = 33, label = Loc))  +
  # geom_text(data = carl_sa[12:16,], aes(x = DistEntrance, y = 27, label = Loc))  +
  # geom_text(data = carl_sa[16,], aes(x = DistEntrance, y = 35, label = Loc))  +
  theme_bw() +
  theme(rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        # panel.background = element_blank(),
        axis.line = element_line(size = 0.4, linetype = "solid",
                                 colour = "black")) +
  xlim(min(cline_fit_sp$position), max(cline_fit_sp$position)) +
  labs(y='Average salinity', x='Distance (km)')
sal_img
# ggsave(file=paste0("baltic_sp_div/figures/baltic_div_", max(as.integer(factor(tar_sp))), "salinity.svg"),
#        plot=sal_img, width=10, height=8, bg = "transparent")

cline_sal_img = clines_img + sal_img + plot_layout(ncol = 1)
cline_sal_img
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_", max(as.integer(factor(tar_sp))), "species_cline_sal_v", vers, ".svg"),
       plot=cline_sal_img, width=12, height=8)
ggsave(file=paste0("baltic_sp_div/figures/baltic_div_", max(as.integer(factor(tar_sp))), "species_cline_sal_v", vers, ".pdf"),
       plot=cline_sal_img, scale = 0.7, dpi = "screen")

# carl$km
# carl_sa$DistEntrance
# setdiff(carl$km, carl_sa$DistEntrance)
# intersect(carl$km, carl_sa$DistEntrance)

# Using the cowplot package
legend <- cowplot::get_legend(clines_img)
library(grid)
library(gridExtra) 
pdf(file = paste0("baltic_sp_div/figures/baltic_div_", max(as.integer(factor(tar_sp))), "species_cline_legends_", vers, ".pdf"),
    width = 12, height = 4)
grid.newpage()
grid.draw(legend)
dev.off()
