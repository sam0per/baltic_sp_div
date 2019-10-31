rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo", "extrafont", "ggthemes", "bbmle", "svglite")
.packagesdev = "thomasp85/patchwork"
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

# carl = read.csv("data/cline data template.csv", sep = ";")
carl = read.csv("Baltic_divergence/data/Baltic_clines.csv")
carl_sa = read.csv("Baltic_divergence/data/Baltic_clines_sal.csv")
head(carl)
with(data = carl, plot(km, rel_fst, col=species, pch=19))

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

theta.init <- list(mytilus=list(centre=50,w=150,left=0.1,right=0.8,sl=0.1,sc=0.1,sr=0.1),
                   cod=list(centre=400,w=150,left=0.005,right=0.9,sl=0.05,sc=0.1,sr=0.1))

# mle2(cline, theta.init, data=list(position=carl[carl$species=="cod", ]$km,
#                                   phen=carl[carl$species=="cod", ]$rel_fst),
#      control=list(parscale=abs(unlist(theta.init))))

mle.cline.m <- lapply(as.character(unique(carl$species)), function(sp) {
  mle2(cline, theta.init[[sp]], data=list(position=carl[carl$species==sp, ]$km,
                                    phen=carl[carl$species==sp, ]$rel_fst),
       control=list(parscale=abs(unlist(theta.init[[sp]]))))
})
lapply(mle.cline.m, summary)
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
# range(carl$km)
# cline_fit = cline_pl(phen = carl$FST, position = carl$km, centre = 24.5, w = 54.9, crab = 0.02, wave = 0.61,
#                      sc = 0.01, sh = 0.04, sw = 0.07)
# cline_fit = cline_pl(position = seq(from = -1440, to = 500, by = 10), centre = 24.5, w = 54.9, crab = 0.02, wave = 0.61,
#                      sc = 0.01, sh = 0.04, sw = 0.07)
# cline_pl(position = 25, centre = 24.5, w = 54.9, crab = 0.02, wave = 0.61,
#          sc = 0.01, sh = 0.04, sw = 0.07)
tar_sp = as.character(unique(carl$species))
cline_fit = lapply(seq_along(mle.cline.coef), function(x) {
  cline_pl(position = seq(from = min(carl[carl$species==tar_sp[x], ]$km),
                          to = max(carl[carl$species==tar_sp[x], ]$km), by = 10), centre = mle.cline.coef[[x]]['centre'],
           w = mle.cline.coef[[x]]['w'], left = mle.cline.coef[[x]]['left'], right = mle.cline.coef[[x]]['right'],
           sl = mle.cline.coef[[x]]['sl'], sc = mle.cline.coef[[x]]['sc'], sr = mle.cline.coef[[x]]['sr'])
})

lines(x = cline_fit[[1]]$position, y = cline_fit[[1]]$phen_cline, col='red', lwd=3)
lines(x = cline_fit[[2]]$position, y = cline_fit[[2]]$phen_cline, col='black', lwd=3)
abline(v = mle.cline.coef[[1]]['centre'], col='red', lty=2)
abline(v = mle.cline.coef[[2]]['centre'], col='black', lty=2)

cline_fit_sp = rbindlist(lapply(seq_along(tar_sp), function(y) {
  mutate(cline_fit[[y]], species = tar_sp[y])
}))
cline_coef_sp = lapply(seq_along(tar_sp), function(y) {
  df_tmp = data.frame(pars = names(mle.cline.coef[[y]]), val = mle.cline.coef[[y]])
  mutate(df_tmp, species = tar_sp[y])
})
clines_img = ggplot(data = cline_fit_sp, aes(col=species)) +
  geom_vline(xintercept = cline_coef_sp[[1]]$val[cline_coef_sp[[1]]$pars == 'centre'], col='blue', size = 1,
             linetype = 'dashed', alpha = 0.7) +
  geom_vline(xintercept = cline_coef_sp[[2]]$val[cline_coef_sp[[2]]$pars == 'centre'], col='red', size = 1,
             linetype = 'dashed', alpha = 0.7) +
  scale_color_manual(values = c('red', 'blue')) +
  geom_point(data = carl, aes(x = km, y = rel_fst), size = 2) +
  geom_line(aes(x = position, y = phen_cline), size = 1.2, alpha = 0.7) +
  theme(legend.position="top")

dir.create(path = "Baltic_divergence/figures")
dir.create(path = "Baltic_divergence/tables")
ggsave(file=paste0("Baltic_divergence/figures/baltic_div_", max(as.integer(factor(tar_sp))), "species.svg"), plot=clines_img, width=10, height=8)
write.csv(x = cline_fit_sp,
          file = paste0("Baltic_divergence/figures/baltic_div_", max(as.integer(factor(tar_sp))), "species.csv"),
          row.names = FALSE)

# plot(x = cline_fit$position, y = cline_fit$phen_cline, type = "l")
# with(data = carl, points(km, FST))
