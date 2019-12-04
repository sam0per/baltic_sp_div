rm(list = ls())
.packages = c("ggplot2", "dplyr", "tibble", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo", "extrafont", "ggthemes", "bbmle", "svglite", "gdata", "RColorBrewer", "segmented")
.packagesdev = "thomasp85/patchwork"
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

vers = 4

##########################
######## sim data ########
##########################
set.seed(15)
a.sim<-0 # intercept for segment 1
b.sim<-0.02 # slope for segment 1
c.sim<-0.2 # intercept for segment 2
d.sim<-0.01 # slope for segment 2
e.sim<-0.2 # intercept for segment 3
f.sim<-0.03 # slope for segment 3
n<-50
brk1<-10 # breakpoint 1
brk2<-35 # breakpoint 2
x <- c(1:n)
y <- numeric(n)
y[1:brk1]<-a.sim+b.sim*x[1:brk1]+rnorm(brk1, 0, .01)
y[(brk1+1):brk2] <- (c.sim+d.sim*brk1) + rnorm((brk2-brk1), 0, .01) # flat segment
y[brk1]
y[brk1+1]
y[(brk1+1):brk2] = y[(brk1+1):brk2] - (y[brk1+1]-y[brk1])
y[(brk2+1):n] <- e.sim+f.sim*x[(brk2+1):n]+rnorm(n-brk2, 0, .01)
y[(brk2+1):n] = y[(brk2+1):n] - (y[(brk2+1)]-y[brk2])
# y[n]<-y[n]-.30*y[n] #subtract from last point, as these are often outside of the normal pattern
# wt<-c(rep(50, n-4), c(40, 40, 35, 5)) #weight variable 

dat<-as.data.frame(cbind(x, y)) # make dataframe 
# dat.expand <- dat[ rep( seq(dim(dat)[1]), dat$wt),]# second dataset with rows repeated based on weight

with(dat, plot(x,y, ylim=c(0, max(y)), pch=16))### plot
abline(v = c(10,36))
# y[brk2+1]

# mod<-function(par, data){
#   a <- par[1]
#   a <- 0.2
#   b <-par[2]
#   b = 0.4
#   x.star <-par[3]
#   x.star = 3.8
#   data = dat
#   yfit<-c(NA,length(data$y))
#   small.x<-I(data$x<=x.star)
#   yfit[small.x==T]<-(a+b*data$x[small.x==T]) 
#   yfit[small.x!=T]<-(a+b*x.star)
#   return(yfit)
#   sum((data$y-yfit)^2)
# }
# fit1<-optim(par=c(a=.5, b=.5, x.star=3), mod, data=dat, hessian=T)
# abline(a = fit1$par[1], b = fit1$par[2], col = "red")
# mod(par = c(a=0.2365734, b=0.4400691, x.star=3.8877919), data = dat)
# lines(mod(par = c(a=0.2365734, b=0.4400691, x.star=3.8877919), data = dat))

stairs_mod = function(x, y.obs, b0, b1, b4, b5, x1.star, x2.star, sigma) {
  # if (sum(x < 0) > 0) {
  #   idx_x1.star = which(abs(x-x1.star)==min(abs(x-x1.star)))
  #   idx_x2.star = which(abs(x-x2.star)==min(abs(x-x2.star)))
  #   x = 1:length(x)
  #   x1.star = x[idx_x1.star]
  #   x2.star = x[idx_x2.star]
  # }
  # plot(x, y.obs)
  # yfit <- c(NA,length(y.obs))
  yfit <- rep(NA,length(y.obs))
  small.x <- I(x<=x1.star)
  large.x <- I(x>=x2.star)
  mid.x <- I(x<x2.star & x>x1.star)
  
  yfit[small.x==T] <- (b0 + b1 * x[small.x==T])
  # yfit[mid.x==T] <- (b2 + b3 * x1.star)
  yfit[mid.x==T] <- max(yfit[small.x==T])
  yfit[large.x==T] <- (b4 + b5 * x[large.x==T])
  
  # yfit[mid.x==T] <- yfit[mid.x==T] - (yfit[x1.star+1]-yfit[x1.star])
  # sum(is.na(yfit))
  # yfit[mid.x==T] <- rep(yfit[x1.star], length(yfit[mid.x==T]))
  # yfit[large.x==T] <- yfit[large.x==T]-(yfit[x2.star]-yfit[x1.star])
  # lines(yfit-0.1)
  
  minusll <- -sum(dnorm(y.obs, mean = yfit, sd = sigma, log=TRUE))
  if(x1.star > x2.star){minusll <- minusll+1000}
  if(b0 < 0){minusll <- minusll+1000}
  return(minusll)
  # return(yfit)
}
b0<-0       # intercept for segment 1
b1<-0.0005  # slope for segment 1
b2<-0.2     # intercept for segment 2
b3<-0.0004  # slope for segment 2
b4<-0.3     # intercept for segment 3
b5<-0.0005  # slope for segment 3

sim_fit = mle2(minuslogl = stairs_mod, start = list(b0=0.01, b1=0.02, b4=0.02, b5=0.05, x1.star=10, x2.star=36, sigma=0.4),
               data=list(x=dat$x, y.obs=dat$y))
summary(sim_fit)
(sim_coef = round(coef(sim_fit), 3))

stairs_mod_simtest = function(x, y.obs, b0, b1, b4, b5, x1.star, x2.star, sigma) {
  # if (sum(x < 0) > 0) {
  #   idx_x1.star = which(abs(x-x1.star)==min(abs(x-x1.star)))
  #   idx_x2.star = which(abs(x-x2.star)==min(abs(x-x2.star)))
  #   x = 1:length(x)
  #   x1.star = x[idx_x1.star]
  #   x2.star = x[idx_x2.star]
  # }
  # plot(x, y.obs)
  # yfit <- c(NA,length(y.obs))
  yfit <- rep(NA,length(y.obs))
  small.x <- I(x<=x1.star)
  large.x <- I(x>=x2.star)
  mid.x <- I(x<x2.star & x>x1.star)
  
  yfit[small.x==T] <- (b0 + b1 * x[small.x==T])
  # yfit[mid.x==T] <- (b2 + b3 * x1.star)
  yfit[mid.x==T] <- max(yfit[small.x==T])
  yfit[large.x==T] <- (b4 + b5 * x[large.x==T])
  
  # yfit[mid.x==T] <- yfit[mid.x==T] - (yfit[x1.star+1]-yfit[x1.star])
  # sum(is.na(yfit))
  # yfit[mid.x==T] <- rep(yfit[x1.star], length(yfit[mid.x==T]))
  # yfit[large.x==T] <- yfit[large.x==T]-(yfit[x2.star]-yfit[x1.star])
  # lines(yfit-0.1)
  
  # minusll <- -sum(dnorm(y.obs, mean = yfit, sd = sigma, log=TRUE))
  # if(x1.star > x2.star){minusll <- minusll+1000}
  # if(b0 < 0){minusll <- minusll+1000}
  # return(minusll)
  return(yfit)
}
sim_dat = data.frame(x=dat$x,
                     sim_yfit=stairs_mod_simtest(x = dat$x, y.obs = dat$y, b0 = sim_coef[1],
                                                 b1 = sim_coef[2], b4 = sim_coef[3],b5 = sim_coef[4],
                                                 x1.star = sim_coef[5], x2.star = sim_coef[6], sigma = sim_coef[7]))
lines(sim_dat$x, sim_dat$sim_yfit)
# abline(v = 35)
# fit3<-mle2(mod2, start=list(a=0,b=0.5,x.star=brk), data=dat)
# ci<-confint(fit3) 


###########################
####### baltic data #######
###########################
rm(list = ls())
spp = 8
carl = rbindlist(lapply(1:spp, function(x) {
  read.xls(list.files(path = getwd(), pattern = "edit", full.names = TRUE), sheet = x, header = TRUE)[, 1:6]
}))
carl = carl[rowSums(is.na(carl)) == 0, ]
colnames(carl) = c("loc", "wrong_km", "Fst", "km", "rel_fst", "species")
carl_sa = read.xls(list.files(path = getwd(), pattern = "edit", full.names = TRUE), sheet = spp+1, header = TRUE)
(tar_sp = as.character(unique(carl$species)))
lapply(1:spp, function(x) {
  ggplot(data = carl[carl$species==tar_sp[x], ]) +
    geom_point(aes(x = km, y = rel_fst, col=species), size=2) +
    scale_colour_brewer(palette = "Set1") +
    # scale_color_viridis_d(option = "D") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank())
})

############################
##### test baltic data #####
############################
test_dt = carl[carl$species=="Flounder",]
test_dt = carl[carl$species=="Cod",]
plot(test_dt$km, test_dt$rel_fst)
# test_dt$km = 1:nrow(test_dt)

# plot(test_dt$km, test_dt$rel_fst)
# lines(test_dt$km, yfit)
# lm(rel_fst ~ km, data = test_dt)
# test_fit = data.frame(km=test_dt$km, fit=fitted(lm(rel_fst ~ km, data = test_dt)))
# lines(test_fit$km, test_fit$fit)

plot((test_dt$km - test_dt$km[1]), test_dt$rel_fst)
test_fit = mle2(minuslogl = stairs_mod,
                start = list(b0=0.01, b1=0.02, b4=0.02, b5=0.05, x1.star=650, x2.star=900, sigma=0.4),
                data=list(x=(test_dt$km - test_dt$km[1]), y.obs=test_dt$rel_fst))

summary(test_fit)
(test_coef = round(coef(test_fit), 3))
test_dat = data.frame(x=(test_dt$km - test_dt$km[1]),
                      test_yfit=stairs_mod_simtest(x = (test_dt$km - test_dt$km[1]), y.obs = test_dt$rel_fst,
                                                   b0 = test_coef[1], b1 = test_coef[2], b4 = test_coef[3],b5 = test_coef[4],
                                                   x1.star = test_coef[5], x2.star = test_coef[6], sigma = test_coef[7]))
lines(test_dat$x, test_dat$test_yfit)



# lines(x, yfit)
stairs_lm = function(x, y.obs, b0, b1, sigma) {
  yfit <- b0 + b1 * x
  # minusll <- -sum(dnorm(y.obs, mean = yfit, sd = sigma, log=TRUE))
  # return(minusll)
  return(yfit)
}
stairs_lm_mod = mle2(minuslogl = stairs_lm, start = list(b0=0.01, b1=0.02, sigma=0.2),
                     data = list(x=(test_dt$km - test_dt$km[1]), y.obs=test_dt$rel_fst))
summary(stairs_lm_mod)
(lm_coef = round(coef(stairs_lm_mod), 3))
lm_yfit = stairs_lm(x = (test_dt$km - test_dt$km[1]), y.obs = test_dt$rel_fst,
                    b0 = lm_coef[1], b1 = lm_coef[2], sigma = lm_coef[3])
plot((test_dt$km - test_dt$km[1]), test_dt$rel_fst)
lines(x = (test_dt$km - test_dt$km[1]), y = lm_yfit)


###########################
#### segmented package ####
###########################
carl_seg = lapply(tar_sp, function(x) {
  transmute(carl[carl$species==x, ], x = km - min(km), y = rel_fst)
})
lapply(1:spp, function(x) {
  ggplot(data = carl_seg[[x]]) +
    geom_point(aes(x = x, y = y), size=3) +
    scale_colour_brewer(palette = "Paired") +
    # scale_color_viridis_d(option = "D") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank())
})
# dat = data.frame(x=(test_dt$km - test_dt$km[1]), y=test_dt$rel_fst)
# plot(dat$x, dat$y)
out.lm = lapply(1:spp, function(l) {
  lm(y ~ x, data = carl_seg[[l]])
})
carl_brk = list(c(500,800), c(400,800), c(1000,2000), c(1000, 1500), c(600), c(300, 700), c(800, 1300),
                c(400,800), c(250), c(300), c(700,1200), c(450,900))
# out.seg = lapply(1:spp, function(l) {
#   segmented(out.lm[[l]], seg.Z = ~x, psi = list(x = carl_brk[[l]]), npsi = 2,
#             control = seg.control(display = FALSE))
# })
l = 12
with(data = carl_seg[[l]], plot(x = x, y = y))
carl_brk[[l]]
abline(v = c(carl_brk[[l]][1], carl_brk[[l]][2]))
out.seg = segmented(out.lm[[l]], seg.Z = ~x, psi = list(x = carl_brk[[l]]), npsi = 2,
                    control = seg.control(display = FALSE))
cat("For species", tar_sp[l], "was not possible to estimate breakpoints.\n")
summary(out.seg)
logLik(out.seg)
carl_seg[[l]]
plot(out.seg$residuals, out.seg$fitted.values)
hist(out.seg$residuals)
round(residuals(out.seg), 3)
# root mean squared error (RMSE). R calls this quantity the residual standard error
RMSE_seg = sqrt(sum(residuals(out.seg)^2) / df.residual(out.seg))
ci95_seg = data.frame(lci=-2*RMSE_seg, uci=2*RMSE_seg)

dat2 = data.frame(x = carl_seg[[l]]$x, y = broken.line(out.seg)$fit)
ggplot(carl_seg[[l]], aes(x = x, y = y)) +
  geom_point() +
  geom_line(data = dat2, color = 'blue')

AIC(out.lm[[l]])
AIC(out.seg)
# length(mle.cline.m)
AIC(mle.cline.m[[l]])
(dAIC = AIC(mle.cline.m[[l]]) - AIC(out.seg))

seg_np = 7
(segAIC = (2 * seg_np) - (2 * logLik(out.seg)))
cli_np = 7
(cliAIC = (2 * cli_np) - (2 * logLik(mle.cline.m[[l]])))

if (dAIC > 2) {
  cat("For species", tar_sp[l], "segmented model is better than cline model.\n")
} else if (dAIC < -2) {
  cat("For species", tar_sp[l], "cline model is better than segmented model.\n")
} else {
  cat("For species", tar_sp[l], "cline and segmented model do not differ.\n")
}

