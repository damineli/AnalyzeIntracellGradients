################################################################################
# Gradient Analysis Script as presented in Damineli, Portes & Feijó (2020) #####
# "Analyzing intracellular gradients in pollen tubes." #########################
# in Pollen and Pollen Tube Biology: Methods and Protocols: 201-210. ###########
#  ################
# For general usage, substitute information labeles as "USER DEFINED" ##########
################################################################################
# Load required packages, download if absent ###################################
pkgs <- c('fields', 'mixtools', 'tidyverse')

missing <- setdiff(pkgs,
                   rownames(installed.packages()))

if (length(missing) != 0) {
  install.packages(missing, dependencies = TRUE)
}
invisible(lapply(pkgs, require, character.only = TRUE))

# Custom functions used here ###################################################
#-------------------------------------------------------------------------------
AlignByTip <- function(tip.loc, imaj){
  exceed <- which(tip.loc > dim(imaj)[2])
  if(length(exceed) > 0){
    row.limit <- min(which(tip.loc > dim(imaj)[2])) - 1
    message(paste("Kymograph frames", (row.limit + 1), "to", dim(imaj)[2], 
                  "excluded"))
  } else{row.limit <- length(tip.loc)}
  
  imaj.nw <- t(sapply(1:row.limit, GetNewFracRow, imaj = imaj,
                      tip.loc = tip.loc))
  imaj.nw <- imaj.nw[, dim(imaj.nw)[2]:1]
  return(imaj.nw[, 1:ceiling(max(tip.loc))])
}
#-------------------------------------------------------------------------------
GetNewFracRow <- function(row.n, imaj, tip.loc){
  tip.round <- floor(tip.loc[row.n])
  tip.frac <- (tip.loc%%1)[row.n]
  #if(length(2:tip.round) != 0){}
  ind.vec <- 2:tip.round
  zvals.vec <- imaj[row.n, ]
  new.row <- sapply(ind.vec, CalcNewFracVal, frac = tip.frac, vec = zvals.vec)
  val <- c(rep(NA, (dim(imaj)[2]) - length(new.row)), new.row)
  
  return(val)
}
#-------------------------------------------------------------------------------
CalcNewFracVal <- function(i, frac, vec){
  return(weighted.mean(c(vec[i], vec[i - 1]), c(frac, 1 - frac)))
}
#-------------------------------------------------------------------------------

# LOAD DATA: kymograph(s) and tip location series #####################

# Load kymographs obtained in ImageJ / Fiji 
kymo2 <- as.matrix(read.delim("~/Downloads/Kymograph_C1.txt", header=FALSE)) # USER DEFINED: file name / path in ""
kymo1 <- as.matrix(read.delim("~/Downloads/Kymograph2.txt", header=FALSE)) # USER DEFINED: file name / path in ""

# Load time series obtained in feijolab.shinyapps.io/CHUK/
tip.dat <- read.csv("~/Downloads/TipFinding/All_time_series.txt", sep="") # USER DEFINED: file names / path in ""

# If pixel size and time step were not defined in the web app, define below
# otherwise set Both to "1"
pixel.size <- 0.2160014 # USER DEFINED
time.step <- 4 / 60 # USER DEFINED

# Convert growth rate to proper units
tip.dat$growth * pixel.size / time.step

# SECTION 3.2 STEP 2 - INVERT COLUMNS ##########################################
# Employ function defined herein to align kymographs
k1 <- AlignByTip(tip.dat$tip.loc.raw, kymo1)
k2 <- AlignByTip(tip.dat$tip.loc.raw, kymo2)

# Create sequence of spacial and temporal coordinates with real units
xk <- (1:dim(k1)[1] - 1) * time.step
yk <- (1:dim(k1)[2] - 1) * pixel.size

# Plots
pdf("~/Desktop/Fig4A_YFP.pdf", height = 4, width = 8) # USER DEFINED: file names / path in ""
par(mar = par("mar"))
image.plot(x = xk, y = yk, k2, 
           col = tim.colors(256), 
           ylab = expression(paste("Length from tip (", mu, "m)")), 
           xlab = "Time (min)", legend.lab = "YFP fluorescence (AU)", 
           smallplot= c(.865,.885,0.2,0.98), bigplot = c(0.105,0.85,0.2,0.98), useRaster = TRUE, legend.line = 3.6)
dev.off()

pdf("~/Desktop/Fig4B_CFP.pdf", height = 4, width = 8) # USER DEFINED: file names / path in ""
par(mar = par("mar"))
image.plot(x = xk, y = yk, k1, 
           col = tim.colors(256), 
           ylab = expression(paste("Length from tip (", mu, "m)")), 
           xlab = "Time (min)", legend.lab = "CFP fluorescence (AU)", 
           smallplot= c(.865,.885,0.2,0.98), bigplot = c(0.105,0.85,0.2,0.98), useRaster = TRUE, legend.line = 3.6)
dev.off()

pdf("~/Desktop/Fig4C_YFP_CFP.pdf", height = 4, width = 8)
par(mar = par("mar"))
image.plot(x = xk, y = yk, k2/k1, 
           col = tim.colors(256), 
           ylab = expression(paste("Length from tip (", mu, "m)")), 
           xlab = "Time (min)", legend.lab = "Ratiometric fluorescence (YFP/CFP)", 
           smallplot= c(.865,.885,0.2,0.98), bigplot = c(0.105,0.85,0.2,0.98), useRaster = TRUE, legend.line = 3.6)
dev.off()

# SECTION 3.3 STEP 1 & 2 - FILTER TIME POINTS ##################################
# Gaussian mixture model
gr.gmm <- normalmixEM(as.numeric(tip.dat$growth[-1]))

i = sort(gr.gmm$mu, index.return = TRUE)$ix[1]
j = sort(gr.gmm$mu, index.return = TRUE)$ix[2]

# Use intersect to define boundaries between the two distributions
intrsct <- uniroot(function(x) gr.gmm$lambda[i]*dnorm(x, mean = gr.gmm$mu[i], 
                                                      sd = gr.gmm$sigma[i], 
                                                      log = FALSE) - gr.gmm$lambda[j]*dnorm(x, 
                                                                                            mean = gr.gmm$mu[j], 
                                                                                            sd = gr.gmm$sigma[j], 
                                                                                            log = FALSE), 
                   c(gr.gmm$mu[i],gr.gmm$mu[j]))

print(c(gr.gmm$mu[i],gr.gmm$mu[j]))
thrsh <- intrsct$root

# Get indexes to remove from analysis
rmv.ind <- which(tip.dat$growth <= thrsh)

# Fit Gaussian Mixture Model to define growth arrest
pdf("~/Desktop/Fig5_GrowthGMM.pdf", height = 5, width = 6) # USER DEFINED: file names / path in ""
plot(gr.gmm, which = 2, breaks = 50,
     main2 = "Gaussian Mixture Model",
     xlab2 = expression(paste("Instantaneous growth rate (",
                              mu, "m min"^{-1}, ")", sep = "")))
abline(v = thrsh, col = "magenta", lwd = 3)

legend("topright",
       legend = as.expression(sapply(1:length(gr.gmm$mu),
                                     function(i) bquote(mu[.(i)]~'='~.(round(gr.gmm$mu[i], digits = 2))~', '~sigma[.(i)]~'='~.(round(gr.gmm$sigma[i], digits = 2))~', '~lambda[.(i)]~'='~.(round(gr.gmm$lambda[i], digits = 2))))),
       text.col = 2:(length(gr.gmm$mu)+1))
dev.off()

# Filter points
k <- (k2/k1)[, 1:min(which(is.na(k1[1,]))) - 1]
dat <- as.data.frame(k[-rmv.ind, ])
colnames(dat) <- ((1:dim(dat)[2]) - 1) * pixel.size

pdf("~/Desktop/Fig4D_YFP_CFP_filt.pdf", height = 4, width = 8) # USER DEFINED: file names / path in ""
par(mar = par("mar"))
image.plot(x = ((1:dim(dat)[1]) - 1) * time.step, 
           y = ((1:dim(dat)[2]) - 1) * pixel.size, 
           as.matrix(dat), 
           col = tim.colors(256), 
           ylab = expression(paste("Length from tip (", mu, "m)")), 
           xlab = "Time (min)", legend.lab = "Ratiometric fluorescence (YFP/CFP)", 
           smallplot= c(.865,.885,0.2,0.98), bigplot = c(0.105,0.85,0.2,0.98), useRaster = TRUE, legend.line = 3.6)
dev.off()

# SECTION 3.4 STEP 1 & 2 - RECONSTRUCT INTRACELLULAR GRADIENT ##################
# Collapse all time points to a 2 column matrix
gd <- gather(dat, key, value)
x <- as.numeric(gd$key)
y <- gd$value

# Plot and apply loess
p <- ggplot(gd,aes(x=as.numeric(key),y=value)) + 
  geom_point(color = "darkgray", alpha = 0.5) + 
  geom_smooth(method = "loess", span = 21 /dim(dat)[2], color = "black", se = T) +
  theme_minimal() 

# Get loess data
loess.fit <- ggplot_build(p)$data[[2]]
loess.max.ind <- which.max(loess.fit$y)
loess.max.x <- loess.fit$x[loess.max.ind]
loess.max <- max(loess.fit$y)
loess.min.ind <- which.min(loess.fit$y[-c(1:loess.max.ind)]) + loess.max.ind
loess.min.x <- loess.fit$x[loess.min.ind]
loess.min <- min(loess.fit$y[-c(1:loess.max.ind)]) 

# Find sharpest decay after maximum
fluo.diff <- diff(loess.fit$y)
fit.n <- 15
fit.vec <- loess.max.ind:(loess.min.ind - (fit.n - 1))

fit.start <- fit.vec[which.min(sapply(fit.vec, function(i) as.numeric(mean(fluo.diff[i:(i + (fit.n - 1))]))))]
fit.end <- fit.start + (fit.n - 1)
x.start <- loess.fit$x[fit.start]
x.end <- loess.fit$x[fit.end]
fit.range <- (x>=x.start) & (x<=x.end)

# Perfor linear regression
fit <- lm(y[fit.range] ~ c(x[fit.range] - loess.max.x))

# Add to plot
p <- p + 
  annotate("text", x = 11, y = 1.71, label = paste("Max =",signif(fit$coef[[1]], 5),
                                                   "\nSlope =",signif(fit$coef[[2]], 5),
                                                   "\np =",signif(summary(fit)$coef[2,4], 5), 
                                                   "\nR2 = ",signif(summary(fit)$r.squared, 5)), color = "red") +
  annotate("text", x = 4, y = loess.max + 0.02, label = "Max", color = "blue") +
  annotate("text", x = 4, y = loess.min + 0.02, label = "Min", color = "blue") +
  geom_vline(xintercept = loess.max.x, color = "green") +
  geom_hline(yintercept = loess.max, color = "blue", linetype = "dotdash") +
  geom_hline(yintercept = loess.min, color = "blue", linetype = "dotdash") +
  geom_text(x = loess.max.x + 5, y = max(y), label = paste("Location =", round(loess.max.x, digits = 3)), color = "green") +
  labs(x= expression(paste("Length from tip (", mu, "m)")),    
       y= "Yellow Cameleon 3.6 fluorescence (YFP/CFP)") + 
  stat_smooth(data=subset(gd, 
                          (x>=x.start) & (x<=x.end)),
              method = "lm", col = "red") 

# Plot
pdf("~/Desktop/Fig6_Gradient.pdf", height = 5, width = 6) # USER DEFINED: file names / path in ""
print(p)
dev.off()