# Sensitivity Analyses
library(Matrix)
library(GenerativeEmotion)
source("aux_functions.R")
# source("model_functions.R")
# source("helper_functions.R")

# I want to vary my parameter settings and check whether the resulting system
# still produces the phenomena of interest
# Two aims:
# A: Show that results are not sensitive to arbitrary choices (specific values)
# B: Show that results ARE sensitive to theory-violating choices (don't always produce phenom)

# There are three things we can vary

# 1: Transition Matrix
# 2: Factor Loadings (Lambda)
# 3: Residual Variances
# 4: Number of states

#---------------------------------------------------------------#
#---------------------  1: Transition Matrix -------------------#
#---------------------------------------------------------------#

# number of situations and emotions
s <- 5
p <- 6

# define transition matrix probabilities
stay_prob <- .6
switch_prob <- (1-stay_prob)/(s-1)

# means of eomtions in each situation
a <- 90
b <- 60
c <- 45


means <- rbind(rep(0,p),
               c(a, b, c,  0  ,  0  ,  0  ),
               c(c, a, b,  0  ,  0  ,  0  ),
               c( 0  ,  0  ,  0  , a, b, c),
               c( 0  ,  0  ,  0  , c, a, b))

theta <- 16^2

#---------------------  Analysis A -------------------#
# Keeping the structure of the transition matrix fixed

# vary staying probability from .2 (all transitions equally likely)


sps <- c(seq(0,.95,.05),.99)
mats <- list()
out <- matrix(NA, nrow = length(sps), ncol = 6)
colnames(out) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

for(i in 1:length(sps)){
  stay_prob <- sps[i]
  switch_prob <- (1-stay_prob)/(s-1)
  
  mm <- model_matrices(s = s,
                       p = p,
                       stay_prob = stay_prob,
                       switch_prob = switch_prob,
                       means = means,
                       theta = theta)
  mats[[i]] <- mm
  out[i,] <- check_phenom(phi = round(mm$Phi,16), om = round(mm$Omega_pc,16))
}


legtext <- c(paste0("Phenomena ", 1:6), "Decision Threshold")

plotout <- out
plotout[,1] <- out[,1]+.008
plotout[,3] <- out[,3]-.008
plotout[6:length(sps),4] <- out[6:length(sps),4]-.014

plotout[,5] <- out[,5]-.009

sc <- 1.25
pdf("figures/tmat_sensitivity.pdf",7*sc,5*sc)
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
axis(1, at = seq(0, 1, 0.1)); axis(2, at = seq(0, 1, 0.2))
cols <- RColorBrewer::brewer.pal(6,"Set1")
for(j in 1:6){
  lines(x = sps, y = plotout[,j], lty = j, col = cols[j], lwd = 3)
}
# legend("bottomright", lty = c(rep(1,6),2), col = c(cols, "black"), legend = legtext, 
#        lwd = c(rep(2,6),1), bty = "n")
legend(x = .725, y = .33, lty = c(1:6, 2), col = c(cols, "black"), legend = legtext, 
       lwd = c(rep(2,6),1), bty = "n", cex = 0.90)
lines(x=c(0,1), y = c(.5,.5),lty = 2)
title(
  main = "Sensitivity Analysis Results", xlab = expression("Staying Probability" ~ alpha),
  ylab = "Phenomena Present (Proportion)", cex.lab = 1.25, cex.main = 1.5, font.main = 1
)
dev.off()



# # check - what's happening when sp = 0.2 with the ar/cl thing?
# qgraph(t(mats[[82]]$Phi), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .3, mar=(rep(6,4)))
# t(mats[[82]]$Phi)

# plot(unlist(lapply(mats, function(l) mean(abs(l$Phi)))))
# 
# lapply(mats, function(l) l$Pi) # equal everywhere....
# 
# plot(unlist(lapply(mats, function(l) mean(abs(l$Gamma)))))
# lapply(mats, function(l) l$Gamma)


#---------------------  Analysis B -------------------#
# Enforce that diagonals bigger than off-diagonals, but otherwise allow anything


# vary staying probability from .2 (all transitions equally likely)

nsim <- 1000

mats <- list()
out2 <- matrix(NA, nrow = nsim, ncol = 6)
colnames(out2) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")


set.seed(12354)
for(i in 1:nsim){
  tmat <- gen_tmat(p = 5, type = "diagheavy")
  mm <- model_matrices(s = s,
                       p = p,
                       stay_prob = NULL,
                       switch_prob = NULL,
                       means = means,
                       theta = theta,
                       tmat = tmat)
  mats[[i]] <- mm
  out2[i,] <- check_phenom(phi = mm$Phi, om = mm$Omega_pc)
}
# what phenomena always present?
apply(out2,2,function(col) all(col==1))
notall <- which(apply(out2,2,function(col) all(col==1))!= T)

for(j in notall){
  hist(out2[,j], main = colnames(out2)[j], xlim = c(0,1),
       xlab = "Degree of phenomena endorsement")
  abline(v = .5, col = "red", lty = 2)
  abline(v = mean(out2[,j]), col = "blue")
}

which(out2[,])
apply(out2[,notall],2,function(col) (sum(col < .5)/nsim)*100)
apply(out2[,notall],2,function(col) which(col < .5))


# Make a table of results

tab1 <- t(rbind(apply(out2, 2, mean),apply(out2, 2, function(vec) sum(vec<.5)/length(vec))))
tab1 <- round(tab1,3)
colnames(tab1) <- c("Mean", "% absent")
rownames(tab1) <- c(paste0("Phenomena ", c(1,2,3.1,3.2,4.1,4.2)))


# mostly fine - what's going on in the weird cases?
# hard to see
# sel <- apply(out[,notall], 2, function(col) which(col < .5))
# sel
# mats[[19]]$tmat
# mats[[792]]$Phi
# # 
#  j <- 282
# qgraph(t(mats[[j]]$Phi), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .3, mar=(rep(6,4)))
# qgraph(t(mats[[j]]$tmat), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .8, mar=(rep(6,4)))
# 
# mats[[892]]$tmat
# mats[[892]]$Phi
# 
# 
# qgraph(t(mats[[892]]$Phi), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .3, mar=(rep(6,4)))
# qgraph(t(mats[[892]]$tmat), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .8, mar=(rep(6,4)))
# 
# pdf("tmat_sensitivity.pdf",8,8)
# par(mfrow = c(2,2))
# plot.new()
# plot.window(xlim = c(0,1), ylim = c(0,1))
# axis(1); axis(2)
# cols <- RColorBrewer::brewer.pal(6,"Set2")
# for(j in 1:6){
#   if(j == 1){ lines(x = sps, y = (out[,j] +.01), type = "b", col = cols[j], lwd = 2, )}
#   else if(j == 3){ lines(x = sps, y = (out[,j] -.01), type = "b", col = cols[j], lwd = 2, )
#   }else{ lines(x = sps, y = out[,j], type = "b", col = cols[j], lwd = 2, )}
# }
# legend("bottomright", lty = 1, col = cols, legend = colnames(out), lwd = 2)
# abline(h = .5, lty = 2)
# title(main = "Effect of changing staying probability", xlab = "staying probability",
#       ylab = "phenomena present proportion")
# 
# for(j in notall){
#   hist(out2[,j], main = colnames(out2)[j], xlim = c(0,1),
#        xlab = "Degree of phenomena endorsement")
#   abline(v = .5, col = "red", lty = 2)
#   abline(v = mean(out2[,j]), col = "blue")
# }
# dev.off()


# to do : plot ratio of ar and cl against phenomena
# check if closesness to .2 staying probability correlates with ar > cl being absent
# check eigenvales of phi matrices
#---------------------------------------------------------------#
#---------------------  2: Lambda's ----------------------------#
#---------------------------------------------------------------#

# We want to see if the choice of specific values for the factor loadings 
# makes any difference to whether or not we produce the phenomena

# To do this, we will keep the basic structure of the lambda matrix the same
# (same zero and non-zero parameters)
# but freely choose the positive free parameter values within some bounds


nsim <- 1000

mats <- list()
out3 <- matrix(NA, nrow = nsim, ncol = 6)
colnames(out3) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

stay_prob <- .6
switch_prob <- (1-stay_prob)/(s-1)

set.seed(123)
for(i in 1:nsim){
  la <- runif(12, min = 0.1, max = 5)
  
  # factor loading matrix 
  means <- rbind(rep(0,p),
                 c(la[1], la[2],la[3],  0  ,  0  ,  0  ),
                 c(la[4], la[5], la[6],  0  ,  0  ,  0  ),
                 c( 0  ,  0  ,  0  , la[7], la[8], la[9]),
                 c( 0  ,  0  ,  0  , la[10], la[11], la[12]))
  
  mm <- model_matrices(s = s,
                       p = p,
                       stay_prob = stay_prob,
                       switch_prob = switch_prob,
                       means = means,
                       theta = theta)
  mats[[i]] <- mm
  out3[i,] <- check_phenom(phi = mm$Phi, om = mm$Omega_pc)
}
# what phenomena always present?
apply(out3,2,function(col) all(col==1))
notall <- which(apply(out3,2,function(col) all(col==1))!= T)

apply(out3[,notall], 2, mean)


# make a table

tab2 <- t(rbind(apply(out3, 2, mean),
                apply(out3, 2, function(vec) sum(vec<.5)/length(vec))))

tab2 <- round(tab2,3)
colnames(tab2) <- c("Mean", "% absent")

tab <- cbind(tab1,tab2)
tab[,c(2,4)] <- tab[,c(2,4)]*100
tab


# ---------------------------------------------------------------------------------
# ------------------------- Lambda Contra-Indication ------------------------------
# ---------------------------------------------------------------------------------

nsim <- 1000

mats <- list()
out4 <- matrix(NA, nrow = nsim, ncol = 6)
colnames(out3) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

stay_prob <- .6
switch_prob <- (1-stay_prob)/(s-1)

set.seed(123)
for(i in 1:nsim){
  la <- runif(12, min = -5, max = 5)
  
  # factor loading matrix 
  means <- rbind(rep(0,p),
                 c(la[1], la[2],la[3],  0  ,  0  ,  0  ),
                 c(la[4], la[5], la[6],  0  ,  0  ,  0  ),
                 c( 0  ,  0  ,  0  , la[7], la[8], la[9]),
                 c( 0  ,  0  ,  0  , la[10], la[11], la[12]))
  
  mm <- model_matrices(s = s,
                       p = p,
                       stay_prob = stay_prob,
                       switch_prob = switch_prob,
                       means = means,
                       theta = theta)
  mats[[i]] <- mm
  out4[i,] <- check_phenom(phi = mm$Phi, om = mm$Omega_pc)
}
# what phenomena always present?
apply(out4,2,function(col) all(col==1))
notall <- which(apply(out4,2,function(col) all(col==1))!= T)

apply(out4[,notall], 2, mean)


tab3 <- t(rbind(apply(out4, 2, mean),
                apply(out4, 2, function(vec) sum(vec<.5)/length(vec))))

tab3 <- round(tab3,3)
colnames(tab3) <- c("Mean", "% absent")
tab3

tab <- cbind(tab1,tab2,tab3)

# transform to percentages

tab[,c(2,4,6)] <- tab[,c(2,4,6)]*100
tab

# 
# pdf("lambda_sensitivity.pdf",8,8)
# for(j in notall){
#   hist(out[,j], main = colnames(out)[j], xlim = c(0,1))
#   abline(v = .5, col = "red", lty = 2)
#   abline(v = mean(out[,j]), col = "blue")
# }
# dev.off()
# 
#  qgraph(t(mats[[792]]$Phi), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .3, mar=(rep(6,4)))
#  mats[[792]]$Lambda
#  
# sel <- which(out[,"PhiBetNeg"] == min(out[,"PhiBetNeg"]))
# for(i in sel){
#  qgraph(t(mats[[i]]$Phi), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#         maximum = .3, mar=(rep(6,4)))
# print(mats[[i]]$Lambda)
# }
# 
# qgraph(t(mats[[  which(out[,"ResBetNeg"] == min(out[,"ResBetNeg"]))]]$Omega_pc), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
#        maximum = .3, mar=(rep(6,4)))
# 
# 
# # this is working fine - sometimes we get between-valence effects which are positive
# # but these are very close to zero
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   geom_density_ridges(rel_min_height = 0.005) +
#   scale_y_discrete(expand = c(0.01, 0)) +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   theme_ridges()
#  
# outm <- data.frame(Phenomena = as.vector(sapply(colnames(out), rep, nrow(out))),
#            Frequency.Endorsed = c(out[,1], out[,2], out[,3],out[,4], out[,5], out[,6]) )
# 
# pdf("lambda_sensitivity.pdf",10,8)
# ggplot(outm, aes(x = Frequency.Endorsed, y = Phenomena)) +
#    geom_density_ridges(rel_min_height = 0.005) +
#    scale_y_discrete(expand = c(0.01, 0)) +
#    scale_x_continuous(expand = c(0.01, 0), limits = c(0.3,1)) +
#   theme_ridges() + geom_vline(xintercept = 0.5, color = "red")
# dev.off()
# what if we choose the same structure, but allow negative values?
# then the meaning of "same valence" emotions disappears.... but ok
#---------------------------------------------------------------#
#----------------  3: Residual Variances -----------------------#
#---------------------------------------------------------------#

# NOTE 17 MARCH; below must be changed to new model specification

#---------------------  Analysis A -------------------#
# Keep the residual variance fixed, but scale it up

# define means in each state
a <- 1
b <- .5

# factor loading matrix 
means <- rbind(rep(0,p),
               c(a, b, a,  0  ,  0  ,  0  ),
               c(b, a, b,  0  ,  0  ,  0  ),
               c( 0  ,  0  ,  0  , a, b, a),
               c( 0  ,  0  ,  0  , b, a, b))

# keep tmat fixed to the usual
stay_prob <- .6
switch_prob <- (1-stay_prob)/(s-1)

# sequence of variances to try out
ths <- seq(.1,2,.1)
mats <- list()
out <- matrix(NA, nrow = length(ths), ncol = 6)
colnames(out) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

for(i in 1:length(ths)){
theta <- ths[i]
  
  mm <- model_matrices(s = s,
                       p = p,
                       stay_prob = stay_prob,
                       switch_prob = switch_prob,
                       means = means,
                       theta = theta)
  mats[[i]] <- mm
  out[i,] <- check_phenom(phi = mm$Phi, om = mm$Omega_pc)
}

plot.new()
plot.window(xlim = c(0,2.1), ylim = c(0,1))
axis(1); axis(2)
cols <- RColorBrewer::brewer.pal(6,"Set2")
for(j in 1:6){
  lines(x = ths, y = jitter(out[,j], amount = .01), type = "b", col = cols[j], lwd = 2, )
}
legend("bottomright", lty = 1, col = cols, legend = colnames(out), lwd = 2)
abline(h = .5, lty = 2)
title(main = "Effect of changing uniform error variance", xlab = "error sd",
      ylab = "phenomena present proportion")


#---------------------  Analysis B -------------------#
# Vary residual sd across variables

nsim <- 1000

mats <- list()
out <- matrix(NA, nrow = nsim, ncol = 6)
colnames(out) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

set.seed(127)
for(i in 1:nsim){
  theta <- diag(runif(6, .1, 1))
  
  mm <- model_matrices(s = s,
                       p = p,
                       stay_prob = stay_prob,
                       switch_prob = switch_prob,
                       means = means,
                       theta = theta)
  mats[[i]] <- mm
  out[i,] <- check_phenom(phi = mm$Phi, om = mm$Omega_pc)
}
# what phenomena always present?
apply(out,2,function(col) all(col==1))
notall <- which(apply(out,2,function(col) all(col==1))!= T)

for(j in notall){
  hist(out[,j], main = colnames(out)[j], xlim = c(0,1))
  abline(v = .5, col = "red", lty = 2)
  abline(v = mean(out[,j]), col = "blue")
}

# looks good!
sel <- which(out[,2]==min(out[,2]))
for(i in sel){
  qgraph(t(mats[[i]]$Phi), directed = T, edge.labels = T, diag = T, fade = F, layout = "circle",
         maximum = .3, mar=(rep(6,4)))
}

#---------------------  Analysis C -------------------#
# Allow residual covariances within-valence.... not sure how to do this in a good way
# library(corpcor)
# create_Sigma <- function(p, noise_sd = 1, corr = 0) {
#   
#   Sigma <- diag(noise_sd, p)
#   Sigma[upper.tri(Sigma)] <- corr*sample(c(-1, 1), p*(p - 1) / 2, replace = TRUE)
#   Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
#   
#   while (!is.positive.definite(Sigma)) {
#     Sigma <- diag(noise_sd, p)
#     Sigma[upper.tri(Sigma)] <- corr*sample(c(-1, 1), p *(p - 1) / 2, replace = TRUE)
#     Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
#   }
#   
#   Sigma
# }
# create_Sigma(5, cor = .3)
