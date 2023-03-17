# Sensitivity Analyses
# devtools::install_github("ryanoisin/GenerativeEmotion")
library(Matrix)
library(GenerativeEmotion)
source("aux_functions.R")

# Vary parameter settings and check whether the resulting system still produces the phenomena of interest
# Two aims:
# A: Show that results are not sensitive to arbitrary choices (specific values)
# B: Show that results ARE sensitive to theory-violating choices (don't always produce phenom)

# We assess changes in:
# 1: Transition Matrix
# 2: Factor Loadings (Lambda)

#---------------------------------------------------------------#
#---------------------  1: Transition Matrix -------------------#
#---------------------------------------------------------------#
# number of situations and emotions
s <- 5
p <- 6

# define transition matrix probabilities
stay_prob <- .6
switch_prob <- (1-stay_prob)/(s-1)

# means of emotions in each situation
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
# vary staying probability from 0 to 0.99 (all transitions equally likely)

sps <- c(seq(0,.95,.05), .99)
mats <- list()
out <- matrix(NA, nrow = length(sps), ncol = 6)
colnames(out) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

set.seed(1)
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

sc <- 1.25
pdf("figures/Fig7_transition_matrix_sensitivity.pdf",7*sc,5*sc)
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
axis(1, at = seq(0, 1, 0.1)); axis(2, at = seq(0, 1, 0.2))
cols <- RColorBrewer::brewer.pal(6,"Set1")
for(j in 1:6){
  lines(x = sps, y = out[,j], lty = j, col = cols[j], lwd = 3)
}
legend(x = .725, y = .33, lty = c(1:6, 2), col = c(cols, "black"), legend = legtext,
       lwd = c(rep(2,6),1), bty = "n", cex = 0.90)
lines(x=c(0,1), y = c(.5,.5),lty = 2)
title(
  main = "Sensitivity Analysis Results", xlab = expression("Staying Probability" ~ alpha),
  ylab = "Phenomena Present (Proportion)", cex.lab = 1.25, cex.main = 1.5, font.main = 1
)
dev.off()


#---------------------  Analysis B -------------------#
# Enforce that diagonals bigger than off-diagonals, but otherwise allow anything;
# draw values from a Uniform(0, 1) distribution

nsim <- 1000
mats <- list()
out2 <- matrix(NA, nrow = nsim, ncol = 6)
colnames(out2) <- c("ARpos", "AR>CL", "PhiWithPos", "PhiBetNeg", "ResWithPos", "ResBetNeg")

set.seed(2)
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

# Make a table of results
tab1 <- t(rbind(apply(out2, 2, mean),apply(out2, 2, function(vec) sum(vec<.5)/length(vec))))
tab1 <- round(tab1,3)
colnames(tab1) <- c("Mean", "% absent")
rownames(tab1) <- c(paste0("Phenomena ", c(1,2,3.1,3.2,4.1,4.2)))

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
# Here we vary the factor loadings drawn from a Uniform(-5, 5) distribution, which is
# against our theoretical assumptions. We should thus see phenomena not being reproduced.

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

# make a table
tab3 <- t(rbind(apply(out4, 2, mean),
                apply(out4, 2, function(vec) sum(vec<.5)/length(vec))))

tab3 <- round(tab3,3)
colnames(tab3) <- c("Mean", "% absent")
tab3

tab <- cbind(tab1,tab2,tab3)
# transform to percentages
tab[,c(2,4,6)] <- tab[,c(2,4,6)]*100
tab
