#  File to simulate empirical time series and compare to empirical time series
# this file creates the Figures in the main text and Appendix A

# ------------------------------------------------------------
# --------------------- Load Packages  & Data-----------------
# ------------------------------------------------------------

# Install the generative model itself as a package

# devtools::install_github("ryanoisin/GenerativeEmotion")
library(GenerativeEmotion)

# load other necessary packages
library(dplyr)
library(RColorBrewer)
library(mlVAR) # used to fit MLVAR model to empirical data for Figure 4
library(qgraph)
library(ks)

source("aux_functions.R")
# source('model_functions.R')
# source('helper_functions.R')


rowids <- c(47, 72, 97)

# Note; this is empirical data originally published and shared open source by
# Rowland, Z., & Wenzel, M. (2020).
#Mindfulness and affect-network density: Does mindfulness facilitate disengagement from affective experiences in daily life?.
# Mindfulness, 11(5), 1253-1266.
# The original data file can be found here:
# The RDS is a processed version of this data, shared at https://github.com/jmbh/EmotionTimeSeries
dat <- readRDS('files/data_Rowland2020.RDS')

# ------------------------------------------------------------
# ----------------------- Figure 3 ---------------------------
# ------------------------------------------------------------
# empirical vs simulated time-series

# ----------------- Select Empirical Data --------------------

datc <- dat %>%
  dplyr::select(
    subj_id, happy, relaxed, satisfied, sad, anxious, angry
  ) %>%
  rename(id = subj_id)

# Pick first subject
datp <- datc[datc$id == 1, -1]
dat_emp <- as.matrix(datp)


# ----------- Simulate from Generative Model -----------------

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
N <- 500


# create model matrices
mm <- model_matrices(s = s,
                     p = p,
                     stay_prob = stay_prob,
                     switch_prob = switch_prob,
                     means = means,
                     theta = theta #, tmat = tmat
)

# simulate data
set.seed(1234)
dat_sim <- datagen(mm, 500, sout = TRUE)


# --- Create Figure 3: Empirical vs Simulated Time Series ----


# Specify plotting parameters
cols <- brewer.pal(3, 'Set1')
cols <- c('#2eb062', '#e08224')
bluecol <- adjustcolor(brewer.pal(6, 'Set1')[2], alpha = 0.75)

bluecols <- brewer.pal(5, 'Blues')
redcols <- brewer.pal(5, 'Reds')
cols_states <- c('white', bluecols[c(4, 5)], redcols[c(4, 5)])

# helper fucntion for plotting
draw_rect <- function(situation, situations, max = 120) {
  situations <- situations[seq(max)]
  rect(
    which(situations == situation), -25, which(situations == situation) + 1, -20,
    col = cols_states[situation], border = NA
  )
}

plot_timeseries <- function(dat, main, legend = TRUE, max = 120, second_plot = FALSE) {
  plot(
    dat[seq(max), 1], type = 'l', axes = FALSE, col = cols[1],
    xlab = '', ylab = '',
    xlim = c(0, max), ylim = c(-30, 130),
    main = main, font.main = 1, cex.main = 1.5
  )

  lines(dat[seq(max), 1], col = cols[1], lwd = 2)
  lines(dat[seq(max), 5], col = cols[2], lwd = 2)

  if (legend) {
    legend(
      'topleft',
      lty = c(1,1), lwd = 2, ncol = 1,
      col = c('#2eb062', '#e08224'), cex = 1.25,
      legend = c('Happy', 'Anxious'), bty = 'n'
    )
  }

  if (!second_plot) {
    mtext('Emotion Intensity', side = 2, line = 2.6, cex = 1.25)
    axis(2, at = seq(0, 100, length = 5), las = 2, cex.axis = 1.25)
  }

  mtext('Days', side = 1, line = 2, cex = 1.25)
  axis(1, at = seq(0, max, length.out = 6), labels = seq(0, max / 6, length.out = 6), cex.axis = 1.25)
}


# create figure

pdf('figures/Figure-simulated-empirical-timeseries.pdf', 9, 4.5)

par(mfrow = c(1, 2))
par(mar = c(3, 4, 2, 0) + 0.1)


plot_timeseries(dat_sim*0.8+15, 'Simulated Data', max = 120, legend = FALSE)
draw_rect(1, dat_sim[, 7])
draw_rect(2, dat_sim[, 7])
draw_rect(3, dat_sim[, 7])
draw_rect(4, dat_sim[, 7])
draw_rect(5, dat_sim[, 7])
legend(
  'topleft',
  fill = cols_states, ncol = 5, cex = 1.1,
  legend = paste0('S', seq(5)), bty = 'n'
)

par(mar = c(3, 2, 2, 2) + 0.1)
plot_timeseries(dat_emp, 'Empirical Data', max = 120, second_plot = TRUE)

dev.off()



# ------------------------------------------------------------
# ----------------------- Figure 4 ---------------------------
# ------------------------------------------------------------
# Empirical vs model-implied networks




# select variables
names <- colnames(dat)[c(5,7,8,9,10,12)]

# fit ML-VAR model
# out <- mlVAR(data = dat,
#              vars = names,
#              idvar = "subj_id",
#              dayvar = "dayno",
#              beepvar = "beep",
#              lags = 1)

# save output
# saveRDS(out, file="files/Rowland2020_mlVAR_model.RDS")
out <- readRDS("files/Rowland2020_mlVAR_model.RDS")

# extract parameters in format required for figure
out_pars <- f_getPars(out, names)

#---------------- Make Figure Empirical --------------------------------

# Make color vector indicating PE/NE
cols <- c("#2eb062", "#e08224")[c(1, 1, 1, 2, 2, 2)]

# Change layout so PE=left; NE=right
Npoints = 6
points = exp(pi * 1i * seq(0, 2, length.out = Npoints+1)[-1])
points.Cartesian = data.frame(x=Re(points), y=Im(points))
points.Cartesian <- as.matrix(points.Cartesian)
points.Cartesian <- rbind(points.Cartesian[-1, ], points.Cartesian[1, ])

PlotLabels <- function(x, srt=0) {

  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, labels=x, srt=srt)

  } # end

pdf("figures/Fig_4_Empirical_vs_Implied.pdf", 10, 10)

# Layout
lmat <- rbind(c(0, 1, 2),
              c(3, 5, 6),
              c(4, 7, 8))
layout(mat=lmat, widths = c(.2, 1, 1), heights = c(.2, 1, 1))

# Plot labels
PlotLabels("Lagged Effects")
PlotLabels("Residual Effects")
PlotLabels("Empirical")
PlotLabels("Model-Implied")

# --- empirical ---
mar <- c(5,3,5,3)+1
vsize <- 11
# Temporal
qgraph(out_pars$phi,
       edge.labels=TRUE,
       layout=points.Cartesian,
       fade=FALSE,
       color = cols,
       labels = names,
       theme = "colorblind",
       mar=mar,
       vsize=vsize)
mtext("Temporal Relations", line=2.8, cex=1.5)

# Residual
qgraph(out_pars$psi,
       edge.labels=TRUE,
       directed=FALSE,
       layout=points.Cartesian,
       fade=FALSE,
       color = cols,
       labels = names,
       theme = "colorblind",
       mar=mar,
       vsize=vsize)
mtext("Residual Partial Correlations", line=2.8, cex=1.5)

# --- model-implied ---
# Temporal
qgraph(round(t(mm$Phi),2),
       edge.labels=TRUE,
       layout=points.Cartesian,
       fade=FALSE,
       color = cols,
       labels = names,
       theme = "colorblind",
       mar=mar,
       vsize=vsize,
       directed = TRUE)
mtext("Temporal Relations", line=2.8, cex=1.5)

# Residual
qgraph(t(mm$Omega_pc),
       edge.labels=TRUE,
       directed=FALSE,
       layout=points.Cartesian,
       fade=FALSE,
       color = cols,
       labels = names,
       theme = "colorblind",
       mar=mar,
       vsize=vsize)
mtext("Residual Partial Correlations", line=2.8, cex=1.5)


dev.off()


# pdf("figures/Fig_Rowland2020_p6.pdf", 10, 5)
#
# par(mfrow=c(1,2))
#
# mar <- c(5,3,5,3)+1
# vsize <- 11
#
# # Temporal
# qgraph(out_pars$phi,
#        edge.labels=TRUE,
#        layout=points.Cartesian,
#        fade=FALSE,
#        color = cols,
#        labels = names,
#        theme = "colorblind",
#        mar=mar,
#        vsize=vsize)
# mtext("Temporal Relations", line=2.8, cex=1.5)
#
# # Residual
# qgraph(out_pars$psi,
#        edge.labels=TRUE,
#        directed=FALSE,
#        layout=points.Cartesian,
#        fade=FALSE,
#        color = cols,
#        labels = names,
#        theme = "colorblind",
#        mar=mar,
#        vsize=vsize)
# mtext("Residual Partial Correlations", line=2.8, cex=1.5)
#
# dev.off()


#---------------- Make Figure Model-Implied --------------------------------


pdf("figures/Fig_ModelImpliedNetworks_p6.pdf", 10, 5)

par(mfrow=c(1,2))


dev.off()


# ------------------------------------------------------------
# ----------- Figure 5 (Empirical vs Simulated Wedge) --------
# ------------------------------------------------------------

# Select empirical particpants for visualization

rowids <- c(47, 72, 97)

datc <- dat %>%
  dplyr::select(
    subj_id, happy, relaxed, satisfied, sad, anxious, angry
  ) %>%
  rename(id = subj_id)

datm <- datc %>%
  group_by(id) %>%
  mutate(
    time = seq(n())
  ) %>%
  ungroup() %>%
  group_by(id, time) %>%
  mutate(
    PE = mean(c(happy, relaxed, satisfied), na.rm = TRUE),
    NE = mean(c(sad, anxious, angry), na.rm = TRUE)
  ) %>%
  filter(!is.nan(PE) && !is.nan(NE)) # remove cases where whole rows had missing data

datm <- filter(datm, id %in% rowids)
ns <- table(datm$id)


# -------- Simulated data; Create composite Positve and Negative Emotion scores

Pos_o <- rowMeans(dat_sim[,1:3])
Neg_o <- rowMeans(dat_sim[,4:6])


Pos <- scales::rescale(Pos_o, from = range(Pos_o), to = c(0, 100))
Neg <- scales::rescale(Neg_o, from = range(Neg_o), to = c(0, 100))


# for plotting
cols <- brewer.pal(3, 'Set1')


pdf('Wedge_almost.pdf', 1.25*8, 1.25*5)

mat <- rbind(
  c(0, 7, 0, 8),
  c(2, 1, 5, 4),
  c(0, 3, 0, 6)
)
sc <- 3
nf <- layout(mat, widths = sc*c(1, 4, 1, 4), heights = sc*c(1,4,1), TRUE)


# ------------ Start EMPIRICAL ----------------------
par(mar = c(4.25,4.25,2,2), oma= c(0,0,0,0),  cex.axis = 1.5, cex.lab = 1.75)
plot(
  datm$PE, datm$NE,
  pch = 20, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE,
  xlab = 'Positive Emotion', ylab = 'Negative Emotion', cex = 1.25,
  font.main = 1, col = adjustcolor(cols[rep(1:3, ns)], alpha = 0.600)
)
axis(1)
axis(2, las = 2)

legend(
  'topright',
  legend = paste0('Participant: ', rowids),
  pch = 20, col = adjustcolor(cols, alpha = 0.60), bty = 'n'
)

par(mar = c(4.75, 0, 2.75, .5))
plot.new()
xlim <- c(0, 100)
ylim <- c(0, 0.040)
plot.window(xlim = -rev(ylim), ylim = xlim)
for (i in seq(1, 3)) {
  datm1 <- filter(datm, id == rowids[i])
  d1 <- density(datm1$NE)
  v <- cbind(d1$x, d1$y)
  q <- matrix(apply(v, 1, function(r) r%*%matrix(c(0,1,-1,0),2,2,byrow = T)),nrow(v), ncol(v), byrow  = T)
  lines(q, main = '', col = cols[i], lwd = 2)
}

par(mar = c(0, 4.75, 0.5, 2.75))
plot.new()

xlim <- c(0, 100)
ylim <- c(0, 0.040)
plot.window(xlim = -xlim, ylim = -rev(ylim))
for (i in seq(1, 3)) {
  datm1 <- filter(datm, id == rowids[i])
  d1 <- density(datm1$PE)
  v <- cbind(d1$x, d1$y)
  q <- matrix(apply(v, 1, function(r) r%*%matrix(c(-1,0,0,-1),2,2,byrow = T)),nrow(v), ncol(v), byrow  = T)
  lines(q, main = '', col = cols[i], lwd = 2)
}
# --------------- END EMPIRICAL --------------------------


# par(mar = c(5, 4, 4, 2) + 0.1)
par(mar = c(4.25,2,2,2), oma= c(0,0,0,0), cex.axis = 1.5, cex.lab = 1.75)
plot(
  Pos, Neg,
  pch = 20, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE,
  font.main = 1,xlab = 'Positive Emotion', ylab = 'Negative Emotion',
  col = adjustcolor('black', alpha = 0.60), cex = 1.25
)
axis(1)
axis(2, las = 2)

# par(mar = c(5, 2, 4, 1.50) + 0.1)
par(mar = c(4.75, 0, 2.75, .5))
# oma= c(0,0,0,0),  cex.axis = 1.5, cex.lab = 1.75)
plot.new()
xlim <- c(0, 100)
ylim <- c(0, 0.025)
plot.window(xlim = -rev(ylim), ylim = xlim)
dpos <- density(Pos)
dneg <- density(Neg)

dposx <- scales::rescale(dpos$x, from = range(dpos$x), to = c(0, 100))
dnegx <- scales::rescale(dneg$x, from = range(dneg$x), to = c(0, 100))

v1 <- cbind(dpos$x, dpos$y)
q1 <- matrix(apply(v1, 1, function(r) r%*%matrix(c(0,1,-1,0),2,2,byrow = T)),nrow(v), ncol(v), byrow  = T)
lines(q1, main = '', col = adjustcolor('black', alpha = 0.60), lwd = 2)

# bottom density
# par(mar = c(5, 4, 0, 2) + 0.1)
par(mar = c(0, 4.75, 0.5, 2.75))
plot.new()
plot.window(xlim = -xlim, ylim = -rev(c(0.00, 0.025)))
v2 <- cbind(dneg$x, dneg$y)
q2 <- matrix(apply(v2, 1, function(r) r%*%matrix(c(-1,0,0,-1),2,2,byrow = T)),nrow(v), ncol(v), byrow  = T)
lines(q2, main = '', col = adjustcolor('black', alpha = 0.60), lwd = 2)
# add titles

plot.new()
mtext("(a) Model-Implied", line = -5, cex = 1.25)

plot.new()
mtext("(b) Empirical", line = -5, cex = 1.25)

dev.off()


# ------------------------------------------------------------
# ----- Create Figure 6 (appendix) Minute vs Hr timescale ----
# ------------------------------------------------------------


alphas_hours <- c(0.30, 0.60, 0.90)
alphas_minutes <- alphas_hours^(1/60)
hours <- seq(0, 12)
minutes <- seq(0, 12, 0.01) * 60
cols <- brewer.pal(3, 'Set1')

pdf('figures/Figure-staying-probability.pdf', 9, 4)
par(mfrow = c(1, 2))
par(mar = c(3, 4, 2, 0) + 0.1)
plot(
  hours, alphas_hours[1]^hours, axes = FALSE, col = cols[1],
  ylab = '', xlab = '', pch = 20, cex = 1.25,
  main = 'Staying probability (hours)', font.main = 1, cex.main = 1.5
)
points(hours, alphas_hours[2]^hours, lwd = 2, col = cols[2], pch = 20)
points(hours, alphas_hours[3]^hours, lwd = 2, col = cols[3], pch = 20)

mtext('Probability', side = 2, line = 2.6, cex = 1.25)
mtext('Hours', side = 1, line = 2, cex = 1.25)

axis(1)
axis(2, las = 2)
legend(
  'topright',
  col = cols, cex = 1.2,
  lwd = 2,
  legend = c(
    expression(alpha ~ ' = 0.30'),
    expression(alpha ~ ' = 0.60'),
    expression(alpha ~ ' = 0.90')
  ), bty = 'n'
)

par(mar = c(3, 2, 2, 2) + 0.1)
plot(
  minutes, alphas_minutes[1]^minutes, lwd = 2, type = 'l', axes = FALSE, col = cols[1],
  ylab = '', xlab = '',
  main = 'Staying probability (minutes)', font.main = 1, cex.main = 1.5
)
lines(minutes, alphas_minutes[2]^minutes, lwd = 2, col = cols[2])
lines(minutes, alphas_minutes[3]^minutes, lwd = 2, col = cols[3])

mtext('Minutes', side = 1, line = 2, cex = 1.25)

axis(1, at = seq(0, 720, 120))
legend(
  'topright',
  col = cols, cex = 1.2,
  lwd = 2,
  legend = c(
    expression(alpha ~ ' = 0.9801'),
    expression(alpha ~ ' = 0.9915'),
    expression(alpha ~ ' = 0.9982')
  ), bty = 'n'
)

dev.off()
