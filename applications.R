#Codes for application-----------------------
##HoPS+ study 
#function of calculating required number of clusters
#continuous outcomes 
#a is signifcance level, b is 1-power, beta is log OR, k is cluster size, ga is rank ICC
size_cont <- function(a, b, beta, k, ga){
  Dg <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta^2 * (1 + (k - 1) * ga)
  ceiling(sqrt(1/4 + Dg^2/4) + Dg/2)
}
#ordinal outcomes 
#ln is a vector of category proportions 
size_ord <- function(a, b, beta, k, ga, ln){
  Dg <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta^2 * (1 + (k - 1) * ga)
  ceiling(Dg / (1 - sum(ln ^ 3)))
}

#functions of calculating required cluster sizes 
#continuous outcomes 
sizek_cont <- function(a, b, beta, n, ga){
  D <- beta ^ 2 / 6 / (qnorm(1-a/2) + qnorm(1-b))^2
  Dg <- (n * ga - n) / (D * n ^ 2 - n * ga) / 2
  ceiling(sqrt(D / (4 * D * n ^ 2 - 4 * n * ga) + Dg ^ 2) - Dg)
}
#oridnal outcomes
sizek_ord <- function(a, b, beta, n, ga, ln){
  Dg <- 6 * (qnorm(1-a/2) + qnorm(1-b))^2 / beta^2 / (1 - sum(ln ^ 3))
  ceiling(Dg * (ga - 1) / (Dg * ga - n))
}


####try different cluster sizes and treatment effect when treating the outcome as continuous
thetas <- seq(0.56, 0.66, 0.01)
a <- 0.05;b <- 0.15
ga <- (0.072+0.077)/2
ks <- c(25, 45, 65)
m_cont <- matrix(NA, ncol = length(ks), nrow = length(thetas))
betas <- rep(NA, length(thetas))
for(i in seq_along(thetas)){
  beta <- approxBeta(thetas[i])
  betas[i] <- beta
  for(j in seq_along(ks)){
    m_cont[i,j] <- size_cont(a, b, beta, ks[j], ga) / ks[j]
  }
}

####try different cluster sizes and treatment effect when treating the outcome as ordinal
xi <- factor(hops_f$po_prop_days_adherent_12m)
lvs <- unique(xi)
ln1 <- table(xi[hops_f$study_site_type.f == 'Control'])[lvs] / length(xi[hops_f$study_site_type.f == 'Control'])
ln2 <- table(xi[hops_f$study_site_type.f == 'Intervention'])[lvs] / length(xi[hops_f$study_site_type.f == 'Intervention'])
ln <- (ln1 + ln2) / 2
a <- 0.05;b <- 0.15
ga <- 0.07#(0.072+0.077)/2
m_ord <- matrix(NA, ncol = length(ks), nrow = length(thetas))
betas <- rep(NA, length(thetas))
for(i in seq_along(thetas)){
  beta <- approxBeta(thetas[i])
  betas[i] <- beta
  for(j in seq_along(ks)){
    m_ord[i,j] <- size_ord(a, b, beta, ks[j], ga, ln) / ks[j]
  }
}
m_cont <- m_cont * 2 
par(mar = c(4.5, 4.5, 1.1, 2.1))
plot(exp(betas), m_cont[,1], type = "b", ylim = c(min(m_cont), max(m_cont)), pch=0, lty=1, ylab = "Total number of clusters", xlab="Odds ratios",main = " ", cex.axis=1.4, cex.lab=1.5, cex.main=1.5, lwd=1.8, cex=1.5)
for(j in 2:length(ks)){
  lines(exp(betas), m_cont[,j], pch=j, lty=j, type="b",lwd=1.8, cex=1.5)
}
legend(2.1, 80, lty=c(1, 4, 5), pch = c(0, 2, 3), legend = c("k=25", "k=45", "k=65"), box.lty = 0, cex=1.5, lwd=1.8, seg.len=4)


####try different numbers of cluster and treatment effect when treating the outcome as continuous 
thetas <- seq(0.56, 0.66, 0.01)
a <- 0.05;b <- 0.15
ga <- (0.072+0.077)/2
ns <- c(6, 12, 24)
m_cont <- matrix(NA, ncol = length(ns), nrow = length(thetas))
betas <- rep(NA, length(thetas))
for(i in seq_along(thetas)){
  beta <- approxBeta(thetas[i])
  betas[i] <- beta
  for(j in seq_along(ns)){
    m_cont[i,j] <- sizek_cont(a, b, beta, ns[j], ga)
  }
}

####try different numbers of cluster and treatment effect when treating the outcome as ordinal 
xi <- factor(hops_f$po_prop_days_adherent_12m)
lvs <- unique(xi)
ln1 <- table(xi[hops_f$study_site_type.f == 'Control'])[lvs] / length(xi[hops_f$study_site_type.f == 'Control'])
ln2 <- table(xi[hops_f$study_site_type.f == 'Intervention'])[lvs] / length(xi[hops_f$study_site_type.f == 'Intervention'])
ln <- (ln1 + ln2) / 2
a <- 0.05;b <- 0.15
ga <- (0.072+0.077)/2
m_ord <- matrix(NA, ncol = length(ns), nrow = length(thetas))
betas <- rep(NA, length(thetas))
for(i in seq_along(thetas)){
  beta <- approxBeta(thetas[i])
  betas[i] <- beta
  for(j in seq_along(ns)){
    m_ord[i,j] <- sizek_ord(a, b, beta, ns[j], ga, ln)
  }
}

par(mar = c(4.5, 4.5, 1.1, 2.1))
idx <- which(m_cont[,1] > 0)
plot(exp(betas)[idx], m_cont[,1][idx], type = "b", ylim = c(0, max(m_cont, m_ord)), pch=0, lty=1, ylab = "Cluster sizes", xlab="Odds ratios", main=" ",cex.axis=1.4, cex.lab=1.5, cex.main=1.5, lwd=1.8, cex=1.5, xlim = c(exp(betas[1]), exp(betas[length(betas)])))
for(j in 2:length(ns)){
  idx <- which(m_cont[,j] > 0)
  lines(exp(betas)[idx], m_cont[,j][idx], pch=j, lty=j, type="b",lwd=1.8, cex=1.5)
}
legend(2.1, 190, lty=c(1, 4, 5), pch = c(0, 2, 3), legend = c("m=12", "m=24", "m=48"), box.lty = 0, cex=1.5, lwd=1.8, seg.len=4)



#BRIDGE study-------------------------------------------
#functions of calculating required number of clusters for one-sided tests
#continuous 
size_cont_one <- function(a, b, beta, k, ga){
  Dg <- 6 * (qnorm(1-a) + qnorm(1-b))^2 / beta^2 * (1 + (k - 1) * ga)
  ceiling(sqrt(1/4 + Dg^2/4) + Dg/2)
}
#ordinal 
size_ord_one <- function(a, b, beta, k, ga, ln){
  Dg <- 6 * (qnorm(1-a) + qnorm(1-b))^2 / beta^2 * (1 + (k - 1) * ga)
  ceiling(Dg / (1 - sum(ln ^ 3)))
} 
#functions of calculating required cluster size for one-sided tests
sizek_cont_one <- function(a, b, beta, n, ga){
  D <- beta ^ 2 / 6 / (qnorm(1-a) + qnorm(1-b))^2
  Dg <- (n * ga - n) / (D * n ^ 2 - n * ga) / 2
  ceiling(sqrt(D / (4 * D * n ^ 2 - 4 * n * ga) + Dg ^ 2) - Dg)
}

sizek_ord_one <- function(a, b, beta, n, ga, ln){
  Dg <- 6 * (qnorm(1-a) + qnorm(1-b))^2 / beta^2 / (1 - sum(ln ^ 3))
  ceiling(Dg * (ga - 1) / (Dg * ga - n))
}


####try different cluster sizes and treatment effect when treating the outcome as continuous 
betas <- log(seq(1.2, 2, 0.1))
thetas <- exp(betas) * (exp(betas) - betas - 1) / (exp(betas) - 1)^2
a <- 0.05;b <- 0.2
ga <- 0.14
ks <- c(20, 30, 40)
m_cont <- matrix(NA, ncol = length(ks), nrow = length(thetas))
for(i in seq_along(thetas)){
  beta <- betas[i]
  for(j in seq_along(ks)){
    m_cont[i,j] <- size_cont_one(a, b, beta, ks[j], ga) / ks[j]
  }
}

####try different cluster sizes and treatment effect when treating the outcome as ordinal
lbs_count <- unique(SeizureTrack_GTC$seize_count)
SeizureTrack_GTC$seize_count_ord <- factor(SeizureTrack_GTC$seize_count, levels = lbs_count, labels = lbs_count)
ni <- table(SeizureTrack_GTC$study_arm)
ln_tsc <- table(SeizureTrack_GTC[SeizureTrack_GTC$study_arm == 'TSC', "seize_count_ord"]) / ni['TSC']
ln_euc <- table(SeizureTrack_GTC[SeizureTrack_GTC$study_arm == 'EUC', "seize_count_ord"]) / ni['EUC']
ln <- (ln_tsc + ln_euc) / 2
m_ord <- matrix(NA, ncol = length(ks), nrow = length(thetas))
for(i in seq_along(thetas)){
  beta <- betas[i]
  for(j in seq_along(ks)){
    m_ord[i,j] <- size_ord_one(a, b, beta, ks[j], ga, ln) / ks[j]
  }
}

m_cont <- m_cont * 2
m_ord <- m_ord * 2
par(mfrow=c(1,3), mar = c(5.1, 4.6, 4.1, 1.6))

plot(exp(betas), m_cont[,1], type = "b", ylim = c(0, max(m_ord)), pch=0, lty=1, ylab = "Total number of clusters", xlab="Odds ratios", main="k=20",cex.axis=1.6, cex.lab=1.6, cex.main=1.6, lwd=1.6, cex=1.6)
lines(exp(betas), m_ord[,1], pch=3, lty=4, type="b",lwd=1.6, cex=1.6)
legend(1.5, 400, lty=c(1, 4), pch = c(0, 3), legend = c("Continuous", "Ordinal"), box.lty = 0, cex=1.6, lwd=1.6, seg.len=4)

plot(exp(betas), m_cont[,2], type = "b", ylim = c(0, max(m_ord)), pch=0, lty=1, ylab = "Total number of clusters", xlab="Odds ratios", main="k=30",cex.axis=1.6, cex.lab=1.6, cex.main=1.6, lwd=1.6, cex=1.6)
lines(exp(betas), m_ord[,2], pch=3, lty=4, type="b",lwd=1.6, cex=1.6)

plot(exp(betas), m_cont[,3], type = "b", ylim = c(0, max(m_ord)), pch=0, lty=1, ylab = "Total number of clusters", xlab="Odds ratios", main="k=40",cex.axis=1.6, cex.lab=1.6, cex.main=1.6, lwd=1.6, cex=1.6)
lines(exp(betas), m_ord[,3], pch=3, lty=4, type="b",lwd=1.6, cex=1.6)


####try different numbers of clusters and treatment effect when treating the outcome as continuous 
betas <- log(seq(1.2, 2, 0.1))
thetas <- exp(betas) * (exp(betas) - betas - 1) / (exp(betas) - 1)^2
a <- 0.05;b <- 0.2
ga <- 0.14
ns <- c(15, 30, 45)
m_cont <- matrix(NA, ncol = length(ns), nrow = length(thetas))
for(i in seq_along(thetas)){
  beta <- betas[i]
  for(j in seq_along(ks)){
    m_cont[i,j] <- sizek_cont_one(a, b, beta, ns[j], ga)
  }
}

####try different numbers of clusters and treatment effect when treating the outcome as ordinal
lbs_count <- unique(SeizureTrack_GTC$seize_count)
SeizureTrack_GTC$seize_count_ord <- factor(SeizureTrack_GTC$seize_count, levels = lbs_count, labels = lbs_count)
ni <- table(SeizureTrack_GTC$study_arm)
ln_tsc <- table(SeizureTrack_GTC[SeizureTrack_GTC$study_arm == 'TSC', "seize_count_ord"]) / ni['TSC']
ln_euc <- table(SeizureTrack_GTC[SeizureTrack_GTC$study_arm == 'EUC', "seize_count_ord"]) / ni['EUC']
ln <- (ln_tsc + ln_euc) / 2
m_ord <- matrix(NA, ncol = length(ks), nrow = length(thetas))
for(i in seq_along(thetas)){
  beta <- betas[i]
  for(j in seq_along(ns)){
    m_ord[i,j] <- sizek_ord_one(a, b, beta, ns[j], ga, ln)
  }
}

par(mfrow=c(1,3), mar = c(5.1, 4.6, 4.1, 1.6))

idx <- which(m_cont[,1] > 0)
plot(exp(betas)[idx], m_cont[,1][idx], type = "b", ylim = c(0, max(m_ord)), pch=0, lty=1, ylab = "Cluster sizes", xlab="Odds ratios", main="m=30",cex.axis=1.6, cex.lab=1.6, cex.main=1.6, lwd=1.6, cex=1.6, xlim = c(exp(betas[1]), exp(betas)[length(betas)]))
idx <- which(m_ord[,1] > 0)
lines(exp(betas)[idx], m_ord[,1][idx], pch=3, lty=4, type="b",lwd=1.6, cex=1.6)
legend(1.3, 100, lty=c(1, 4), pch = c(0, 3), legend = c("Continuous", "Ordinal"), box.lty = 0, cex=1.6, lwd=1.6, seg.len=4)

idx <- which(m_cont[,2] > 0)
plot(exp(betas)[idx], m_cont[,2][idx], type = "b", ylim = c(0, max(m_ord)), pch=0, lty=1, ylab = "Cluster sizes", xlab="Odds ratios", main="m=60",cex.axis=1.6, cex.lab=1.6, cex.main=1.6, lwd=1.6, cex=1.6, xlim = c(exp(betas[1]), exp(betas)[length(betas)]))
idx <- which(m_ord[,2] > 0)
lines(exp(betas)[idx], m_ord[,2][idx], pch=3, lty=4, type="b",lwd=1.6, cex=1.6)

idx <- which(m_cont[,3] > 0)
plot(exp(betas)[idx], m_cont[,3][idx], type = "b", ylim = c(0, max(m_ord)), pch=0, lty=1, ylab = "Cluster sizes", xlab="Odds ratios", main="m=90",cex.axis=1.6, cex.lab=1.6, cex.main=1.6, lwd=1.6, cex=1.6, xlim = c(exp(betas[1]), exp(betas)[length(betas)]))
idx <- which(m_ord[,3] > 0)
lines(exp(betas)[idx], m_ord[,3][idx], pch=3, lty=4, type="b",lwd=1.6, cex=1.6)
