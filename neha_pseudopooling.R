#!/usr/bin/env Rscript

library(dada2); packageVersion("dada2")
library(ggplot2)
library(vegan)
theme_set(theme_bw())

#Install profvis package
install.packages("profvis")
library(profvis)


#instead of spaces in your path add _, spaces can lead to errors in some programs
#path <- "C:\\Users\\krrav\\Desktop\\PreMiEr REU Material\\Callahan Lab\\miseqsopdata\\MiSeq_SOP"
path <- "/home5/jrabasc/amp_seq_sim/StabilityNoMetaG_from_dada2_paper/"

list.files(path)




################ FORWARD READS
#Get forward read files and visualize quality profile
fnFs <- list.files(path, pattern="_R1", full.names = TRUE)
plotQualityProfile(fnFs[c(3,13,33)])
#Filter and truncate to 240 based on quality profile
filtFs <- file.path(path, "filtered", basename(fnFs))
trackFs <- filterAndTrim(fnFs, filtFs, truncLen=240, maxEE=2, multithread=FALSE) #Set multithread = FALSE on Windows
head(trackFs)
#Learn error rates and visualize
errF <- learnErrors(filtFs, multithread=TRUE, verbose=0)
plotErrors(errF, nominalQ=TRUE)
#If visualized error rates look normal, read in dereplicated sequences
drpF <- derepFastq(filtFs)
#Parse sample names out of filenames, and add to drp
samplenames.F <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
names(drpF) <- samplenames.F

#Profiling pseudo-pooling for forward reads
p_ddF <- profvis({
system.time(ddF <- dada(drpF, err=errF, multithread=TRUE, pool=FALSE, verbose=0))
})
ddF[[1]]
prof_output(/home5/jrabasc/neha_pseudopooling/pp_prof.output/)
htmlwidgets::saveWidget(p_ddF, "ddF_profile.html")
# Can open in browser from R
#browseURL("profile.html")

p_ddF.pseudo <- profvis({
system.time(ddF.pseudo <- dada(drpF, err=errF, multithread=TRUE, pool="pseudo", verbose=0))
})
ddF.pseudo[[1]]
prof_output(/home5/jrabasc/neha_pseudopooling/pp_prof.output/)
htmlwidgets::saveWidget(p_ddF.pseudo, "ddF_pseudo_profile.html")

p_ddF.pool <- profvis({
system.time(ddF.pool <- dada(drpF, err=errF, multithread=TRUE, pool=TRUE, verbose=0))
})
ddF.pool[[1]]
prof_output(/home5/jrabasc/neha_pseudopooling/pp_prof.output/)
htmlwidgets::saveWidget(p_ddF.pool, "ddF_pool_profile.html")

#Make sequence table and remove bimeras
stF <- removeBimeraDenovo(makeSequenceTable(ddF), multithread=TRUE, verbose=TRUE)
stF.pseudo <- removeBimeraDenovo(makeSequenceTable(ddF.pseudo), multithread=TRUE, verbose=TRUE)
stF.pool <- removeBimeraDenovo(makeSequenceTable(ddF.pool), multithread=TRUE, verbose=TRUE)

#Visualizing number of observed ASVs in each sample for each processing mode
nsamF <- length(filtFs)
df.obsF <- data.frame(observed=c(rowSums(stF>0), rowSums(stF.pseudo>0), rowSums(stF.pool>0)),
                      mode=rep(c("independent", "pseudo", "pooled"), each=nsamF),
                      rank=rank(rowSums(stF.pool>0)), times=3)
keep <- !grepl("Mock", samplenames.F) & rowSums(stF)>1000

#save this guy
ASV_table_F <- ggplot(data=df.obsF, aes(x=rank, y=observed, color=mode)) + geom_point() +
  xlab("Samples") + ylab("Observed ASVs")
ggsave(filename = "/home5/jrabasc/neha_pseudopooling/pp_plotanalysis/asv_table_F.png", plot = ASV_table_F)

### Bray-curtis distances between samples processed in each mode
# Give unique sample names to each mode
rownames(stF) <- paste0(samplenames.F, ".ind")
rownames(stF.pseudo) <- paste0(samplenames.F, ".pseudo")
rownames(stF.pool) <- paste0(samplenames.F, ".pool")
# Merge into one big sequence table
staF <- mergeSequenceTables(stF, stF.pseudo, stF.pool) #
# Calculate distances
get.distsF <- function(i, staF) {
  distsF <- as(vegdist(staF[c(i,i+nsamF,i+nsamF+nsamF),]), "matrix")
  c(ind.pseudoF=distsF[1,2], ind.poolF=distsF[1,3], pseudo.poolF=distsF[2,3])
}
distsF <- as.data.frame(t(sapply(seq(nsamF), get.distsF, sta=staF)))
nonmetricF <- apply(distsF, 1, function(xx) max(xx) > 0.5*sum(xx)); sum(nonmetricF) # 1

# Calculate MDS positions by hand, fixing the pooled point to lie at (0,0)
x.poolF <- rep(0, nrow(distsF)); y.poolF <- rep(0, nrow(distsF))
x.indF <- -distsF$ind.poolF; y.indF <- rep(0, nrow(distsF))
# Trig
x.pseudoF <- (-distsF$ind.poolF^2 + distsF$ind.pseudoF^2 - distsF$pseudo.poolF^2)/(2*distsF$ind.poolF)
y.pseudoF <- sqrt(distsF$pseudo.poolF^2 - x.pseudoF^2)

# Check our arithmetic
all.equal(sqrt(x.pseudoF^2 + y.pseudoF^2)[!nonmetricF], distsF$pseudo.poolF[!nonmetricF]) # TRUE

all.equal(sqrt((x.pseudoF-x.indF)^2 + (y.pseudoF-y.indF)^2)[!nonmetricF], distsF$ind.pseudoF[!nonmetricF]) # TRUE

# Make plotting data.frame
df.F <- data.frame(x = c(x.indF, x.pseudoF, x.poolF), y=c(y.indF, y.pseudoF, y.poolF),
                   modeF=rep(c("independent", "pseudo", "pooled"), each=nsamF),
                   type="all", include="all", sampleF=c(samplenames.F, samplenames.F, samplenames.F),
                   stringsAsFactors=FALSE)
plotdf.F <- df.F[!grepl("Mock", samplenames.F) & rowSums(stF)>1000,] # Recyling the values 3x here
# Randomly pick 3 example samples to plot separately
set.seed(100); NEX <- 3
examplesF <- sample(unique(plotdf.F$sampleF), NEX)
plotdf.F <- rbind(plotdf.F, plotdf.F[plotdf.F$sampleF %in% examplesF,])
example.rowsF <- seq(nrow(plotdf.F)-3*NEX+1,nrow(plotdf.F))
plotdf.F$type[example.rowsF] <- "single"
plotdf.F$include[example.rowsF] <- plotdf.F$sampleF[example.rowsF]
# Plot!
library(ggplot2)
bray_curt_plot_F <- ggplot(data=plotdf.F, aes(x=x, y=y, color=modeF, alpha=type)) + geom_point() +
  facet_grid(include~.) + scale_alpha_manual(values=c("all"=0.2, "single"=1)) +
  theme_bw() + coord_fixed(ratio=1) + guides(alpha=FALSE) + theme(panel.grid=element_blank()) +
  xlab("Distance (Bray-Curtis)") + ylab("Distance")
ggsave(filename = "/home5/jrabasc/neha_pseudopooling/pp_plotanalysis/bc_plot_F.png", plot = bray_curt_plot_F)




################ REVERSE READS
#Get reverse read files and visualize quality profile
fnRs <- list.files(path, pattern="_R2", full.names = TRUE)
plotQualityProfile(fnRs[c(3,13,33)])
#Filter and truncate to 160 based on quality profile
filtRs <- file.path(path, "filtered", basename(fnRs))
trackRs <- filterAndTrim(fnRs, filtRs, truncLen=160, maxEE=2, multithread=FALSE) #Set multithread = FALSE on Windows
head(trackRs)
#Learn error rates and visualize
errR <- learnErrors(filtRs, multithread=TRUE, verbose=0)
plotErrors(errR, nominalQ=TRUE)
#If visualized error rates look normal, read in dereplicated sequences
drpR <- derepFastq(filtRs)
#Parse sample names out of filenames, and add to drp
samplenames.R <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)
names(drpR) <- samplenames.R

#Profiling pseudo-pooling for reverse reads
p_ddR <- profvis({
system.time(ddR <- dada(drpR, err=errR, multithread=TRUE, pool=FALSE, verbose=0))
}) 
ddR[[1]]
prof_output(/home5/jrabasc/neha_pseudopooling/pp_prof.output/)
htmlwidgets::saveWidget(p_ddR, "ddR_profile.html")

p_ddR.pseudo <- profvis({
system.time(ddR.pseudo <- dada(drpR, err=errR, multithread=TRUE, pool="pseudo", verbose=0))
})
ddR.pseudo[[1]]
prof_output(/home5/jrabasc/neha_pseudopooling/pp_prof.output/)
htmlwidgets::saveWidget(p_ddR.pseudo, "ddR_pseudo_profile.html")

p_ddR.pool <- profvis({
system.time(ddR.pool <- dada(drpR, err=errR, multithread=TRUE, pool=TRUE, verbose=0))
})
ddR.pool[[1]]
prof_output(/home5/jrabasc/neha_pseudopooling/pp_prof.output/)
htmlwidgets::saveWidget(p_ddR.pool, "ddR_pool_profile.html")

#Make sequence table and remove bimeras
stR <- removeBimeraDenovo(makeSequenceTable(ddR), multithread=TRUE, verbose=TRUE)
stR.pseudo <- removeBimeraDenovo(makeSequenceTable(ddR.pseudo), multithread=TRUE, verbose=TRUE)
stR.pool <- removeBimeraDenovo(makeSequenceTable(ddR.pool), multithread=TRUE, verbose=TRUE)

#Visualizing number of observed ASVs in each sample for each processing mode
nsamR <- length(filtRs)
df.obsR <- data.frame(observed=c(rowSums(stR>0), rowSums(stR.pseudo>0), rowSums(stR.pool>0)),
                      mode=rep(c("independent", "pseudo", "pooled"), each=nsamR),
                      rank=rank(rowSums(stR.pool>0)), times=3)
keep <- !grepl("Mock", samplenames.R) & rowSums(stR)>1000

#save this guy
ASV_table_R <- ggplot(data=df.obsR, aes(x=rank, y=observed, color=mode)) + geom_point() +
  xlab("Samples") + ylab("Observed ASVs")
ggsave(filename = "/home5/jrabasc/neha_pseudopooling/pp_plotanalysis/asv_table_R.png", plot = ASV_table_R)

### Bray-curtis distances between samples processed in each mode
# Give unique sample names to each mode
rownames(stR) <- paste0(samplenames.R, ".ind")
rownames(stR.pseudo) <- paste0(samplenames.R, ".pseudo")
rownames(stR.pool) <- paste0(samplenames.R, ".pool")
# Merge into one big sequence table
staR <- mergeSequenceTables(stR, stR.pseudo, stR.pool) #
# Calculate distances
get.distsR <- function(i, staR) {
  distsR <- as(vegdist(staR[c(i,i+nsamR,i+nsamR+nsamR),]), "matrix")
  c(ind.pseudo=distsR[1,2], ind.pool=distsR[1,3], pseudo.pool=distsR[2,3])
}
distsR <- as.data.frame(t(sapply(seq(nsamR), get.distsR, sta=staR)))
nonmetricR <- apply(distsR, 1, function(xx) max(xx) > 0.5*sum(xx)); sum(nonmetricR) # 1

# Calculate MDS positions by hand, fixing the pooled point to lie at (0,0)
x.poolR <- rep(0, nrow(distsR)); y.poolR <- rep(0, nrow(distsR))
x.indR <- -distsR$ind.poolR; y.indR <- rep(0, nrow(distsR))
# Trig
x.pseudoR <- (-distsR$ind.poolR^2 + distsR$ind.pseudoR^2 - distsR$pseudo.poolR^2)/(2*distsR$ind.poolR)
y.pseudoR <- sqrt(distsR$pseudo.poolR^2 - x.pseudoR^2)

# Check our arithmetic
all.equal(sqrt(x.pseudoR^2 + y.pseudoR^2)[!nonmetricR], distsR$pseudo.poolR[!nonmetricR]) # TRUE

all.equal(sqrt((x.pseudoR-x.indR)^2 + (y.pseudoR-y.indR)^2)[!nonmetricR], distsR$ind.pseudoR[!nonmetricR]) # TRUE


# Make plotting data.frame
df.R <- data.frame(x = c(x.indR, x.pseudoR, x.poolR), y=c(y.indR, y.pseudoR, y.poolR),
                   modeR=rep(c("independent", "pseudo", "pooled"), each=nsamR),
                   type="all", include="all", sampleR=c(samplenames.R, samplenames.R, samplenames.R),
                   stringsAsFactors=FALSE)
plotdf.R <- df.R[!grepl("Mock", samplenames.R) & rowSums(stR)>1000,] # Recyling the values 3x here
# Randomly pick 3 example samples to plot separately
set.seed(100); NEX <- 3
examplesR <- sample(unique(plotdf.R$sampleR), NEX)
plotdf.R <- rbind(plotdf.R, plotdf.R[plotdf.R$sampleR %in% examplesR,])
example.rowsR <- seq(nrow(plotdf.R)-3*NEX+1,nrow(plotdf.R))
plotdf.R$type[example.rowsR] <- "single"
plotdf.R$include[example.rowsR] <- plotdf.R$sampleR[example.rowsR]
# Plot!
library(ggplot2)
#save this guy
bray_curt_plot_R <- ggplot(data=plotdf.R, aes(x=x, y=y, color=modeR, alpha=type)) + geom_point() +
  facet_grid(include~.) + scale_alpha_manual(values=c("all"=0.2, "single"=1)) +
  theme_bw() + coord_fixed(ratio=1) + guides(alpha=FALSE) + theme(panel.grid=element_blank()) +
  xlab("Distance (Bray-Curtis)") + ylab("Distance")
ggsave(filename = "/home5/jrabasc/neha_pseudopooling/pp_plotanalysis/bc_plot_R.png", plot = bray_curt_plot_R)
