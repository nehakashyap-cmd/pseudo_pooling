library(dada2); packageVersion("dada2")
library(ggplot2)
library(vegan)


exp_run<- function(file_path) { # create a function with the name my_function
  print("Hello World!")
}

fnF <- list.files("/Users/jrabasc/Desktop/PP_data/full/pilot_data", pattern="_R1", full.names = TRUE)

plotQualityProfile(fnF[c(3,13,33)])

exp_run<- function(file_path) { # create a function with the name my_function
  print("Hello World!")
}


# Filter and truncate to 240 based on quality profile
filtF <- file.path("/Users/jrabasc/Desktop/psuedo_pooling/full/pilot_data/filtered", basename(fnF))
track <- filterAndTrim(fnF, filtF, truncLen=100, maxEE=2, multithread=TRUE)
# Learn error rates
err <- learnErrors(filtF, multithread=TRUE, verbose=0)
plotErrors(err, nominalQ=TRUE) # sanity check

# Sanity check looks normal, so go ahead and read in dereplicatd sequences
drp <- derepFastq(filtF)
# Parse sample names out of filenames, and add to drp
sams <- sapply(strsplit(basename(fnF), "_"), `[`, 1)
names(drp) <- sams

system.time(dd <- dada(drp, err=err, multithread=TRUE, pool=FALSE, verbose=0)) # default
##     user   system  elapsed 
## 1062.113   24.500  185.825
system.time(dd.pseudo <- dada(drp, err=err, multithread=TRUE, pool="pseudo", verbose=0))
##     user   system  elapsed 
## 3192.272   59.492  506.946
system.time(dd.pool <- dada(drp, err=err, multithread=TRUE, pool=TRUE, verbose=0))


st <- removeBimeraDenovo(makeSequenceTable(dd), multithread=TRUE, verbose=TRUE)
## Identified 412 bimeras out of 942 input sequences.
st.pseudo <- removeBimeraDenovo(makeSequenceTable(dd.pseudo), multithread=TRUE, verbose=TRUE)
## Identified 439 bimeras out of 944 input sequences.
st.pool <- removeBimeraDenovo(makeSequenceTable(dd.pool), multithread=TRUE, verbose=TRUE)


nsam <- length(filtF)
df.obs <- data.frame(observed=c(rowSums(st>0), rowSums(st.pseudo>0), rowSums(st.pool>0)),
                     mode=rep(c("independent", "pseudo", "pooled"), each=nsam),
                     rank=rank(rowSums(st.pool>0)), times=3)
#keep <- !grepl("Mock", sams) & rowSums(st)>1000
keep <- !grepl("H2O", sams) & rowSums(st)>1000 & !grepl("NC", sams) & !grepl("RB", sams)
ggplot(data=df.obs, aes(x=rank, y=observed, color=mode)) + geom_point() +
  xlab("Samples") + ylab("Observed ASVs")


# Give unique sample names to each mode
rownames(st) <- paste0(sams, ".ind")
rownames(st.pseudo) <- paste0(sams, ".pseudo")
rownames(st.pool) <- paste0(sams, ".pool")
# Merge into one big sequence table
sta <- mergeSequenceTables(st, st.pseudo, st.pool) # 
# Calculate distances
get.dists <- function(i, sta) {
  dists <- as(vegdist(sta[c(i,i+nsam,i+nsam+nsam),]), "matrix")
  c(ind.pseudo=dists[1,2], ind.pool=dists[1,3], pseudo.pool=dists[2,3])
}
dists <- as.data.frame(t(sapply(seq(nsam), get.dists, sta=sta)))
nonmetric <- apply(dists, 1, function(xx) max(xx) > 0.5*sum(xx)); sum(nonmetric) # 1
## [1] 1
# Calculate MDS positions by hand, fixing the pooled point to lie at (0,0)
x.pool <- rep(0, nrow(dists)); y.pool <- rep(0, nrow(dists))
x.ind <- -dists$ind.pool; y.ind <- rep(0, nrow(dists))
# Trig
x.pseudo <- (-dists$ind.pool^2 + dists$ind.pseudo^2 - dists$pseudo.pool^2)/(2*dists$ind.pool)
y.pseudo <- sqrt(dists$pseudo.pool^2 - x.pseudo^2)
## Warning in sqrt(dists$pseudo.pool^2 - x.pseudo^2): NaNs produced
# Check our arithmetic
all.equal(sqrt(x.pseudo^2 + y.pseudo^2)[!nonmetric], dists$pseudo.pool[!nonmetric]) # TRUE
## [1] TRUE
all.equal(sqrt((x.pseudo-x.ind)^2 + (y.pseudo-y.ind)^2)[!nonmetric], dists$ind.pseudo[!nonmetric]) # TRUE
## [1] TRUE

# Make plotting data.frame
df <- data.frame(x = c(x.ind, x.pseudo, x.pool), y=c(y.ind, y.pseudo, y.pool),
                 mode=rep(c("independent", "pseudo", "pooled"), each=nsam), 
                 type="all", include="all", sample=c(sams, sams, sams),
                 stringsAsFactors=FALSE)
#plotdf <- df[!grepl("Mock", sams) & rowSums(st)>1000,] # Recyling the values 3x here
plotdf  <- df[!grepl("H2O", sams) & rowSums(st)>1000 & !grepl("NC", sams) & !grepl("RB", sams),]
# Randomly pick 3 example samples to plot separately
set.seed(100); NEX <- 3
examples <- sample(unique(plotdf$sample), NEX)
plotdf <- rbind(plotdf, plotdf[plotdf$sample %in% examples,])
example.rows <- seq(nrow(plotdf)-3*NEX+1,nrow(plotdf))
plotdf$type[example.rows] <- "single"
plotdf$include[example.rows] <- plotdf$sample[example.rows]
# Plot!
library(ggplot2)
ggplot(data=plotdf, aes(x=x, y=y, color=mode, alpha=type)) + geom_point() + 
  facet_grid(include~.) + scale_alpha_manual(values=c("all"=0.2, "single"=1)) +
  theme_bw() + coord_fixed(ratio=1) + guides(alpha=FALSE) + theme(panel.grid=element_blank()) +
  xlab("Distance (Bray-Curtis)") + ylab("Distance")