# https://2-bitbio.com/2017/04/clustering-rnaseq-data-making-heatmaps.html

TPM <- read.table("~/TPM_Pg653_Cambium_TRD_4april2017_nr.txt", header=TRUE, sep="\t",
                  row.names=1,
                  as.is=TRUE)

class(TPM)
# [1] "data.frame"

dim(TPM)
# [1] 90145    30

TPM_MJ <- dplyr::select(TPM, starts_with("M"))

# Centers and scales data.
TPM_MJ_scaled <- t(scale(t(TPM_MJ)))
str(TPM_MJ_scaled)

TPM_MJ_scaled <- TPM_MJ_scaled[complete.cases(TPM_MJ_scaled),]


# Clusters columns by Spearman correlation.
hc <- hclust(as.dist(1-cor(TPM_MJ_scaled, method="spearman")), method="complete")

TreeC = as.dendrogram(hc, method="average")

plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")





wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))

for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     
                                     centers=i)$withinss)

plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



library(cluster)

sil <- rep(0, 20)

#repeat k-means for 1:20 and extract silhouette:

for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])
}



# Plot the  average silhouette width

plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")

abline(v = which.max(sil), lty = 2)





#clustering the data



set.seed(20)

kClust <- kmeans(scaledata, centers=4, nstart = 1000, iter.max = 20)

kClusters <- kClust$cluster



# function to find centroid in cluster i

clust.centroid = function(i, dat, clusters) {
  
  ind = (clusters == i)
  
  colMeans(dat[ind,])
  
}

kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)





library(ggplot2)

library(reshape)

#get in long form for plotting

Kmolten <- melt(kClustcentroids)

colnames(Kmolten) <- c('sample','cluster','value')



#plot

p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) +
  
  geom_point() +
  
  geom_line() +
  
  xlab("Time") +
  
  ylab("Expression") +
  
  labs(title= "Cluster Expression by Time",color = "Cluster")

p1










#------------------------------------------------------------------------



#filter transcripts



tmp <- filter.std(data,min.std=1)  # otherwise transcript with no variation will produce an error in the clustering step



data.r <- filter.NA(tmp, thres=1)

data.f <- fill.NA(data.r,mode="mean")

data.s <- standardise(data.f)



cl <- mfuzz(data.s,c=16,m=1.05)



windows(11,11)

mfuzz.plot(data.s,cl=cl,mfrow=c(4,4))





mfuzz.plot(data.s,cl=cl,mfrow=c(4,4),time.labels=seq(0,160,10))





m1 <- mestimate(data.s)

m1

#1.05









###------------------------------------------------------------------------------



#clustering option 2

library(SummarizedExperiment)



counts <- assays(data)$counts

y <- as.matrix((counts))

y <- DGEList(counts = y, group=c(1,2,3,4,5,6,7,8,9,10))

y <- calcNormFactors(y)

z <- cpm(y, normalized.lib.size=TRUE)



# https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html



# Filtering



#mean/variance calculations

z_var <- apply(z, 1, var)

z_mean <- apply(z, 1, mean)

plot(log2(z_mean), log2(z_var), pch='.')

abline(h=log2(50), col='red')

abline(v=log2(50), col='red')

text(x=13,y=23, labels="variance > 50 &\n mean > 50", col='red')
