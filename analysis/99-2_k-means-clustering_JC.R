# https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html

TPM <- read.table("~/TPM_Pg653_Cambium_TRD_4april2017_nr.txt", header=TRUE, sep="\t",
                  row.names=1,
                  as.is=TRUE)

class(TPM)
# [1] "data.frame"

dim(TPM)
# [1] 90145    30

head(TPM)


TPM_MJ <- dplyr::select(TPM, starts_with("M"))

dim(TPM_MJ)
# [1] 90145    15

head(TPM_MJ)

TPM_MJ <- dplyr::transmute(TPM_MJ,
                           "2H" = (M1_2H + M2_2H + M3_2H)/3,
                           "1D" = (M1_1D + M2_1D + M3_1D)/3,
                           "2D" = (M1_2D + M2_2D + M3_2D)/3,
                           "4D" = (M1_4D + M2_4D + M3_4D)/3,
                           "8D" = (M1_4D + M2_4D + M3_4D)/3) # MY: error in calculation

TPM_MJ.m <- as.matrix(TPM_MJ)

dim(TPM_MJ.m)
# [1] 90145     5


# mean/variance calculations

# MY: generate 2 named numeric vectjor
TPM_MJ_var <- apply(TPM_MJ, 1, var)
length(TPM_MJ_var)
# [1] 90145

TPM_MJ_mean <- apply(TPM_MJ, 1, mean)
length(TPM_MJ_mean)
# [1] 90145


plot(log2(TPM_MJ_mean), log2(TPM_MJ_var), pch='.')

abline(h=log2(50), col='red')

abline(v=log2(50), col='red')

text(x=13, y=23, labels="variance > 50 &\n mean > 50", col='red')


TPM_MJ_top50 <- TPM_MJ[which(TPM_MJ_var > 50 & TPM_MJ_mean > 50), 1:5]

str(TPM_MJ_top50)
# 'data.frame':	2410 obs. of  5 variables:


# first get the time point data together:
time <- c(2,24,48,96,192)

# bind that to the dataframe
TPM_MJ_scaled <- rbind(time, TPM_MJ_scaled)

write.table(TPM_MJ_scaled,file="~/TPM_MJ_scaled.txt", sep='\t', quote = F, col.names=NA)



# read it back in as an expression set
data <- table2eset(file="~/TPM_MJ_scaled.txt")

m1 <- mestimate(data)

m1

Dmin(data, m = m1, crange = seq(2,22,1), repeats=3, visu=TRUE)


# clust=10
clust=6

c <- mfuzz(data,c=clust,m=m1)


# make a plot with the clusters identified
mfuzz.plot(data,cl=c,mfrow=c(1,1),time.labels=c(2,24,48,96,192), new.window=FALSE)


mfuzz.plot2(data,cl=c,mfrow=c(4,3),time.labels=c(2,24,48,96,192),
            new.window=FALSE, min.mem=0.5, colo = "fancy")

