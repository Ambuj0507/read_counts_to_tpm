#set working directory
setwd("/home/........")

# Creating a sample dataframe
df <- data.frame(
  Gene=c("A_2kb","B_4kb","C_1kb","D_10kb"),
  Rep1=as.numeric(c(12,24,6,0)),
  Rep2=as.numeric(c(14,30,10,0)),
  Rep3=as.numeric(c(40,80,20,1))
)

#get the length of respective genes in kb
length<-c(2,4,1,10)

#read counts from data
counts<-df[,-1]

#Convert counts to RPK
rpk<-counts
for (i in 1:nrow(counts)) {
  rpk[i,] <- rpk[i,]/length[i]
}

# Calculate the total number of reads in each sample
total_reads <- colSums(rpk)

##########################################################
#############---------IMPORTANT-------####################
#############--CHANGE 10 to 10e6 below while dealing with real datasets--(here 10 is used as data is less)
##########################################################
scaling_factors <- total_reads/10  #CHANGE 10 TO 10e6

# Convert  to TPM
tpm<-rpk

#Convert RPK to TPM
for (i in 1:ncol(tpm)) {
  tpm[, i] <- rpk[, i]/scaling_factors[i]
}

#change column names of rpkm values
colnames(tpm)<-c("Rep1_tpm","Rep2_tpm","Rep3_tpm")

#merge rpm dataframe with genes
tpm_df<-cbind(df,tpm)

#save the dataframe
write.csv(tpm,"tpm.csv",row.names = F)
