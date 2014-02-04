# 20131209 - R code to inspect the genotype files

#setwd("../../Phenotype_files/")

setwd("Desktop/ROTATIONS/ROTATION_2/Phenotype_files/")

# Calculate mean
pheno_mean <- read.delim("Phenotypes_M.csv", header=T, sep=",")
# 374 * 188 (374 characteristics for 188 strains)

pheno_res  <- read.delim("Phenotypes_Resid.csv", header=T, sep=",")
# 374 * 188 (374 characteristics for 188 strains)

# calculate standard deviation
pheno_sd   <- read.delim("Phenotypes_SD.csv", header=T)
# 374 st. dev

# Great, now everything is fine and it makes sense.
# We will go through all rows (traits) and we will find the correlation
# for every two traits (each row with all the others). We will report
# back the characteristics compared and the resulting Pearson
# coefficient.

char1_mean_matrix <- matrix(rep(0,374), nrow=374, ncol=187)
char2_mean_matrix <- matrix(rep(0,374), nrow=374, ncol=187)
correl_matrix <- matrix(rep(0,139876), nrow=374, ncol=374)

for (i in 1:length(row.names(pheno_mean)))#{print (i)}
	#pheno_mean[i,]

    # to get a numeric vector with all the values of mean for each characteristic
	# for all strains 
	char1_mean_matrix[i,] <- as.numeric(pheno_mean[i,2:188]) # as you have 188 strains and start with the second column
	#char1_name[i] <- pheno_mean[i,1] # This gives the content of the first column which is the name
	
	for (n in 1:length(row.names(pheno_mean)))
		char2_mean_matrix[n,] <- as.numeric(pheno_mean[n,2:188])
		#char2_name[n] <- pheno_mean[n,]
		
		r <- cor(char1_mean_matrix[i,],char2_mean_matrix[n,],use="pairwise.complete.obs",method="pearson")
		#comp_name <- paste(char1_name,char2_name,sep="_")
		
		correl_matrix[i,n] <- r
		
head(correl_matrix)

# Do it with toy dataset

matrix_1 <- matrix(rep(0,100), nrow=0, ncol=0)
matrix_1 <- matrix(rep(0,100), nrow=0, ncol=0)


   

