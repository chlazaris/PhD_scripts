# 20140130. This is to check if there is
# any kind of problem with range restriction
# The idea is to check if between-strain correlation
# falls within one, two or more than two
# deviations from the mean within-strain standard deviation 

# Let's calculate the between-strain standard deviation
# for each phenotype based on the Phenotypes_M
# file.

#square_calc <- function(x){x**2}


Phenotypes_M <- read.csv("/Users/chlazaris/Desktop/Sackler_PhD/ROTATIONS/ROTATION_2/Data/Data_files/Phenotypes_M.csv", sep=",", header=T)

head(Phenotypes_M)
dim(Phenotypes_M)
# 374 which is the number of strains
# times 188 (and because the first column)
# is the strains, there are 187 phenotypes 
# in total

# So to get the phenotypes only
phenotypes <- Phenotypes_M[,2:188]

# to get the between-strain standard 
# deviations for all phenotypes:
between_strain_sd <- apply(phenotypes,2,sd)

# Now we will make a dataframe containing
# the phenotypes and the corresponding
# standard deviations
between_strain_sd_data <- data.frame(between_strain_sd)




