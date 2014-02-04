# Analysis of between and within correlation
# for the remainder of the Raw.data.Rfile

# we set the working directory
setwd("/Users/chlazaris/Desktop/Sackler_PhD/ROTATIONS/ROTATION_2/Data/Data_files")

# we load the data file
data <- load("Raw.data.Rfile")

# check what the data are
str(data)
>chr "Assembled"

# check what "Assembled" is
class(Assembled)
[1] "list"

# check the length of the list
length("Assembled")
[1] 3

# Thus the list has three elements
# which based on what Kerry said, correspond 
# to 1. unbudded cells (A), 2. small buds (A1B),
# 3. large buds (C).

# Now in order to find the number of columns


names(Assembled[[2]])
names(Assembled[[3]])

length(names(Assembled[[2]]))
[1] 63

length(names(Assembled[[3]]))
[1] 104

# We want though to remove the first columns
# that do not correspond to phenotypes
# (namely plate, replicate, condition)
ass2_vars <- names(Assembled[[2]]) %in% c("plate", "replicate", "condition")
ass2 <- Assembled[[2]][!ass2_vars]

ass3_vars <- names(Assembled[[3]]) %in% c("plate", "replicate", "condition")
ass3 <- Assembled[[3]][!ass3_vars]

# So to find the final number of phenotypes, we have:
length(names(ass2))
[1] 60 # which means 60 phenotypes are going to be checked in this case

length(names(ass3))
[1] 101 # 101 phenotypes are going to be checked in this case

# Calculation of unique strains in each of the cases
uniq_strain_names <- unique(ass2$strain)
uniq_strain_names <- unique(ass3$strain)

# So as expected the number of strains is the same
# in both cases and equals 392.

# Now we are ready to proceed with the calculation of
# between and within correlation.

# the within correlation is based on the
# raw data file while the between is based
# on Phenotypes M file. We will make sure
# that we take the same strains into account
# in our calculations. So:

ass_1 <- Assembled[[1]]
unique(ass_1$strain)
# 392 strains

# while
phenotypes_M <- read.delim("Phenotypes_M.csv", sep=",", header=T)
unique(phenotypes_M$X)

# in order to find out how many strains are common between
# the two:
common_strains <- intersect(unique(phenotypes_M$X),unique(ass_1$strain))
length(common_strains)
# 374 are the common ones F2_1 to F2_374

# So now from ass_1 we create a new dataset for the calculation
# of median within correlation, which contains data only from
# the strains F2_1 to F2_374 
new_ass_1 <- ass_1[ass_1$strain %in% common_strains,]
ass1_vars <- names(ass_1) %in% c("plate", "replicate", "condition")
ass1 <- ass_1[!ass1_vars]

ass <- ass1
uniq_strain_names <- unique(common_strains)

# Prepare the necessary matrices first:
strain_name_matrix <- matrix(rep(0,length(uniq_strain_names)),nrow=length(uniq_strain_names),ncol=1)
phen_corr_matrix <- matrix(rep(0,(length(names(ass))-1)*length(uniq_strain_names)), nrow=length(uniq_strain_names), ncol=length(names(ass))-1)
no_cells_matrix <- matrix(rep(0,length(uniq_strain_names)),nrow=length(uniq_strain_names),ncol=1)

# Calculation of median within and between correlation
# for the 2nd element of the list (small buds)

# median within correlation
for (i in 1:length(uniq_strain_names)) { 
	strain_data <- ass[which(ass$strain==uniq_strain_names[i]),]
	# Get the strain name 
	strain_name <- uniq_strain_names[i]
	# print (strain_name)
	# Get the correlation between phenotype 1 and all the others
	num_cols <- length(names(ass))
	for (z in 2:num_cols) {
		phen_corr <- cor(strain_data[,2],strain_data[,z])
		#print (phen_corr)
		phen_corr_matrix[i,(z-1)] <- phen_corr
	}  
	
	# Get the number of cells for strain
	no_cells <- length(strain_data[,1])
	
	# Fill in the rest matrices
	strain_name_matrix[i] <- strain_name
	no_cells_matrix[i] <- no_cells
}

print (phen_corr_matrix)

# Combine all matrices into one final matrix
final_matrix <- cbind(strain_name_matrix,phen_corr_matrix,no_cells_matrix)
# turn it to dataframe
final_matrix <- data.frame(final_matrix)

# Check the dimensions of this matrix
dim(final_matrix)

# Calculate the medians for the correlations 
median_corr <- apply(phen_corr_matrix,2,median)

# So finally we have the median within correlation
median_within_correlation <- median_corr

# Get the absolute values
abs_median_within_cor <- abs(median_within_correlation)

# Calculation of between correlation
# ==================================

# Correction: Based on what Kerry said
# I have to use the Phenotypes_M file
# for this one.

# The new matrix must have two columns
# one for the phenotypes and one
# for the correlation coeffiecient

# From the Phenotypes_M file I subset the data that 
# correspond to A (28 first phenotypes)

# This is the important line
phen <- phenotypes_M[,89:188]


# Everything else can be written in a script
between_corr_matrix <- matrix(rep(0,2*length(names(phen))), ncol=2)

# Create a loop to report the pair
# of phenotypes (column 1) and calculate the
# between-strain correlations of every 
# pair of P1 with the other phenotypes (column 2)

for (i in 1:length(names(phen))) {
	between_corr_matrix[i,1] <- paste(names(phen)[1], names(phen)[i],sep=":")
	between_corr_matrix[i,2] <- cor(phen[1],phen[i])
}

# Now we turn the matrix to dataframe and we have the following:

between_corr_data <- data.frame(between_corr_matrix)

abs_between_corr_data <- between_corr_data

# get the absolute values for the between_corr_data as well

# to get the numeric values of a factor
abs_between_corr_data$X2 <- abs(as.numeric(levels(between_corr_data$X2)[between_corr_data$X2]))

# So for the between_corr_data
abs_between_corr_data

# Now we can plot median within and between correlation
final_corr_data <- cbind(abs_between_corr_data,abs_median_within_cor)

# remove missing values
final_corr_data <- na.omit(final_corr_data)
final_corr_data_2 <- final_corr_data[,2:3]

# We create the final dataframe (we put proper descriptive names
# to the columns)
colnames(final_corr_data) <- c("Phenotypes","abs_between_cor","abs_median_within_cor")

# now we use ggplot2 to plot median within and between correlation
library(ggplot2)

# Calculate slope and intercept of line of best fit
# line_values <- coef(lm(abs_median_within_cor ~ abs_between_cor, data = final_corr_data))

# Function for line equation
# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: http://goo.gl/K4yh

# lm_eqn <- function(df){
#     m = lm(y ~ x, df);
#     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#          list(a = format(coef(m)[1], digits = 2), 
#               b = format(coef(m)[2], digits = 2), 
#              r2 = format(summary(m)$r.squared, digits = 3)))
#     as.character(as.expression(eq));                 
# }

line_values <- coef(lm(abs_median_within_cor ~ abs_between_cor, data = final_corr_data))

corr_plot <- ggplot(final_corr_data, aes(x=abs_between_cor,y=abs_median_within_cor)) +
theme(text = element_text(size=10),axis.text.x = element_text()) + xlab("Between correlation") +
ylab("Median within correlation") + geom_point(fill="purple",alpha=0.7) + geom_text(aes(label=Phenotypes), size=2, hjust=0) + geom_abline(intercept=0, slope=1, color="blue") 

cp2 <- corr_plot + geom_abline(intercept=line_values[1], slope=line_values[2], color="red")

# save the plot
ggsave("C_Final_between_within_plot.pdf",dpi=300)





# Calculation of median within and between correlation
# for the 3rd element of the list (large buds)
