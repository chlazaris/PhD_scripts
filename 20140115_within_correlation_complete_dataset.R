# 20140115_within_strain_correlation_for_complete_dataset

# Change to the directory where the complete dataset is:
load("/Users/chlazaris/Desktop/Sackler_PhD/ROTATIONS/ROTATION_2/Data/20140113_between_strain_cor_for_25_penotypes/Raw.data.Rfile")

length(Assembled) #to see the size of the corresponding list which should be 3

# From these we care about the first element of the list which is Assembled[[1]]

# Inspection of the element
head(Assembled[[1]])
tail(Assembled[[1]])

# We get rid of the columns we do not want
# namely "plate","replicate","condition"
# We do this as shown below

names(Assembled[[1]])
# Now we will keep the strain information as we want to calculate
# within strain correlation
myvars <- names(Assembled[[1]]) %in% c("plate", "replicate", "condition")
A_data_new <- Assembled[[1]][!myvars]

#now we are going to find how many unique strains we have:
uniq_strain_names <- unique(A_data_new$strain)

#By calculating the length of uniq_strain_names
#we can see how many unique strains there are
#in the dataset
length(uniq_strain_names)
[1] 392

# We generate the matrix that we will need to store the data.
# As there are 25 phenotypes and there will be one column
# for the correlation of phenotype 1 with each of the other
# phenotypes, the matrix will be 29*392, where strains=392
# and strain column + number of cells column + 27 phen. correlation
# columns = 29 columns in total

# I will make three individual ones that I will combine later
# one with all the unique strain names
# one for the correlations 
# one for the number of cells

strain_name_matrix <- matrix(rep(0,length(uniq_strain_names)),nrow=length(uniq_strain_names),ncol=1)
phen_corr_matrix <- matrix(rep(0,28*392), nrow=length(uniq_strain_names), ncol=28)
no_cells_matrix <- matrix(rep(0,length(uniq_strain_names)),nrow=length(uniq_strain_names),ncol=1)


for (i in 1:length(uniq_strain_names)) { 
	strain_data <- A_data_new[which(A_data_new$strain==uniq_strain_names[i]),]
	# Get the strain name 
	strain_name <- uniq_strain_names[i]
	# print (strain_name)
	# Get the correlation between phenotype 1 and all the others
	for (z in 2:29) {
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

# Check the dimensions of this matrix
dim(final_matrix)

# Calculate the medians for the correlations 
median_corr <- apply(phen_corr_matrix,2,median)
median_within_correlation <- median_corr

# final correlation data
final_corr_data <- cbind(corr_data,median_within_correlation) 
#to be sure that we have a dataframe as required by ggplot2
final_corr_data <- data.frame(final_corr_data)

# Now we get rid of the non-available values
# and more over we get the absolute values
# of the correlations

# to omit the non-available values
final_corr_data_2 <- na.omit(final_corr_data)

# to replace with absolute values
final_corr_data_2[,2] <- abs(as.numeric(final_corr_data_2[,2]))
final_corr_data_2[,3] <- abs(as.numeric(final_corr_data_2[,3]))

# to plot median_within_correlation with between_correlation
# with ggplot2
# First we load ggplot2
library(ggplot2)
> corr_plot <- ggplot(final_corr_data_2, aes(x=X2,y=median_within_correlation)) +
+ theme(text = element_text(size=10),axis.text.x = element_text(angle=90, vjust=1)) + xlab("Between correlation") +
+ ylab("Median within correlation") + geom_point(shape=2, color="black") + geom_text(aes(label=X1), size=1.5, hjust=1)
> corr_plot + geom_abline(intercept=0.01133664, slope=0.99517976, color="red") + annotate("text", label="r^2 == 0.9992", parse = TRUE, x=0.625, y=0.875, size=3)
> ggsave("final_between_within_plot.pdf",dpi=300)
> 
 