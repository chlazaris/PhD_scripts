# set the working directory

setwd("~/Desktop/20131213")

#read in the data
cell_data <- read.delim("CellData.csv", sep=",", stringsAsFactors=FALSE)

#inspect the data
head(cell_data)
tail(cell_data)

#To isolate the data that belong to F2_6
cell_data_F2_6 <- cell_data[which(cell_data$strain=="F2_6"),]

#To find out the total number of strains:
unique(cell_data$strain)
# The conclusion is that in this dataset, there are 225 strains

# The idea is to select two certain phenotypes, then
# go through all strains that we have available
# and report the strain name, calculate the correlation coefficient
# and the numbers of cells (of this strain) where the phenotypes
# where measured

# Thus the outcome will be a table:

# Strain #   r^2    #   No cells   
################################### 
#        #          #
#        #          #
#        #          #

#Let's start with a loop that gets from the cell_data, all
#the strain names and prints them

# This table is basically a matrix
corr_matrix <- matrix(rep(0,675),nrow=225,ncol=3)

# First we find the unique strains:
strains <- unique(cell_data$strain)

for (i in 1:length(strains)) {
	#print (i)
	#print (factor(strains[i]))
	#data_name <- paste(strains[i],"strain_data",sep="_") 
	strain_data <- cell_data[which(cell_data$strain==strains[i]),]
	# Get the strain name 
	strain_name <- strains[i]
	# Get the correlation between phenotype 1 and 2
	phen_cor <- cor(strain_data[,2],strain_data[,3])
	# Get the number of cells for strain
	no_cells <- length(strain_data[,1])
	# Now fill in the matrix
	# each row is a unique strain
	corr_matrix[i,1] <- strain_name
 	corr_matrix[i,2] <- phen_cor
 	corr_matrix[i,3] <- no_cells
} 

corr_data <- data.frame(corr_matrix)

# Get the first one as character (strain names)
strain_names <- as.character(levels(corr_data$X1))[corr_data$X1]

# Get the second column as numeric(corr coefficients)
corr_data_num <- as.numeric(levels(corr_data$X2))[corr_data$X2]

# Get the third one as numbers of cells
no_cells <-  as.numeric(levels(corr_data$X3))[corr_data$X3]

# Now make new data frame
P1_P2_data <- cbind(strain_names,corr_data_num,no_cells)
P1_P2_data <- data.frame(P1_P2_data)

# to make it numeric
corr_data_num <- as.numeric(levels(P1_P2_data$corr_data_num))[P1_P2_data$corr_data_num]

# produce .pdf files
pdf("P1_P2_within_strain_cor_1.pdf")
hist(corr_data_num, xlab="Within Strain Correlation", ylab="Frequency", main="P1_P2_within_strain_correlation", xlim=c(-1,1), ylim=c(0,100))
dev.off()

pdf("P1_P2_within_strain_cor_2.pdf")
hist(corr_data_num, xlab="Within Strain Correlation", ylab="Frequency", main="P1_P2_within_strain_correlation")
dev.off()



# To set the X coordinates
# coord_cartesian(xlim = c(-1, 1)) 
# To set the Y coordinates
# coord_cartesian(ylim = c(0, 100)) 



