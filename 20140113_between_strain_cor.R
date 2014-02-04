# 20130113_Calculate the between-strain correlation using the Raw.data.Rfile
# The file contains a list which is called "Assembled". By typing
# length(Assembled) I should get 3 (which corresponds to the 3 cell types
#(A,A1B and C respectively. A are the most cells with the least number 
#of phenotypes, A1B are the cells with small buds and C the ones with
#large buds)).

#setwd("~/Desktop/Sackler_PhD/ROTATIONS/ROTATION_2/Data/20140113_between_strain_cor_for_25_penotypes")
#raw_data <- read.delim("Raw.data.Rfile",sep=",",header=T)

load("/Users/chlazaris/Desktop/Sackler_PhD/ROTATIONS/ROTATION_2/Data/20140113_between_strain_cor_for_25_penotypes/Raw.data.Rfile")

length(Assembled) #to see the size of the corresponding list which should be 3

# From these we care about the first element of the list which is Assembled[[1]]

# Inspection of the element
head(Assembled[[1]])
tail(Assembled[[1]])

# We get rid of the columns we do not want
# namely "plate","replicate","condition","strain"
# We do this as shown below

names(Assembled[[1]])
myvars <- names(Assembled[[1]]) %in% c("plate", "replicate", "condition","strain")
A_data <- Assembled[[1]][!myvars]

#Then calculate all the coefficients and store them
#in a matrix (the correlation coefficients for between
# correlation).

# Create the correlation matrix
corr_matrix <- matrix(rep(0,length(names(A_data))), ncol=1)

# Create a loop to calculate the
# between-strain correlations of every 
# pair of P1 with the other phenotypes

for (i in 1:length(names(A_data))) {
	corr_matrix[i] <- cor(A_data[1],A_data[i])
}

# The results are given below:
> corr_matrix
              [,1]
 [1,]  1.000000000
 [2,]  0.986750690
 [3,] -0.408753725
 [4,]  0.948534634
 [5,]  0.944500852
 [6,] -0.074810588
 [7,]  0.248291518
 [8,]  0.149856917
 [9,]  0.425600214
[10,]  0.361278099
[11,]  0.034858762
[12,]  0.002914680
[13,]  0.847280392
[14,]  0.089145145
[15,]  0.376165392
[16,]  0.844318412
[17,]  0.404765814
[18,]           NA
[19,]           NA
[20,]           NA
[21,]           NA
[22,]  0.442660356
[23,]  0.443634944
[24,]  0.426496044
[25,]           NA
[26,]  0.216022931
[27,] -0.004382141
[28,]  0.173783614

# Save all the values for between correlation in a vector
between_correlation <- corr_matrix[,1]

# to remove the non-available values
corr_matrix_new <- corr_matrix[!is.na(corr_matrix)]

length(corr_matrix_new)
[1] 23

# All this work has been saved in the file: "20140113_within_strain_correlation_A"

# 20140114 

# Modify the matrix that is shown above in order to
# include the pairs of phenotypes for which the between
# correlation was calculated:

# The new matrix must have two columns
# size: 28 * 2

corr_matrix_2 <- matrix(rep(0,2*length(names(A_data))), ncol=2)

# Create a loop to report the pair
# of phenotypes (column 1) and calculate the
# between-strain correlations of every 
# pair of P1 with the other phenotypes (column 2)

for (i in 1:length(names(A_data))) {
	corr_matrix_2[i,1] <- paste(names(A_data)[1], names(A_data)[i],sep=":")
	corr_matrix_2[i,2] <- cor(A_data[1],A_data[i])
}

# Now we turn the matrix to dataframe and we have the following:

corr_data <- data.frame(corr_matrix_2)

corr_data

               X1                   X2
1  C11.1_A:C11.1_A                    1
2  C11.1_A:C12.1_A    0.986750689912998
3    C11.1_A:C13_A   -0.408753725029642
4   C11.1_A:C103_A    0.948534633938131
5   C11.1_A:C104_A     0.94450085198142
6   C11.1_A:C115_A  -0.0748105877199252
7   C11.1_A:C126_A    0.248291518309769
8   C11.1_A:C127_A    0.149856916836651
9  C11.1_A:D14.1_A    0.425600213843058
10 C11.1_A:D15.1_A    0.361278099093148
11 C11.1_A:D16.1_A   0.0348587618820034
12 C11.1_A:D17.1_A  0.00291467987577507
13  C11.1_A:D102_A    0.847280391513621
14  C11.1_A:D105_A   0.0891451445636542
15  C11.1_A:D117_A    0.376165391559611
16  C11.1_A:D127_A    0.844318411995605
17  C11.1_A:D135_A    0.404765813841391
18  C11.1_A:D147_A                 <NA>
19  C11.1_A:D148_A                 <NA>
20  C11.1_A:D154_A                 <NA>
21  C11.1_A:D155_A                 <NA>
22  C11.1_A:D173_A    0.442660356224744
23  C11.1_A:D176_A    0.443634943561417
24  C11.1_A:D179_A    0.426496044173185
25  C11.1_A:D182_A                 <NA>
26  C11.1_A:D188_A    0.216022930770395
27  C11.1_A:D191_A -0.00438214098503342
28  C11.1_A:D194_A    0.173783614138688

# Now remove the non-available values
# In order to retain the structure
# of the data frame we will use na.omit

corr_data_new <- na.omit(corr_data)
corr_data_new

                X1                   X2
1  C11.1_A:C11.1_A                    1
2  C11.1_A:C12.1_A    0.986750689912998
3    C11.1_A:C13_A   -0.408753725029642
4   C11.1_A:C103_A    0.948534633938131
5   C11.1_A:C104_A     0.94450085198142
6   C11.1_A:C115_A  -0.0748105877199252
7   C11.1_A:C126_A    0.248291518309769
8   C11.1_A:C127_A    0.149856916836651
9  C11.1_A:D14.1_A    0.425600213843058
10 C11.1_A:D15.1_A    0.361278099093148
11 C11.1_A:D16.1_A   0.0348587618820034
12 C11.1_A:D17.1_A  0.00291467987577507
13  C11.1_A:D102_A    0.847280391513621
14  C11.1_A:D105_A   0.0891451445636542
15  C11.1_A:D117_A    0.376165391559611
16  C11.1_A:D127_A    0.844318411995605
17  C11.1_A:D135_A    0.404765813841391
22  C11.1_A:D173_A    0.442660356224744
23  C11.1_A:D176_A    0.443634943561417
24  C11.1_A:D179_A    0.426496044173185
26  C11.1_A:D188_A    0.216022930770395
27  C11.1_A:D191_A -0.00438214098503342
28  C11.1_A:D194_A    0.173783614138688

#corr_data_new <- corr_data[!is.na(corr_data)]

# Create a histogram to visualise the relationship
# of between correlation with the phenotypes 

library(ggplot2)
m <- ggplot(corr_data_new, aes(x=X1,y=X2)) +
    theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, vjust=1),
		axis.text.y = element_blank(),
		axis.ticks = element_blank()) + xlab("Phenotypes") + ylab("Between Correlation")
m + geom_histogram()
ggsave("Pheno_vs_between_cor.pdf",dpi=300)





