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

# Calculation of median within and between correlation
# for the 2nd element of the list (small buds)


# Calculation of median within and between correlation
# for the 3rd element of the list (large buds)
