## This code creates a table with the less variant and high
#  expressed 500 genes in ecery time point



# Get the actual matrix with values

rlt <- assay(rld1)

# Compute row sums and invert 

rs <- 1 / rowSums(rlt)

# Compute row variances
rv <- rowVars(rlt)

# Get order of rows, higher sums are better, smaller variances are better

o<- order(rs, rv)

# Reorder table

rlt <- rlt[o, ]

# Look at the table

head(rlt)

# Subset for the first 500 genes

newrlt <- rlt[1:500,]

write.csv(newrlt, file = "E19_less_var.csv")

