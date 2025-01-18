#check if two nucleotide sequences are equal
A <- c("")
B <- c("")
all(A==B)



A <- c("sequence 1")
B <- c("sequence 2")
all(A==B)


#nucleotide pruning
y <- gsub("\n" , "", c("Enter the sequence"))
substr(y, start = 22 , stop = 720)


y <- gsub("\n" , "", c(""))
substr(y, start = 30 , stop = 652)
