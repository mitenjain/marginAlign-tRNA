#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

#thanks SO: http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
raw1 <- strsplit(readLines(args[1]), "[[:space:]]+")
data1 <- lapply(raw1, tail, n = -1)
names(data1) <- lapply(raw1, head, n = 1)
data1 <- lapply(data1, as.numeric)

pdf(args[2])


library(lattice)

data <- list()
data$MappedReadLengths <- data1$ValuesReadLength
data$MappedReadCoverage <- data1$ValuesReadCoverage
data$MismatchesPerReadBase <- data1$ValuesMismatchesPerAlignedBase
data$ReadIdentity <- data1$ValuesReadIdentity
data$AlignmentIdentity <- data1$ValuesAlignmentIdentity
data$DeletionsPerBase <- data1$ValuesDeletionsPerReadBase
data$InsertionsPerBase <- data1$ValuesInsertionsPerReadBase

p1 <- smoothScatter(data$AlignmentIdentity~data$MappedReadLengths, main=list("Alignment Identity vs. Read Length"), xlab="Read Length (bases)", ylab="Read Identity", xlim=c(0,5000), ylim=c(0.40, 1.0), cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
p2 <- hist(data$MappedReadLengths, breaks="FD", xlim=c(0,5000), xlab="Read Length (bases)", ylab="Frequency", main="Read Length Histogram", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
p3 <- hist(data$MappedReadCoverage, xlab="Read Coverage", main="Read Coverage Histogram", breaks="FD", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)

alignedReadLength <- data$MappedReadCoverage * data$MappedReadLengths
p4 <- hist(alignedReadLength, breaks="FD", xlim=c(0,5000), xlab="Aligned Read Length", ylab="Frequency", main="Aligned Read Length Histogram", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)

#print the graphs
print(p1)
print(p3)
    
dev.off()
