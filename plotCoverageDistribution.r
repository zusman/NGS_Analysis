#prints coverage distribution in a given coverage file and
#gives a coverage range

cmd_args = commandArgs()
data<-read.table(cmd_args[4], header=T)
pdf(cmd_args[5])
ref_length<-as.numeric(cmd_args[6])


library(ggplot2)
print(ref_length)

covered<-data[which(data$coverage>0),]

p <- ggplot(covered, aes(x=position, y=coverage, group=reference))
p + geom_line(aes(colour = reference)) + labs(x="position on reference genome (in nucleotides)", y="coverage (in number of mapped reads)", fill="")+scale_x_continuous(limits=c(1,ref_length))
dev.off()
q()