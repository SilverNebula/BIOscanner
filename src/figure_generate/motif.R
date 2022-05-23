library("ggseqlogo")
require("ggplot2")
pcm <- read.table("motif2.txt")
# pcm <- t(pcm)
pcm <- as.matrix(pcm)
# pcm <- pcm[,0:ncol(pcm)]
rownames(pcm) <- c("A", "G", "C", "T")
# motif <- table(pcm,rownames)
# plot(motif, ic.scale=FALSE, ylab="probability")
print(pcm)
g <- ggseqlogo(pcm, method = "custom") +
    theme(axis.text.x = element_blank()) + ggtitle("chr3:142740137-142740263")
png("seqlogo.png", width = 3600, height = 400, units = "px", res = 300)
g
dev.off()