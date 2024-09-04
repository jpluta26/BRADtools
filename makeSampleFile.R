# file containing phenotype definitions
fam <- read.table("pheno.fam", header = F)

# subject info
dat <- read.table("dat.fam", header = F)

# principle components
pc <- read.table("dat.eigenvec", header = F)
dat$V1 <- dat$V2
dat$V6 <- 0

# phenotype assignment example
dat$V6[ dat$V2 %in% fam$V2] <- 1
dat <- dat[,-3]
dat <- cbind(dat, pc[,3:12])


colnames(dat) <- c("ID_1", "ID_2", "missing","sex", "phenotype", "EV1", "EV2", "EV3", "EV4", "EV5", "EV6", "EV7", "EV8", "EV9", "EV10")

line <- c(0,0,0,"D","B","C", "C", "C", "C", "C", "C", "C", "C", "C", "C")
dat <- rbind(line,dat)

write.table(dat, "dat.sample", quote = FALSE, row.names = FALSE, col.names = TRUE)
