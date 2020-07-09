#Allele Count = allele count in genotyprs for each AT allele
#Allele number = total number of alleles in called genotypes

# -----
# col1 is always european frequenncies, col2 is african or asian
# 1:212449403
# african
x1 <- matrix( c(9061, 15368, 7885, 8678), nrow = 2)
fisher.test(x1)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))

# east asian
x2 <- matrix( c(9061, 15368, 539, 1544), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 2:111927379
# african
x1 <- matrix( c(7577, 15368, 3101, 8690), nrow = 2)
fisher.test(x1)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))

# east asian
x2 <- matrix( c(7577, 15368, 673, 1548), nrow = 2)



fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 5:1280128
# african
x1 <- matrix( c(6463, 15388, 5778, 8680), nrow = 2)
fisher.test(x1)
fisher.test(x1)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))

# east asian
x2 <- matrix( c(6463, 15388, 607, 1554), nrow = 2)
fisher.test(x2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 6:32032421
# african
x1 <- matrix( c(1829, 15406, 1938, 8702), nrow = 2)

eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(1829, 15406, 195, 1552), nrow = 2)

fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 6:33533625
# african
x1 <- matrix( c(3918, 15340, 1992, 8660), nrow = 2)


eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)


# east asian
x2 <- matrix( c(3918, 15340, 64, 1554), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 8:120933963
# african
x1 <- matrix( c(6731, 15090, 5919, 8536), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(6731, 15090, 79, 1548), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 9:779507
# african
x1 <- matrix( c(8918, 15386, 2220, 8694), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(8918, 15386, 744, 1546), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 9:127190340
# african
# values from dbsnp 1kg
x1 <- matrix( c(334, 1006, 1148, 1322), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(334, 1006, 998, 1552), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 9:140073294
# african
x1 <- matrix( c(4095, 15384, 3242, 8668), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(4095, 15384, 94, 1556), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 10:7534248
# african
x1 <- matrix( c(5194, 15386, 953, 8690), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(5194, 15386, 898, 1556), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------


# ------
# 11:30351223
# african
x1 <- matrix( c(4497, 15382, 3786, 8688), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(4497, 15382, 898, 1554), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 11:1051495
# african
x1 <- matrix( c(2898, 15406, 254, 8714), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(2898, 15406, 116, 1560), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------




# ------
# 12:51301431
# african
x1 <- matrix( c(5183, 15398, 1112, 8686), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(5183, 15398, 371, 1554), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 12:53793209
# african
x1 <- matrix( c(2653, 15072, 1749, 8598), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(2653, 15072, 308, 1538), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------




# ------
# 17:76691564
# african
x1 <- matrix( c(8591, 15376, 1028, 8692), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(8591, 15376, 655, 1554), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 17:692095 
# african
x1 <- matrix( c(8674, 15358, 6890, 8694), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(8674, 15358, 1072, 1548), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 19:28356614
# african
x1 <- matrix( c(3109, 15410, 2760, 8670), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(3109, 15410, 381, 1550), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# 20:52197366
# african
x1 <- matrix( c(1986, 15188, 1833, 8692), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(1986, 15188, 110, 1556), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# X:24384181
# african
x1 <- matrix( c(1461, 10667, 336, 5778), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(1461, 10667, 117, 961), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------



# ------
# X:66489986
# african
x1 <- matrix( c(8878, 10699, 767, 5793), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(8878, 10699, 987, 989), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------




# ------
# X:100432681
# african
x1 <- matrix( c(4849, 10552, 3470, 5718), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(4849, 10552, 324, 969), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------




# ------
# X:153535143
# african
x1 <- matrix( c(2886, 10246, 3646, 5523), nrow = 2)
eur.maf <- c(x1[1,1] / x1[2,1])
afr.maf <- c(x1[1,2] / x1[2,2])
print(paste0("eur.maf = ", eur.maf, "; afr.maf = ", afr.maf))
fisher.test(x1)

# east asian
x2 <- matrix( c(2886, 10246, 385, 885), nrow = 2)
fisher.test(x2)
eur.maf <- c(x2[1,1] / x2[2,1])
asian.maf <- c(x2[1,2] / x2[2,2])
print(paste0("eur.maf = ", eur.maf, "; asian.maf = ", asian.maf))
# ------