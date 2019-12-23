

if("haven" %in% rownames(installed.packages()) == FALSE)
	{
		install.packages("haven")
	}

library(haven)

file <- read_sas("sasfile.sas7bdat")
write.table(file, "output.csv", col.names=T, row.names=F, append=F, quote=F, sep=",")
