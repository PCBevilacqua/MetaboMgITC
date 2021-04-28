

####Make a vector of data frames to read####

file.path.vector = paste("Data", list.files("Data"), sep = "/")

data.list = lapply(file.path.vector, read.itc)

head(data.list[[1]]$inj)
head(data.list[[5]]$inj)

?nls

MetaboMgITC(data.list[[5]], data.list[[1]], Saturation.threshold = 0.8,
            Save.path.and.prefix = "Test")




