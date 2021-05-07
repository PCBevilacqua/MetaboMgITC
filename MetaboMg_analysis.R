####Load the package####

devtools::load_all()

####Make a vector of data frames to read####

file.path.vector = paste("Data", list.files("Data"), sep = "/")

data.list = lapply(file.path.vector, read.itc)

head(data.list[[1]]$inj)
head(data.list[[5]]$inj)

?read.itc

df = read.itc("Data/B15mMMgintoB25CpH7js2021.itc", print.peak.integration.graph = TRUE)


?MetaboMgITC

fit = MetaboMgITC(data.list[[5]], data.list[[1]], Saturation.threshold = 0.8,
            Save.path.and.prefix = "Test")

fit


