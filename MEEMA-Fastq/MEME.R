meema_metadata <- read.csv("meema_metadata.csv")
summary(meema_metadata)
VO2max <- meema_metadata$VO2max
PercentFat <- meema_metadata$PercentFat
BMC_BMD <- na.omit(meema_metadata$BMC) / na.omit(meema_metadata$BMD)
