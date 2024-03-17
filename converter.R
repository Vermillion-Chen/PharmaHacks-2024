getwd()
setwd("./Coding/PharmaHacks-2024") 

data = readRDS("./stephenson.rds")
data = UpdateSeuratObject(cts)

library(data.table)
data_to_write_out <- as.data.frame(as.matrix(seuratObject@scale.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "stephenson.csv")