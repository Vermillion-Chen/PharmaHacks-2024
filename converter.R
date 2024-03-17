getwd()
setwd("./Coding/PharmaHacks-2024") 

data = readRDS("./reduced.rds")
data = UpdateSeuratObject(data)
str(data[["data"]])
library(data.table)
data_to_write_out <- as.data.frame(as.matrix(seuratObject@scale.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "./reduced.csv")