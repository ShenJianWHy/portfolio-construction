rm(list=ls())
for (file in list.files(path = "code/", full.names = TRUE)) source(file)
