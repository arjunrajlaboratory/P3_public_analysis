library(tidyverse)

workingDir = "~/Dropbox (RajLab)/Projects/cellid/miscellaneous/"

tallTableName = paste0(workingDir, "RNAtagTest4_conditionsPositions_tall.csv")
spreadTableName = paste0(workingDir, "RNAtagTest4_conditionsPositions_wide.csv")

stockPositions = as_tibble(read.table(file = paste0(workingDir, "RNAtagTest4_drugs_stockPositions.tsv"), sep = "\t", header = T)) # TBD

set.seed(629) # seed 629 for test4

# initialize positions
positions = paste0(rep(LETTERS[1:8], each = 12), rep("_", times = 96), rep(1:12, times = 8))
conditions = rep("empty", times = 96)
names(conditions) = positions

conditionsStocks = rep("empty", times = 96)
names(conditionsStocks) = positions

drugs = paste0("d", 34:58)
# stocks = stockPositions$Test3StocksPlatePosition

drugStocks = as.data.frame(list(
  drugIDforPlating = drugs)) %>%
  inner_join(stockPositions, by = "drugIDforPlating")
stocks = as.character(drugStocks$test4StocksPlatePosition)

# each in triplicate
drugConds = rep(drugs, times = 3)
stockConds = rep(stocks, times = 3)

controls = rep("DMSOonly-control", times = 17)
hela = rep("DMSOHeLa-HeLa", times = 4)


# indVec = integer(0)
# for (hp in helaPos) {
#   ind = which(names(conditions) == hp)
#   indVec = c(indVec, ind)
# }

positionsForShuffle = positions
# rawConds = c(drugConds, controls, hela)
rawConds = c(drugConds, controls, hela)
names(rawConds) = sample(positionsForShuffle) # randomize via seed above

# rawCondsStocks = c(stockConds, controls, hela)
rawCondsStocks = c(stockConds, controls, hela)
names(rawCondsStocks) = names(rawConds)

# conditions[helaPos] = hela
conditions[names(rawConds)] = rawConds

# conditionsStocks[helaPos] = hela
conditionsStocks[names(rawCondsStocks)] = rawCondsStocks


tabForPrinting = list()
tabForPrinting$platePosition = names(conditions)
tabForPrinting$condition = as.character(conditions)
tabForPrinting = as_tibble(as.data.frame(tabForPrinting))

write.csv(tabForPrinting, file = tallTableName, quote = F, row.names = F)

tabForPrintingSpread = tabForPrinting %>%
  separate(platePosition, into = c("row", "column"), sep = "\\_") %>%
  spread(column, condition)

tabForPrintingSpread = tabForPrintingSpread[,c(1,2,6:13,3,4,5)]

write.csv(tabForPrintingSpread, file = spreadTableName, quote = F, row.names = F)


tallTableName = paste0(workingDir, "RNAtagTest4_conditionsPositions_stockPositions_tall.csv")
spreadTableName = paste0(workingDir, "RNAtagTest4_conditionsPositions_stockPositions_wide.csv")

tabForPrintingStocks = list()
tabForPrintingStocks$platePosition = names(conditionsStocks)
tabForPrintingStocks$condition = as.character(conditionsStocks)
tabForPrintingStocks = as_tibble(as.data.frame(tabForPrintingStocks))

write.csv(tabForPrintingStocks, file = tallTableName, quote = F, row.names = F)

tabForPrintingStocksSpread = tabForPrintingStocks %>%
  separate(platePosition, into = c("row", "column"), sep = "\\_") %>%
  spread(column, condition)

tabForPrintingStocksSpread = tabForPrintingStocksSpread[,c(1,2,6:13,3,4,5)]

write.csv(tabForPrintingStocksSpread, file = spreadTableName, quote = F, row.names = F)


### Generic well-pool-tag layout
wellTagTab = list(
  rowLet = rep(LETTERS[1:8], times = 12),
  colNum = rep(1:12, each = 8),
  poolID = rep(c("test4-pool1", "test4-pool2", "test4-pool3"), each = 32), #change to poolID names for convenience
  tagID = rep(c(paste0("Tag0", as.character(1:9)), paste0("Tag", as.character(10:32))), times = 3)
)
wellTagTab = as.data.frame(wellTagTab)

write.csv(wellTagTab, file = paste0(projectDir, "metadata/RNAtagTest4//wellIDpoolTagMapping.csv"), quote = F, row.names = F)
