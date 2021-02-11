library(tidyverse)

workingDir = "~/Dropbox (RajLab)/Projects/cellid/miscellaneous/"

tallTableName = paste0(workingDir, "RNAtagTest3_conditionsPositions_tall_d20-100xToDMSO.csv")
spreadTableName = paste0(workingDir, "RNAtagTest3_conditionsPositions_wide_d20-100xToDMSO.csv")

stockPositions = read.table(file = paste0(workingDir, "RNAtagTest3_drugs_stockPositions.tsv"), sep = "\t", header = T) # TBD

set.seed(876) # seed 876 for test3

# initialize positions
positions = paste0(rep(LETTERS[1:8], each = 12), rep("_", times = 96), rep(1:12, times = 8))
conditions = rep("empty", times = 96)
names(conditions) = positions

conditionsStocks = rep("empty", times = 96)
names(conditionsStocks) = positions

drugs = paste0("d", c(2,5,7,13:21,23,24))
# stocks = stockPositions$Test3StocksPlatePosition

drugStocks = as.data.frame(list(
  drugIDforPlating = drugs)) %>%
  inner_join(stockPositions)
stocks = drugStocks$Test3StocksPlatePosition

drugConds = rep(paste0(rep(drugs, times = 2), rep(c("-1x", "-20x", "-100x"), each = 14)), times = 2)
stockConds = rep(paste0(rep(stocks, times = 2), rep(c("-1x", "-20x", "-100x"), each = 14)), times = 2)

# adjust for allopurinol 100x; switch d20-100x to DMSO control
drugConds2 = drugConds
drugConds2[drugConds2 == "d20-100x"] = "DMSOonly-control"
stockConds2 = stockConds
stockConds2[drugConds2 == "DMSOonly-control"] = "DMSOonly-control"

controls = rep("DMSOonly-control", times = 8)
hela = rep("DMSOHeLa-HeLa", times = 4)


# indVec = integer(0)
# for (hp in helaPos) {
#   ind = which(names(conditions) == hp)
#   indVec = c(indVec, ind)
# }

positionsForShuffle = positions
# rawConds = c(drugConds, controls, hela)
rawConds = c(drugConds2, controls, hela)
names(rawConds) = sample(positionsForShuffle) # randomize via seed above

# rawCondsStocks = c(stockConds, controls, hela)
rawCondsStocks = c(stockConds2, controls, hela)
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


tallTableName = paste0(workingDir, "RNAtagTest3_conditionsPositions_stockPositions_tall_d20-100xToDMSO.csv")
spreadTableName = paste0(workingDir, "RNAtagTest3_conditionsPositions_stockPositions_wide_d20-100xToDMSO.csv")

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
  poolID = rep(c("test2-pool1", "test2-pool2", "test2-pool3"), each = 32), #change to poolID names for convenience
  tagID = rep(c(paste0("Tag0", as.character(1:9)), paste0("Tag", as.character(10:32))), times = 3)
)
wellTagTab = as.data.frame(wellTagTab)

write.csv(wellTagTab, file = paste0(projectDir, "miscellaneous/wellIDpoolTagMapping.csv"), quote = F, row.names = F)
