library(tidyverse)

concTab = as_tibble(read.table("~/Dropbox (RajLab)/Projects/cellid/metadata/RNAtagTest6/20171106_RNAtagTest6_drugStockConcentrations - Sheet1.tsv", sep = "\t", header = T, stringsAsFactors = F))
pairTab = as_tibble(read.table("~/Dropbox (RajLab)/Projects/cellid/miscellaneous/pickingDrugs/20171106_drugPairsWithPairID.txt", sep = "\t", header = T, stringsAsFactors = F))
concTab$deliveredWell = "incomplete"

concTab2 = pairTab %>%
  mutate(Product.Name = drug1) %>%
  dplyr::select(Product.Name, pairID, stockWellID) %>%
  inner_join(concTab, by = "Product.Name")

concTab3 = pairTab %>%
  mutate(Product.Name = drug2) %>%
  dplyr::select(Product.Name, pairID, stockWellID) %>%
  inner_join(concTab, by = "Product.Name")

concTabF = bind_rows(concTab2, concTab3)

write.table(concTabF, file = "~/Dropbox (RajLab)/Projects/cellid/metadata/RNAtagTest6/RNAtagTest6_drugStockConcentrationWithPairPositions.txt", sep = "\t", quote = F, row.names = F)
