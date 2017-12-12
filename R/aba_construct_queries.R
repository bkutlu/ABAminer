# Band aid fix to query the ABA

library(tidyjson)
library(dplyr)
library(magrittr)
library(jsonlite)
library(tibble)
library(urltools)


# structures wanted
Pvhd
arc
vmh
lha

#
baseUrl <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56']"


structureCodes <- list(
  grey = "8",
  pvhd = "63",
  vmh = "693",
  lha = "194",
  arc = "223"
)




