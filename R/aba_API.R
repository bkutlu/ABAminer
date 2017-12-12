# https://github.com/sailthru/tidyjson
# devtools::install_github("sailthru/tidyjson")
library(tidyjson)
library(dplyr)
library(magrittr)
library(jsonlite)
library(tibble)
library(urltools)


#### PVHd vs grey
# http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq29]
# json
pvhdVsgrey_1 <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'4,50'][start_row$eq0][num_rows$eq20]"
pvhdVsgrey_2 <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'1,50'][start_row$eq2001][num_rows$eq4000]"

#### PvHD+ HY vs grey
# json
hyPvhdVsGrey <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'63,1097'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

#### Hy vs grey
hyVsgrey <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'1097'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

#### Hy vs pvhd
pvhdVsHy <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'1097'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"
hyVspvhd <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'63'][domain1_threshold$eq'0,50'][domain2$eq'1097'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

#### parse pvhd vs grey
pvhdVsgreyDf1 <- extractISH_de(pvhdVsgrey_1)
pvhdVsgreyDf1R <- mutate(pvhdVsgreyDf1, target2contrast = as.numeric(target_sum)/(as.numeric(contrast_sum)/10^5))
pvhdVsgreyDf1RSort <- arrange(pvhdVsgreyDf1R, as.numeric(target_sum), as.numeric(fold_change))
pvhdVsgreyDf1RSort <- arrange(pvhdVsgreyDf1R, desc(as.numeric(target_sum)),
                              desc(as.numeric(fold_change)),
                              desc(as.numeric(target2contrast)))

pvhdVsgreyDf1RSort2 <- arrange(pvhdVsgreyDf1R, desc(as.numeric(target_sum)),
                              desc(as.numeric(target2contrast)),
                              desc(as.numeric(fold_change)))

pvhdVsgreyDf1RSort3 <- arrange(pvhdVsgreyDf1R,
                               desc(as.numeric(target2contrast)),
                               desc(as.numeric(target_sum)),
                               desc(as.numeric(fold_change)))


# second batch
pvhdVsgreyDf2 <- extractISH_de(pvhdVsgrey_2)

hist(as.numeric(pvhdVsgreyDf$target_sum))

## HY
hyVsgreyDf <- extractISH_de(hyVsgrey)

## HY+PVHD
hyPvhdVsGreyDf <- extractISH_de(hyPvhdVsGrey)

## Hy vs PVHD
hyVspvhdDf <- extractISH_de(hyVspvhd)
pvhdVsHyDf <- extractISH_de(pvhdVsHy)

hy_pvhd_merge <- merge(hyVspvhdDf,pvhdVsHyDf, by = 'object_id', all = F)


# comparing the target sum values between different structures
dd = merge(hyVsgreyDf, pvhdVsgreyDf1, by = 'object_id', all = F)
mutate(dd, foldDiff = target_sum.x/target_sum.y)


extractISH_de <- function (myApiCall){
  tmpChar <- as.character(toJSON(fromJSON(myApiCall)))
  tmpChar %>% enter_object("msg") %>% gather_array() %>%
    spread_values(
      object_id = jstring("id"),
      atlas_gene_id = jstring("gene-id"),
      entrez = jstring("entrez-id"),
      gene_symbol = jstring("gene-symbol"),
      section = jstring("plane-of-section"),
      fold_change = jstring("fold-change"),
      target_sum = jstring("target-sum"),
      contrast_sum = jstring("contrast-sum")
    ) %>%
    select(object_id, atlas_gene_id, entrez, gene_symbol, section, fold_change, target_sum, contrast_sum) -> pvhdVsgreyDf

}#extractRowsFromDE

## documentation
# must read for rules
http://help.brain-map.org/display/api/RESTful+Model+Access+(RMA)
# structures
http://api.brain-map.org/api/v2/structure_graph_download/1.json
# fold change
http://api.brain-map.org/examples/foldchange/index.html

http://help.brain-map.org/download/attachments/2818169/FineStructureAnnotation.pdf?version=1&modificationDate=1319477154436


##### RMA query builder
http://www.brain-map.org/api/examples/examples/rma_builder/index.htmlhttp://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'4,50'][start_row$eq0][num_rows$eq20]



##############################################
######### experimental work ################
###############################################
# simple example - cat
simplJEx <- fromJSON("data/simplJsonExample.json")
simplJExChar <- as.character(toJSON(simplJEx))

simplJExChar %>% prettify()

simplJExChar %>% gather_keys%>% gather_array() %>% json_types()

simplJExChar %>% gather_keys%>% gather_array() %>%


  # enter_object("children")
  # spread_values(name = jstring("name"))
  # enter_object("name")

amnts <- simplJExChar %>% spread_values(name = jstring("name")) %>%
  enter_object("children") %>% gather_array() %>%
  spread_values(companyType = jstring("name")) %>%
  enter_object("children") %>% gather_array() %>%
  spread_values(
    companyName = jstring("name"),
    cost = jnumber("cost"),
    pVal = jnumber("p")
  ) %>%
   mutate(totalCost = cost / 10^5) %>%
select (companyType, companyName, totalCost, pVal)


amntSum <- amnts %>% group_by(companyType) %>%
  summarize(aggreg = sum(totalCost))

# tally to count the numbers

# attack the real world example
# try with smaller sample
abaSimpl <- fromJSON("data/abaDownExample.json")

# exploring the json onbj
# somehow you have to convert to JSON then character
abaSimplChar <- as.character(toJSON(abaSimpl))
abaSimplChar %>% prettify
abaSimplChar %>% gather_keys %>% json_types
abaSimplChar %>% gather_keys %>% json_types %>% json_lengths



abaSimplChar %>% enter_object("children") %>%
  gather_array() %>% gather_keys()






##### some ideas that did not work
# # Fetch JSON data from ABA
#
# mouseStructureGraph <- fromJSON("http://api.brain-map.org/api/v2/structure_graph_download/1.json")
# hadi <- flatten(mouseStructureGraph[[6]], recursive = TRUE)
#
#
# # library(rjson)
# # mouseStructureGraph2 <- rjson::fromJSON(readLines("http://api.brain-map.org/api/v2/structure_graph_download/1.json")[1])
#
#
# library(RJSONIO)
# mouseStructureGraph2 <- RJSONIO::fromJSON(readLines("http://api.brain-map.org/api/v2/structure_graph_download/1.json")[1])

# abaStr <- fromJSON("http://api.brain-map.org/api/v2/structure_graph_download/1.json")
#
# # convert to JSON then to tbl_json
# abaStrJson <- toJSON(abaStr[['msg']])
# abaStrTbl <- abaStrJson %>% as.tbl_json
#
# # convert the df to JSON
# z <- abaStr[['msg']] %>% as.tbl_json(json.column = 'children')
#
# abaStrJson %>% gather_array
#
#
# abaStr <- fromJSON("http://api.brain-map.org/api/v2/structure_graph_download/1.json")
#
# abaRaw <- readLines("http://api.brain-map.org/api/v2/structure_graph_download/1.json", skipNul = TRUE)
# abaRaw %>% as.tbl_json # does not parse


http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,
structure_sets[id$eq2],pipe::list[xstructures$eq'id'],model::Gene[id$eq18311],
rma::include,organism,chromosome,gene_aliases,data_sets(products[id$eq1]),
pipe::list[xdataset$eq'data_sets/*/id'],model::Gene[id$eq18311],
rma::include,organism,chromosome,gene_aliases,data_sets(products[id$eq1]),
model::SectionDataSet[id$in$xdataset],rma::include,genes,plane_of_section,treatments,specimen(donor(age,organism)),
probes(orientation,predicted_sequence,forward_primer_sequence,reverse_primer_sequence),
products[id$eq1],model::StructureUnionize,rma::criteria,section_data_set[id$in$xdataset],
rma::include,structure[id$in$xstructures],rma::options[only$eq'id,name,expression_energy,acronym,red,green,blue'],
model::SectionImage[data_set_id$in$xdataset],
rma::include,associates,alternate_images,rma::options[order$eq'sub_images.section_number$asc'],


# pvhd vs others
http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,structure_sets%5Bid$eq2%5D,rma::options%5Bonly$eq%27id%27%5D,pipe::list%5Bxstructures$eq%27id%27%5D,service::differential_rows%5Bset$eq%27P56%27%5D%5Bdomain1$eq%271129,567,157,141,88,331,515,980,1004,693,946,290,10671,313,1065,512,549%27%5D%5Bdomain1_threshold$eq%270,50%27%5D%5Bdomain2$eq%2763%27%5D%5Bdomain2_threshold$eq%271,50%27%5D%5Bstart_row$eq21%5D%5Bnum_rows$eq21%5D




arhVnucleus <- 'http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,structure_sets%5Bid$eq2%5D,rma::options%5Bonly$eq%27id%27%5D,pipe::list%5Bxstructures$eq%27id%27%5D,service::differential_rows%5Bset$eq%27P56%27%5D%5Bdomain1$eq%278%27%5D%5Bdomain1_threshold$eq%270,50%27%5D%5Bdomain2$eq%27223%27%5D%5Bdomain2_threshold$eq%271,50%27%5D%5Bstart_row$eq0%5D%5Bnum_rows$eq29%5D'

# arhVnucleusDc <- url_decode(arhVnucleus)
arhVnucleusDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'223'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"
arhVnucleusVsHypoDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'1097'][domain1_threshold$eq'0,50'][domain2$eq'223'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

arcVsgrey <- extractISH_de(arhVnucleusDc)
arcVsHy <- extractISH_de(arhVnucleusVsHypoDc)


arcSpec <- merge(arcVsgrey, arcVsHy, by = 'object_id', all = F)
save(arcSpec, file = "C:/Users/bukt/Downloads/arcSPec.RData")

# Fab 28 2017
# 1. Find out how they sort
arcSpecDf <- arcSpec %>%
  as_data_frame %>%
  mutate_each(funs(as.numeric), matches("entrez*|fold*|atlas*|contrast*|target*")) %>%
  mutate(originalOrder = 1:nrow(.))




# 2. Design and implement





# convert secreted and tm from human to mouse
# secreted and tm
library(DBI)
library(RSQLite)
sqlite <- dbDriver("SQLite")
atlasdb <- dbConnect(sqlite,"C:/Users/bukt/My Documents/DataResources/atlasLists.db") # makes a new file
dbListTables(atlasdb)
# test
dbReadTable(atlasdb,"geneAnnots")[1:3,]
# only human
sql <- paste0("SELECT *
                FROM atlasListTable alt
                LEFT JOIN listTable lt
                ON alt.listTableId = lt.listTableId
                WHERE
                lt.lists = 'hcomb_secreted'
                OR
                lt.lists = 'hcomb_pmreceptor'")
# with mouse annotations
sql2 <- paste0("SELECT *
                FROM atlasListTable alt
                LEFT JOIN listTable lt
                ON alt.listTableId = lt.listTableId
                LEFT JOIN orthologs otl
                ON alt.ensembl_id = otl.ensembl_gene_id
                WHERE
                lt.lists = 'hcomb_secreted'
                OR
                lt.lists = 'hcomb_pmreceptor'")

# with mouse annotations
sql3 <- paste0("SELECT *
                FROM atlasListTable alt
                LEFT JOIN listTable lt
                ON alt.listTableId = lt.listTableId
                LEFT JOIN orthologs otl
                ON alt.ensembl_id = otl.ensembl_gene_id
                LEFT JOIN geneAnnots ga
                ON otl.mmusculus_homolog_ensembl_gene = ga.ensembl_id
                WHERE
                lt.lists = 'hcomb_secreted'
                OR
                lt.lists = 'hcomb_pmreceptor'")


# all the extracellular gene list
extracellGenes <- fetch(dbSendQuery(atlasdb,sql), n = -1)
extracellGenesWithOrthologs <- fetch(dbSendQuery(atlasdb,sql2), n = -1)
# myEnsemblGenematrix_ExtC <- merge(myMatrixWithEnsIds,extracellGenes, by.x = 'ensembl_id', by.y = 'ensembl_id', all.y = F)

extracellGenesWithOrthologsGeneIds <- fetch(dbSendQuery(atlasdb,sql3), n = -1)


# merge
arcSpecExt <- merge(arcSpec, extracellGenesWithOrthologsGeneIds, by.x = 'entrez.x', by.y ='gene_id')

library(xlsx)
write.xlsx(arcSpecExt, file = 'results/extracellularSpecific.xlsx', row.names = F, sheetName = "ARC", append = TRUE)

# create function
# arguments
# target structure
# contrast structure (usually within the same mother structure)
# extracellular


# TO DO
# lookup the structure id
# get the chunks of 2000




#### VMH
vmhUrl <- url_decode("http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,structure_sets%5Bid$eq2%5D,rma::options%5Bonly$eq%27id%27%5D,pipe::list%5Bxstructures$eq%27id%27%5D,service::differential_rows%5Bset$eq%27P56%27%5D%5Bdomain1$eq%278%27%5D%5Bdomain1_threshold$eq%270,50%27%5D%5Bdomain2$eq%27693%27%5D%5Bdomain2_threshold$eq%271,50%27%5D%5Bstart_row$eq0%5D%5Bnum_rows$eq28%5D")
vmhDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'693'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"
vmhVsHypoDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'1097'][domain1_threshold$eq'0,50'][domain2$eq'693'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

vmhVsgrey <- extractISH_de(vmhDc)
vmhVsHy <- extractISH_de(vmhVsHypoDc)


vmhSpec <- merge(vmhVsgrey, vmhVsHy, by = 'object_id')
vmhSpecExt <- merge(vmhSpec, extracellGenesWithOrthologsGeneIds, by.x = 'entrez.x', by.y ='gene_id')
write.xlsx(vmhSpecExt, file = 'results/extracellularSpecific.xlsx', row.names = F, sheetName = "VMH", append = TRUE)


### vmh refined search
baseVmhUrl1 <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56']"
baseUrl2 <- "[domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"
baseUrl2_2 <- "[domain2_threshold$eq'1,50'][start_row$eq2001][num_rows$eq4000]"
baseUrl2_3 <- "[domain2_threshold$eq'1,50'][start_row$eq4001][num_rows$eq6000]"
domain1_Contrast <- "[domain1$eq'567,549,313,1065,512']"
domain2_target <- "[domain2$eq'693']"

vmhVsRest <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2))
vmhVsRest2 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_2))
vmhVsRest3 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_3))

vmhVsRestDat <- extractISH_de(vmhVsRest)
vmhVsRestDat_2 <- extractISH_de(vmhVsRest2)
vmhVsRestDat_3 <- extractISH_de(vmhVsRest3)

vmhVsRestDat_all <- rbind(vmhVsRestDat, vmhVsRestDat_2, vmhVsRestDat_3)
vmhVsRestDat_all <- subset(vmhVsRestDat_all, fold_change >= 1)


# # vmh vs rest of hypothalamus
# vmhVsAHN_MBO_PH_PVR_PVZ_LZ_ME <- "[domain1$eq'88,331,946,141,157,290,10671']"
#
# vmhVsAHN_MBO_PH_PVR_PVZ_LZ_ME2k <- extractISH_de(URLencode(paste0(baseVmhUrl1,vmhVsAHN_MBO_PH_PVR_PVZ_LZ_ME,domain2_target, baseUrl2)))
#

# Compare fold changes of different comparisons
# library(ggplot2)
# vmhGreyVsAllOther <- merge(vmhVsgrey, vmhVsChThMbHbCb2k, by = 'object_id')
# qplot(as.numeric(vmhGreyVsAllOther[,"fold_change.x"]), as.numeric(vmhGreyVsAllOther[,"fold_change.y"]), xlim = c(-1, 175), ylim = c(-1, 175), xlab = "VMH vs Grey", ylab = "VMH vs All others", main = "Fold Enrichment of VMH specific transcripts") + geom_abline()

vmhVsRestExt <- merge(vmhVsRestDat_all, extracellGenesWithOrthologsGeneIds, by.x = 'entrez', by.y ='gene_id')

library(xlsx)
library(dplyr)
# write.xlsx(vmhVsChThMbHbCb2kExt, file = 'results/extracellularSpecificVmh.xlsx', row.names = F, sheetName = "VMH", append = TRUE)
vmh_sec_memb <- vmhVsRestExt
vmh_sec_memb <- vmh_sec_memb[,c("object_id","gene_symbol", "section", "fold_change","lists","ensembl_id.1")]
colnames(vmh_sec_memb) <- c("object_id","gene_symbol", "section", "fold_change","Location","ensembl_id")
save(vmh_sec_memb, file = "results/vmh_sec_memb.RData")

vmh_sec_membDf <- as_data_frame(vmh_sec_memb)
save(vmh_sec_membDf, file = "results/vmh_sec_membDf.RData")


### vmh modified 3-20-2017
# add more rows
# merge the results
# save
# study rentrez for inspiration on how to make calls



# PVHD updated
domain2_target <- "[domain2$eq'63']"
pvhdVsR <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2))
pvhdVsR_2 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_2))
pvhdVsR_3 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_3))


pvhdVsRD <-  extractISH_de(pvhdVsR)
pvhdVsR_2 <-  extractISH_de(pvhdVsR_2)
pvhdVsR_3 <-  extractISH_de(pvhdVsR_3)

pvhdVsRAll <- rbind(pvhdVsRD, pvhdVsR_2, pvhdVsR_3)
pvhdVsRAll <- subset(pvhdVsRAll, fold_change >= 1)


pvhdVsRAllExt <- merge(pvhdVsRAll, extracellGenesWithOrthologsGeneIds, by.x = 'entrez', by.y ='gene_id')
pvhd_sec_memb <- pvhdVsRAllExt
pvhd_sec_memb <- pvhd_sec_memb[,c("object_id","gene_symbol", "section", "fold_change","lists","ensembl_id.1")]
colnames(pvhd_sec_memb) <- c("object_id","gene_symbol", "section", "fold_change","Location","ensembl_id")
save(pvhd_sec_memb, file = "results/pvhd_sec_memb.RData")

pvhd_sec_membDf <- as_data_frame(pvhd_sec_memb)
save(pvhd_sec_membDf, file = "results/pvhd_sec_membDf.RData")

# ARC
domain2_target <- "[domain2$eq'223']"
arcVsR <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2))
arcVsR_2 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_2))
arcVsR_3 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_3))


arcVsRD <-  extractISH_de(arcVsR)
arcVsR_2 <-  extractISH_de(arcVsR_2)
arcVsR_3 <-  extractISH_de(arcVsR_3)

arcVsRAll <- rbind(arcVsRD, arcVsR_2, arcVsR_3)
arcVsRAll <- subset(arcVsRAll, fold_change >= 1)


arcVsRAllExt <- merge(arcVsRAll, extracellGenesWithOrthologsGeneIds, by.x = 'entrez', by.y ='gene_id')
arc_sec_memb <- arcVsRAllExt
arc_sec_memb <- arc_sec_memb[,c("object_id","gene_symbol", "section", "fold_change","lists","ensembl_id.1")]
colnames(arc_sec_memb) <- c("object_id","gene_symbol", "section", "fold_change","Location","ensembl_id")
save(arc_sec_memb, file = "results/arc_sec_memb.RData")

# LHA
domain2_target <- "[domain2$eq'194']"
lhaVsR <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2))
lhaVsR_2 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_2))
lhaVsR_3 <- URLencode(paste0(baseVmhUrl1,domain1_Contrast,domain2_target, baseUrl2_3))


lhaVsRD <-  extractISH_de(lhaVsR)
lhaVsR_2 <-  extractISH_de(lhaVsR_2)
lhaVsR_3 <-  extractISH_de(lhaVsR_3)

lhaVsRAll <- rbind(lhaVsRD, lhaVsR_2, lhaVsR_3)
lhaVsRAll <- subset(lhaVsRAll, fold_change >= 1)


lhaVsRAllExt <- merge(lhaVsRAll, extracellGenesWithOrthologsGeneIds, by.x = 'entrez', by.y ='gene_id')
lha_sec_memb <- lhaVsRAllExt
lha_sec_memb <- lha_sec_memb[,c("object_id","gene_symbol", "section", "fold_change","lists","ensembl_id.1")]
colnames(lha_sec_memb) <- c("object_id","gene_symbol", "section", "fold_change","Location","ensembl_id")
save(lha_sec_memb, file = "results/lha_sec_memb.RData")




# PVHd
pvhdUrl <- url_decode("http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,structure_sets%5Bid$eq2%5D,rma::options%5Bonly$eq%27id%27%5D,pipe::list%5Bxstructures$eq%27id%27%5D,service::differential_rows%5Bset$eq%27P56%27%5D%5Bdomain1$eq%278%27%5D%5Bdomain1_threshold$eq%270,50%27%5D%5Bdomain2$eq%2763%27%5D%5Bdomain2_threshold$eq%271,50%27%5D%5Bstart_row$eq0%5D%5Bnum_rows$eq28%5D")
pvhdDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"
pvhdVsHypoDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'1097'][domain1_threshold$eq'0,50'][domain2$eq'63'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

pvhdVsgrey <- extractISH_de(pvhdDc)
pvhdVsHy <- extractISH_de(pvhdVsHypoDc)


pvhdSpec <- merge(pvhdVsgrey, pvhdVsHy, by = 'object_id')
pvhdSpecExt <- merge(pvhdSpec, extracellGenesWithOrthologsGeneIds, by.x = 'entrez.x', by.y ='gene_id')
# what if you keep all
pvhdSpecAll <- merge(pvhdVsgrey, pvhdVsHy, by = 'object_id', all = T)
# create another column with non-NA values
nonNaEG <- apply(pvhdSpecAll, 1, function (s){
  na.exclude(unique(c(s['entrez.x'],s['entrez.y'])))
})
pvhdSpecAll <- cbind(pvhdSpecAll, nonNaEG)
pvhdSpecAllExt <-  merge(pvhdSpecAll, extracellGenesWithOrthologsGeneIds, by.x = 'nonNaEG', by.y ='gene_id')


write.xlsx(pvhdSpecExt, file = 'results/extracellularSpecific.xlsx', row.names = F, sheetName = "PVHd", append = TRUE)


# LHA
lhaUrl <- url_decode("http://mouse.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,structure_sets%5Bid$eq2%5D,rma::options%5Bonly$eq'id'%5D,pipe::list%5Bxstructures$eq'id'%5D,service::differential_rows%5Bset$eq'P56'%5D%5Bdomain1$eq'8'%5D%5Bdomain1_threshold$eq'0,50'%5D%5Bdomain2$eq'194'%5D%5Bdomain2_threshold$eq'1,50'%5D%5Bstart_row$eq0%5D%5Bnum_rows$eq20%5D")
lhaDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'8'][domain1_threshold$eq'0,50'][domain2$eq'194'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"
lhaVsHypoDc <- "http://mouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq2],rma::options[only$eq'id'],pipe::list[xstructures$eq'id'],service::differential_rows[set$eq'P56'][domain1$eq'1097'][domain1_threshold$eq'0,50'][domain2$eq'194'][domain2_threshold$eq'1,50'][start_row$eq0][num_rows$eq2000]"

lhaVsgrey <- extractISH_de(lhaDc)
lhaVsHy <- extractISH_de(lhaVsHypoDc)


lhaSpec <- merge(lhaVsgrey, lhaVsHy, by = 'object_id')
lhaSpecExt <- merge(lhaSpec, extracellGenesWithOrthologsGeneIds, by.x = 'entrez.x', by.y ='gene_id')
write.xlsx(lhaSpecExt, file = 'results/extracellularSpecific.xlsx', row.names = F, sheetName = "LHA", append = TRUE)



