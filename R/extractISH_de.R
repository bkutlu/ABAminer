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
