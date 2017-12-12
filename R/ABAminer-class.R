.ABAminer <- setClass("ABAminer",
                      slots=c(df="DataFrame"))

setGeneric("differential", signature="object",
           function(object,
                    atlas_experiment_id=NA,
                    atlas_gene_id=NA,
                    entrez_id=NA,
                    gene_symbol=NA,
                    section=NA,
                    fold_change=NA,
                    target_sum=NA,
                    contrast_sum=TRUE)
             standardGeneric ("differential"))




