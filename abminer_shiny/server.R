library(shiny)
library(DT)
load('data/pvhd_sec_memb.RData')
load('data/vmh_sec_memb.RData')
load('data/arc_sec_memb.RData')
load('data/lha_sec_memb.RData')


# pvhd_sec_mem <- pvhd_sec_memb[1:20,]
# vmh_sec_memb <- vmh_sec_memb[1:20,]

pvhd_sec_memb <- cbind(pvhd_sec_memb, brain_explorer = pvhd_sec_memb$object_id)
pvhd_sec_memb$brain_explorer <- sapply(pvhd_sec_memb$brain_explorer, function (x){
  toString(tags$a(href=paste0("aibe://mouse.brain-map.org/grid_data/v1/visualize/", x,"?atlas=310"), x))
})

vmh_sec_memb <- cbind(vmh_sec_memb, brain_explorer = vmh_sec_memb$object_id)
vmh_sec_memb$brain_explorer <- sapply(vmh_sec_memb$brain_explorer, function (x){
  toString(tags$a(href=paste0("aibe://mouse.brain-map.org/grid_data/v1/visualize/", x,"?atlas=310"), x))
})

arc_sec_memb <- cbind(arc_sec_memb, brain_explorer = arc_sec_memb$object_id)
arc_sec_memb$brain_explorer <- sapply(arc_sec_memb$brain_explorer, function (x){
  toString(tags$a(href=paste0("aibe://mouse.brain-map.org/grid_data/v1/visualize/", x,"?atlas=310"), x))
})

lha_sec_memb <- cbind(lha_sec_memb, brain_explorer = lha_sec_memb$object_id)
lha_sec_memb$brain_explorer <- sapply(lha_sec_memb$brain_explorer, function (x){
  toString(tags$a(href=paste0("aibe://mouse.brain-map.org/grid_data/v1/visualize/", x,"?atlas=310"), x))
})



pvhd_sec_memb$object_id <- sapply(pvhd_sec_memb$object_id, function (x){
  toString(tags$a(href=paste0("http://mouse.brain-map.org/experiment/show/", x), x))
})

vmh_sec_memb$object_id <- sapply(vmh_sec_memb$object_id, function (x){
  toString(tags$a(href=paste0("http://mouse.brain-map.org/experiment/show/", x), x))
})

lha_sec_memb$object_id <- sapply(lha_sec_memb$object_id, function (x){
  toString(tags$a(href=paste0("http://mouse.brain-map.org/experiment/show/", x), x))
})

arc_sec_memb$object_id <- sapply(arc_sec_memb$object_id, function (x){
  toString(tags$a(href=paste0("http://mouse.brain-map.org/experiment/show/", x), x))
})


pvhd_sec_memb$ensembl_id <- sapply(pvhd_sec_memb$ensembl_id, function (x){
  toString(tags$a(href=paste0("http://bioweb/bioatlas/genes/", x), x))
})

vmh_sec_memb$ensembl_id <- sapply(vmh_sec_memb$ensembl_id, function (x){
  toString(tags$a(href=paste0("http://bioweb/bioatlas/genes/", x), x))
})

lha_sec_memb$ensembl_id <- sapply(lha_sec_memb$ensembl_id, function (x){
  toString(tags$a(href=paste0("http://bioweb/bioatlas/genes/", x), x))
})

arc_sec_memb$ensembl_id <- sapply(arc_sec_memb$ensembl_id, function (x){
  toString(tags$a(href=paste0("http://bioweb/bioatlas/genes/", x), x))
})


pvhd_sec_memb[,'fold_change'] <- as.numeric(pvhd_sec_memb[,'fold_change'])
vmh_sec_memb[,'fold_change'] <- as.numeric(vmh_sec_memb[,'fold_change'])
arc_sec_memb[,'fold_change'] <- as.numeric(arc_sec_memb[,'fold_change'])
lha_sec_memb[,'fold_change'] <- as.numeric(lha_sec_memb[,'fold_change'])


# Define server logic required to summarize and view the selected
# dataset
function(input, output) {

  # By declaring datasetInput as a reactive expression we ensure
  # that:
  #
  #  1) It is only called when the inputs it depends on changes
  #  2) The computation and result are shared by all the callers
  #	  (it only executes a single time)
  #
  datasetInput <- reactive({
    switch(input$structure,
           "VMH" = vmh_sec_memb,
           "PVHD" = pvhd_sec_memb,
           "ARC" = arc_sec_memb,
           "LHA" = lha_sec_memb)
  })

  # The output$view depends on both the databaseInput reactive
  # expression and input$obs, so will be re-executed whenever
  # input$dataset or input$obs is changed.
  # output$view <- renderTable({
  #   head(datasetInput(), n = input$obs)
  # })

  # output$view <- DT::renderDataTable(
  #   DT::datatable(datasetInput())
  # )

  output$view <- DT::renderDataTable(
    view <- DT::datatable(datasetInput(), escape = FALSE, options = list(order = list(list(4, 'desc'))))
  )
}
