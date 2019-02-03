# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

geneIDs <- function(){
  query<-'
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
  PREFIX dc: <http://purl.org/dc/elements/1.1/>
  PREFIX dcterms: <http://purl.org/dc/terms/>
  PREFIX dbpedia2: <http://dbpedia.org/property/>
  PREFIX dbpedia: <http://dbpedia.org/>
  PREFIX foaf: <http://xmlns.com/foaf/0.1/>
  PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

  SELECT DISTINCT ?geneIDURI ?geneIDLabel
  where {
    { ?geneIDURI ?p <http://edamontology.org/data_2295> .
      ?geneIDURI rdfs:label ?geneIDLabel .}
    UNION
    {
      <http://edamontology.org/data_1026> rdfs:label ?geneIDLabel .
      ?geneIDURI rdfs:label ?geneIDLabel .
    }
  }
  '
  endpoint <- 'https://www.ebi.ac.uk/rdf/services/sparql'

  SPARQL(endpoint,query)$results
}

jsonldContext <- function(){
  "http://schema.org/"
}

jsonldCreator <- function(name=NULL, email=NULL ){
  creator = list("@type"="Person")
  if (is.null(name)) {
    name = as.character(Sys.info()['user'])
    creator = c(creator, "name" = name)
  }
  
  if (!is.null(email)){
    creator = c(creator, "email"= email)
  }
  
  #rjson::toJSON(creator)
  creator
}

# type can be matrix or list
# in case of matrix describe what are columns and rows
jsonldData <- function(data, id, creator, description, datePublished, keywords, dataFormat, measurementTechnique, type="Dataset"){
  dataset = list("@type"=type,
                 "@id"=id,
                 "data"=data,
                 "measurementTechnique" = measurementTechnique,
                 "creator"= creator,
                 "fileformat"=dataFormat,
                 "description"=description,
                 "datePublished"=datePublished,
                 "keywords"=keywords)
  # rjson::toJSON(dataset)
  dataset
}

# row and column type must be selected by the geneID functions as geneIDURI
# usually rows represents genes
# dataPath can be an URL, path to file or a path to directory
# in this case matrix contains only one type of data, meta information are not available.
jsonldMatrix <- function(x=NULL, dataPath, rowType, rowLabel, columnType, columLabel,...){
  data = jsonldData()
  data = list("@type"="Matrix",
              "url"="dataPath",
              "rows" = list("@type"=rowType,
                            "name"=rowLabel),
              "columns"= list("@type"=columnType,
                              "name"=columnLabel))
  if (!is.null(x)){
    n=nrow(x)
    m=ncol(x)
    data$rows=c(data$rows, "size"=n)
    data$columns=c(data$columns, "size"=m)
  }
    # rjson::toJSON(data)
  data

}

# row and column type must be selected by the geneID functions as geneIDURI
# usually rows represents genes
# dataPath can be an URL, path to file or a path to directory
jsonldList <- function(dataPath, type, label){
  data = list("url"=dataPath,
              "list"=list("@type"=type,
                          "name"=label))
  # rjson::toJSON(data)
  data
}


jsonldHtmlScript <- function(data){
  script=htmltools::tags$script(jsonlite::prettify(rjson::toJSON(data)) )
  htmltools::tagAppendAttributes(script, "type"="application/ld+json")
}

# Dump the json-ld into a json file to be embedded into html
bhjsonldDump <- function(file){

}
labels=geneIDs()
columnLabel = "Ensembl gene ID"
columnType = labels[which(labels$geneIDLabel==columnLabel),"geneIDURI"]


d<-jsonldData(
  context=
  data=jsonldMatrix(x=dataset_mad_regout_k_nsclc_V2@raw.data,
                    dataPath="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                    rowType="http://edamontology.org/data_3273",
                    rowLabel = "10XCellBarcoding",
                    columnType = columnType,
                    columLabel = columnLabel),
  id="http://mock_schema.org/1.1",
  creator= jsonldCreator(),
  description="Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are 2,700 single cells that were sequenced on the Illumina NextSeq 500.",
  datePublished="2017-09-13",
  keywords="matrix, rnaseq, single cell, 10xgenomics",
  dataFormat = "application/x-gzip",
  measurementTechnique = "10x scRNA-seq")



jsonldHtmlScript(d)
