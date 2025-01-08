# Get list of species that are represented in the Protein Ontology

source("./R/endpoints.R")
#source("endpoints.R")

library(magrittr)

# SPARQL to get a list of species
SPARQL_query <-
  "
  PREFIX pro: <http://purl.org/hpi/patchr#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
  PREFIX obo:   <http://purl.obolibrary.org/obo/>

  SELECT DISTINCT ?species_ID ?species_name ?broad_synonym
                  ?exact_synonym ?related_synonym_2

  # Search Protein Ontology
  FROM <http://purl.obolibrary.org/obo/merged/PR>

  WHERE {
    ?PRO_term rdfs:label ?PRO_name .

    ?PRO_term rdfs:subClassOf ?l .

    ?l owl:someValuesFrom ?species_ID .

  OPTIONAL {
    ?species_ID rdfs:label ?species_name
  }
  OPTIONAL {
    ?species_ID oboInOwl:hasBroadSynonym ?broad_synonym
  }
  OPTIONAL {
    ?species_ID oboInOwl:hasExactSynonym ?exact_synonym
  }
  OPTIONAL {
    ?species_ID oboInOwl:hasRelatedSynonym ?related_synonym .
    BIND(REPLACE(str(?related_synonym), '[. -]','') AS ?related_synonym_2) .
    BIND(strlen(str(?related_synonym_2)) AS ?num_characters)
      }

  FILTER (strstarts(str(?species_ID), 'http://purl.obolibrary.org/obo/NCBITaxon_'))
  FILTER (?num_characters > 2 || !bound(?num_characters))
  }
  order by strlen(str(?related_synonym_2))
"

json_header <- c('Accept'="application/json")

# Run SPARQL on the endpoint
result <- httr::POST(onto_endpoint,
                     body = list(query = SPARQL_query),
                     httr::user_agent(R.version.string),
                     httr::add_headers(json_header))

# Will show a warning/error if there is any
httr::stop_for_status(result)

# Get result in text JSON
x <- httr::content(result, "text", encoding = "UTF-8")

# Convert from JSON to a list
species_in_PRO <- jsonlite::fromJSON(x, flatten = TRUE)

# Extract the dataframe
species_in_PRO <- species_in_PRO$results$bindings

# Remove unneeded info
species_in_PRO <- species_in_PRO %>% dplyr::select(-ends_with(c(".type", ".datatype", "lang")))

# Remove this part that was added onto the column names
colnames(species_in_PRO) <- gsub(".value", "", colnames(species_in_PRO))

# Format name for drop down menu
species_in_PRO$full_name <-
  ifelse(
    is.na(species_in_PRO$exact_synonym) == FALSE &
      is.na(species_in_PRO$related_synonym_2) == FALSE,
    paste0(
      species_in_PRO$species_name,
      " (",
      species_in_PRO$exact_synonym,
      "/",
      species_in_PRO$related_synonym_2,
      ")"
    ),
    ifelse(
      is.na(species_in_PRO$exact_synonym) == TRUE &
        is.na(species_in_PRO$related_synonym_2) == FALSE,
      paste0(
        species_in_PRO$species_name,
        " (",
        species_in_PRO$related_synonym_2,
        ")"
      ),
      ifelse(
        is.na(species_in_PRO$exact_synonym) == FALSE &
          is.na(species_in_PRO$related_synonym_2) == TRUE,
        paste0(
          species_in_PRO$species_name,
          " (",
          species_in_PRO$exact_synonym,
          ")"
        ),
        species_in_PRO$species_name
      )
    )
  )

# Add row for species not specified
species_in_PRO[nrow(species_in_PRO) + 1, ] <-
  c(
    NA,
    "NA/No species specified",
    "NA/No species specified",
    "NA/No species specified",
    "NA/No species specified",
    "NA/No species specified"
  )

# Sort most species by alpbetical order
species_in_PRO <-
  species_in_PRO[order(species_in_PRO$species_name), ]

# Change order so the most commonly used options are at the top
new_order <-
  c(
    intersect(
      c("Homo sapiens", "Mus musculus", "NA/No species specified"),
      species_in_PRO$species_name
    ),
    setdiff(
      species_in_PRO$species_name,
      c("Homo sapiens", "Mus musculus", "NA/No species specified")
    )
  )

# Order by species name
species_in_PRO$species_name <- factor(species_in_PRO$species_name, levels = c(new_order))
species_in_PRO <- species_in_PRO[order(species_in_PRO$species_name), ]

# Just get NCBI taxon ID without URI
species_in_PRO$NCBI_taxon <- species_in_PRO$species_ID
species_in_PRO$NCBI_taxon <- gsub(".*/obo/", "", species_in_PRO$NCBI_taxon)
species_in_PRO$NCBI_taxon <- gsub('.{1}$', '', species_in_PRO$NCBI_taxon)
