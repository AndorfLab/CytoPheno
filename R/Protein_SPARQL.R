#' Function for SPARQL query and processing of the resulting data
#'
#' @param unique_marker_names A vector of marker names
#' @param NCBI_taxon_ID String specifying the NCBI taxon ID
#' @param onto_endpoint String specifying the Ontobee endpoint
#' @param Match_step String indicating step in the pipeline
#'
#' @importFrom magrittr %>%
#' @return Dataframe of query results
#' @export
#'
#' @examples
#' unique_syn_marker_names = c("CD4", "CD16")
#' NCBI_taxon_ID_short = "NCBITaxon_9606"
#' onto_endpoint <- "https://sparql.hegroup.org/sparql"
#' Match_step <- "Step 2"
#' PRO_SPARQL(unique_marker_names = unique_syn_marker_names, NCBI_taxon_ID = NCBI_taxon_ID_short,
#' onto_endpoint = onto_endpoint, Match_step = Match_step)

PRO_SPARQL <- function(unique_marker_names, NCBI_taxon_ID, onto_endpoint, Match_step) {

  all_ids <- list()

  for (marker in unique_marker_names) {

    syn_cols <- c("has_exact_synonym", "has_broad_synonym", "has_related_synonym", "has_narrow_synonym", "Species")

    specific_marker <-  shQuote(marker)

    SPARQL_query<-
      paste0(
        "
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
  PREFIX obo: <http://purl.obolibrary.org/obo/>

  SELECT DISTINCT ?PRO_term ?PRO_name
  ?has_exact_synonym ?has_broad_synonym ?has_related_synonym ?has_narrow_synonym ?Species

  FROM <http://purl.obolibrary.org/obo/merged/PR>

  WHERE {
    ?PRO_term rdfs:label ?PRO_name .

    OPTIONAL {
      ?PRO_term rdfs:subClassOf ?l .
      ?l owl:someValuesFrom ?Species
    }

      FILTER (?Species = obo:",
        NCBI_taxon_ID,
        "|| !bound(?Species) || (strstarts(str(?Species), 'http://purl.obolibrary.org/obo/PR_'))) .

    OPTIONAL {
      ?PRO_term oboInOwl:hasExactSynonym ?has_exact_synonym
    }
    OPTIONAL {
      ?PRO_term oboInOwl:hasBroadSynonym ?has_broad_synonym .
    }
    OPTIONAL {
      ?PRO_term oboInOwl:hasNarrowSynonym ?has_narrow_synonym
    }
    OPTIONAL {
      ?PRO_term oboInOwl:hasRelatedSynonym ?has_related_synonym
    }

    {  SELECT DISTINCT ?PRO_term

      WHERE {
        ?PRO_term rdfs:label ?_PRO_name .
        OPTIONAL {
          ?PRO_term oboInOwl:hasExactSynonym ?has_exact_synonym .
        }
        OPTIONAL {
          ?PRO_term oboInOwl:hasBroadSynonym ?has_broad_synonym .
        }
        OPTIONAL {
          ?PRO_term oboInOwl:hasNarrowSynonym ?has_narrow_synonym .
        }
        OPTIONAL {
          ?PRO_term oboInOwl:hasRelatedSynonym ?has_related_synonym .
        }

         FILTER NOT EXISTS {?PRO_term owl:deprecated TRUE}.

         FILTER (UCASE(REPLACE(str(?_PRO_name),'[ -.]','')) =",
        specific_marker,
        "|| UCASE(REPLACE(str(?has_exact_synonym),'[ -.]','')) =",
        specific_marker,
        "|| UCASE(REPLACE(str(?has_broad_synonym),'[ -.]','')) =",
        specific_marker,
        "|| UCASE(REPLACE(str(?has_narrow_synonym),'[ -.]','')) =",
        specific_marker,
        "|| UCASE(REPLACE(str(?has_related_synonym),'[ -.]','')) =",
        specific_marker, ") .
        }
      }
    }
  "
      )

    # Run SPARQL on the endpoint
    result <- httr::POST(onto_endpoint,
                         body = list(query = SPARQL_query),
                         httr::user_agent(R.version.string))

    # Will show a warning/error if there is any
    httr::stop_for_status(result)

    # Get result in text JSON
    x <- httr::content(result, "text", encoding = "UTF-8")

    # Convert from JSON to a list
    df <- jsonlite::fromJSON(x, flatten = TRUE)

    # Extract the dataframe
    df <- df$results$bindings

    # Put the results in a list of data frames
    if (length(df) != 0) {

      # Remove unneeded info
      df <- df %>% dplyr::select(-dplyr::ends_with(c(".type", ".datatype", "lang")))

      # Remove this part that was added onto the column names
      colnames(df) <- gsub(".value", "", colnames(df))

      # Add columns even if that marker didn't have that type of synoymn listed
      df[syn_cols[!(syn_cols %in% colnames(df))]] <- NA

      # Make everything upper case, remove spaces, dashes
      df_match_type <- dplyr::mutate_all(df, .funs = toupper)
      df_match_type <- as.data.frame(lapply(df_match_type, function(y) gsub('-', '', y)))
      df_match_type <- as.data.frame(lapply(df_match_type, function(y) gsub('_', '', y)))
      df_match_type <- as.data.frame(lapply(df_match_type, function(y) gsub(' ', '', y)))
      df_match_type <- as.data.frame(lapply(df_match_type, function(y) gsub('\\.', '', y)))

      # Remove parenthesis and everything in between ()
      df_match_type[2:ncol(df_match_type)] <- lapply(df_match_type[2:ncol(df_match_type)], function(x) as.character(gsub("\\s*\\([^\\)]+\\)", "", x)))

      unique_PRO_terms <- sort(unique(df_match_type$PRO_term))

      match_types <- c()

      # Indicate what type of match it is
      for (PRO_term in unique_PRO_terms) {
        df_match_type_subsetted <- df_match_type[df_match_type$PRO_term == PRO_term,]
        if (sum(grepl(paste0('^', marker, '$'), df_match_type_subsetted[,"PRO_term"]) == TRUE) +
            sum(grepl(paste0('^', marker, '$'), df_match_type_subsetted[,"PRO_name"]) == TRUE) +
            sum(grepl(paste0('^', marker, '$'), df_match_type_subsetted[,"has_exact_synonym"]) == TRUE) > 0) {
          match_type_exact <- "Primary"
          match_types <- c(match_types, match_type_exact)
        } else if (sum(grepl(paste0('^', marker, '$'), df_match_type_subsetted[,"has_related_synonym"]) == TRUE) +
                   sum(grepl(paste0('^', marker, '$'), df_match_type_subsetted[,"has_narrow_synonym"]) == TRUE) +
                   sum(grepl(paste0('^', marker, '$'), df_match_type_subsetted[,"has_broad_synonym"]) == TRUE) > 0) {
          match_type_secondary <- "Secondary"
          match_types <- c(match_types, match_type_secondary)
        }
      }

      PRO_term_to_match_type <- data.frame(PRO_term = sort(unique(df$PRO_term)), Match_type = match_types)

      # Just get the unique PRO IDs
      matched_IDs <- dplyr::distinct(df, PRO_term, PRO_name, Species)

      final_df <- data.frame(matched_IDs, Match_step)

      final_df <- merge(final_df, PRO_term_to_match_type)

      final_df[is.na(final_df)] <- "Not species specific"

      unique_species_names <- unique(final_df$Species)

      remove_secondary_per_species <- data.frame()

      for(species in unique_species_names) {
        df_species_subsetted <- final_df[final_df$Species == species,]
        if(any(df_species_subsetted$Match_type == "Primary")) {
          df_per_species <- df_species_subsetted[df_species_subsetted$Match_type != "Secondary",]
          remove_secondary_per_species <- rbind(remove_secondary_per_species, df_per_species)
        } else {
          remove_secondary_per_species <- rbind(remove_secondary_per_species, df_species_subsetted)
        }
      }

      remove_secondary_per_species <- remove_secondary_per_species[ , c("PRO_term", "PRO_name", "Species", "Match_type", "Match_step")]

      all_ids[[marker]] <- remove_secondary_per_species

      # If a protein was not found in PRO, indicate that
    } else if (length(df) == 0) {
      all_ids[[marker]] <-
        data.frame(PRO_term = "", PRO_name = "", Species = "", Match_type = "", Match_step = "")
    }
  }
  return(all_ids)
}
