#' Browse the current data available in TFregulomeR
#'
#' This function allows you to get the current data in TFregulomeR
#' @param species The species of interest
#' @param organ The organ of interest
#' @param sample_type The sample type of interest
#' @param cell_tissue_name The name of tissue or cell of interest
#' @param tf The TF of interest
#' @param disease_state The disease state of interest
#' @param source The source of interest
#' @param server server location to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return  data.frame containing the information of the queried TFBSs in TFregulomeR
#' @keywords ChIPseq data
#' @export
#' @examples
#' TFBS_brain <- dataBrowser(organ = "brain")

dataBrowser <- function(
  species, organ, sample_type, cell_tissue_name, tf, disease_state, source,
  server = "ca", TFregulome_url, local_db_path = NULL
) {
  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)

  query_index <- rep(0, 7)
  query_value <- rep("", 7)
  query_terms <- c(
    "species", "organ", "sample_type", "cell_or_tissue_name",
    "tf", "disease_state", "source"
  )
  names(query_index) <- query_terms
  names(query_value) <- query_terms
  if (!missing(species)) {
    query_index["species"] <- 1
    query_value["species"] <- paste0("species=", species)
  }
  if (!missing(organ)) {
    query_index["organ"] <- 1
    query_value["organ"] <- paste0("organ=", organ)
  }
  if (!missing(sample_type)) {
    query_index["sample_type"] <- 1
    query_value["sample_type"] <- paste0("sample_type=", sample_type)
  }
  if (!missing(cell_tissue_name)) {
    query_index["cell_tissue_name"] <- 1
    query_value["cell_tissue_name"] <- paste0(
      "cell_tissue_name=", cell_tissue_name
    )
  }
  if (!missing(tf)) {
    query_index["tf"] <- 1
    query_value["tf"] <- paste0("TF=", tf)
  }
  if (!missing(disease_state)) {
    query_index["disease_state"] <- 1
    query_value["disease_state"] <- paste0("disease_state=", disease_state)
  }
  if (!missing(source)) {
    query_index["source"] <- 1
    query_value["source"] <- paste0("source=", source)
  }

  # make the request
  request_content_df <- apiRequest(api_object, query_index, query_value)
  # check db output
  if (is.null(request_content_df)) {
    return(NULL)
  } else {
    request_content_df_output <- request_content_df[, c(
      "ID", "species", "organ",
      "sample_type",
      "cell_tissue_name",
      "description",
      "disease_state",
      "TF", "source",
      "source_ID",
      "peak_num",
      "peak_with_motif_num",
      "Consistent_with_HOCOMOCO_JASPAR",
      "Ncor_between_MEME_ChIP_and_HOMER"
    )]
    output_message(request_content_df_output)
    return(request_content_df_output)
  }
}


#' A structured format for printing to the console
#'
#' This function facilitates the printing of select metrics calculated from
#' the the queried TFBSs in TFregulomeR. This is an internal helper function
#' not meant to be called on its own.
#' @param requested_df data.frame containing the information of the queried TFBSs in TFregulomeR
#' @return  prints to console select metrics from the queried TFBSs in TFregulomeR
#' @keywords internal
#' @export
#' @examples
#' output_message(request_content_df_output)

output_message <- function(requested_df) {
  tfbs_num <- nrow(requested_df)
  species_result <- unique(requested_df$species)
  species_num <- length(species_result)
  organ_result <- unique(requested_df$organ)
  organ_num <- length(organ_result)
  sample_type_result <- unique(requested_df$sample_type)
  sample_type_num <- length(sample_type_result)
  cell_tissue_name_num <- length(unique(requested_df$cell_tissue_name))
  disease_state_result <- unique(requested_df$disease_state)
  disease_state_num <- length(disease_state_result)
  tf_num <- length(unique(requested_df$TF))
  source_result <- unique(requested_df$source)
  message(paste0(tfbs_num, " record(s) found: ..."))
  message(paste0("... covering ", tf_num, " TF(s)"))
  message(paste0("... from ", species_num, " species:"))
  message(paste0("... ...", paste0(species_result, collapse = ", ")))
  message(paste0("... from ", organ_num, " organ(s):"))
  message(paste0("... ... ", paste0(organ_result, collapse = ", ")))
  message(paste0("... in ", sample_type_num, " sample type(s):"))
  message(paste0("... ... ", paste0(sample_type_result, collapse = ", ")))
  message(paste(
    "... in ", cell_tissue_name_num, " different cell(s) or tissue(s)"
  ))
  message(paste0("... in ", disease_state_num, " type(s) of disease state(s):"))
  message(paste0("... ... ", paste0(disease_state_result, collapse = ", ")))
  message(paste0(
    "... from the source(s): ", paste0(source_result, collapse = ", ")
  ))
}
