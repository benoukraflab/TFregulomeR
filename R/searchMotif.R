#' Search motif PFM and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#'
#' This function allows you to obtain motif PFM matrix and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#' @param id Required. TFregulomeR ID.
#' @param motif_format Motif PFM format, either in MEME by default or TRANSFAC.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return MethMotif class object
#' @keywords MethMotif
#' @export
#' @examples
#' K562_CEBPB <- searchMotif(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_transfac <- searchMotif(id = "MM1_HSA_K562_CEBPB",
#'                                    motif_format = "TRANSFAC")

searchMotif <- function(id, motif_format = "MEME",
                        server = "ca", TFregulome_url,
                        local_db_path = NULL) {
  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)
  results <- .searchMotif(
    id = id,
    motif_format = motif_format,
    api_object = api_object
  )
  return(results)
}

#' Search motif PFM and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#'
#' This function allows you to obtain motif PFM matrix and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#' @param id Required. TFregulomeR ID.
#' @param motif_format Motif PFM format, either in MEME by default or TRANSFAC.
#' @param api_object s4 object for interacting with the TFregulome database
#' @return MethMotif class object
#' @keywords MethMotif
#' @export
#' @examples
#' K562_CEBPB <- searchMotif(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_transfac <- searchMotif(id = "MM1_HSA_K562_CEBPB",
#'                                    motif_format = "TRANSFAC")

.searchMotif <- function(id, motif_format = "MEME", api_object) {
  # check motif_format MEME and TRANSFAC.
  motif_format <- toupper(motif_format)
  if (motif_format != "MEME" && motif_format != "TRANSFAC") {
    stop("Please check motif_format! Currently we only support MEME (default) and TRANSFAC formats!")
  }
  if (missing(id)) {
    stop("Please input TFregulome ID using 'id = '!")
  } else if (!is(api_object, "API") && !validObject(api_object)) {
    stop("Invalid API object!")
  } else {
    # make the request
    request_content_df <- apiRequest(api_object,  id = id)
  }
  # check db output
  if (is.null(request_content_df) || nrow(request_content_df) == 0) {
    return(NULL)
  } else {
    message("There are a matched record exported in a MethMotif object.")
    # using the previous backup TFregulome_url to parse into readMMmotif()
    if (motif_format == "MEME") {
      motif_file_path <- request_content_df[1, "motif_MEME"]
    } else {
      motif_file_path <- request_content_df[1, "motif_TRANSFAC"]
    }
    betaScore_file_path <- request_content_df[1, "beta_score_matrix"]
    num_peak <- request_content_df[1, "peak_with_motif_num"]
    nsites <- request_content_df[1, "TFBS_num"]

    # MEME file path for background extraction when motif formati is TRANSFAC
    motif_file_path_MEME <- request_content_df[1, "motif_MEME"]
    methmotif_item_motif <- readMMmotif(
      motif_file_path = motif_file_path,
      motif_format = motif_format,
      id = id,
      num_peak = as.integer(num_peak),
      nsites = as.integer(nsites),
      motif_file_path_MEME = motif_file_path_MEME
    )
    methmotif_item_betaScore <- readBetaScoreMatrix(
      betaScore_file_path = betaScore_file_path
    )
    methmotif_item <- new("MethMotif")
    methmotif_item <- updateMethMotif(
      methmotif_item,
      MMBetaScore = methmotif_item_betaScore,
      MMmotif = methmotif_item_motif
    )

    return(methmotif_item)
  }
}
