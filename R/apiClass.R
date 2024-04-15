######## S4 class API ########
setClass(
  # Set the name for the class
  "API",
  # define the slots
  slots = c(
    address = "character",
    hg38_gene_name = "character",
    hg38_new_gene_name = "character",
    mm10_gene_name = "character"
  ),
  # Set the default values for the slots.
  prototype = list(
    address = NA_character_,
    hg38_gene_name = NA_character_,
    hg38_new_gene_name = NA_character_,
    mm10_gene_name = NA_character_
  )
)

setGeneric("apiRequest", function(x, query_index, query_value, id) standardGeneric("apiRequest"))

######## S4 class LocalDatabase ########
setClass(
  # Set the name for the class
  "LocalDatabase",
  # inherits
  contains = "API"
)

setValidity("LocalDatabase", function(object) {
  if (is.null(object@address)) {
    "the path to the local database was not provided!"
  } else if (!endsWith(object@address, suffix = "/tfregulome.sqlite")) {
    "local SQLite database should be call 'tfregulome.sqlite'!"
  } else if (!file.exists(object@address)) {
    "cannot locate the SQLite database from the provided path!"
  } else {
    TRUE
  }
})

# helper
LocalDatabase <- function(address) {
  # set genomeAnnotate files
  hg38_gene_name <- gsub("/tfregulome.sqlite", "/TFregulomeR/genomeAnnotate/hg38_UCSC_to_GeneName.txt.gz", address)
  hg38_new_gene_name <- gsub("/tfregulome.sqlite", "/TFregulomeR/genomeAnnotate/hg38_UCSC_to_GeneName_NewVersion.txt.gz", address)
  mm10_gene_name <- gsub("/tfregulome.sqlite", "/TFregulomeR/genomeAnnotate/mm10_UCSC_to_GeneName.txt.gz", address)

  new("LocalDatabase", address = address, hg38_gene_name = hg38_gene_name,
    hg38_new_gene_name = hg38_new_gene_name, mm10_gene_name = mm10_gene_name)
}

setMethod("show", "LocalDatabase", function(object) {
  cat(is(object)[[1]], "\n",
      "  Address: ", object@address, "\n",
      sep = ""
  )
})

setMethod("apiRequest", "LocalDatabase", function(x, query_index, query_value, id) {
  # get the root path of the local tfregulome directory
  base_dir <- paste0(dirname(x@address), "/TFregulome_database")

  # prepare query
  if (!missing(id)) {
    query <- paste0("SELECT * FROM TFBS_table WHERE UPPER(id)=UPPER('", id, "')")
  } else if (!missing(query_index) && !missing(query_value)) {
    if (sum(query_index) == 0) {
      query <- "SELECT * FROM TFBS_table"
    } else {
      query <- paste0("SELECT * FROM TFBS_table WHERE ", paste0("UPPER(", sub("=", ")=UPPER('", query_value[query_index == 1]), "')", collapse = " AND "))
    }
  } else {
    stop("Trying to query database without parameters")
  }

  # query the database
  con <- RSQLite::dbConnect(RSQLite::SQLite(), x@address)
  results <- tryCatch({
    RSQLite::dbGetQuery(con, query)
  },
  error = function(cond) {
    message("There is a warning to connect TFregulomeR SQLite database!")
    message("Advice:")
    message("1) Check the path to the local database;")
    message("2) Check dependent package 'RSQLite';")
    message(paste0("warning: ", cond))
    return(NULL)
  })
  RSQLite::dbDisconnect(con)

  # modify the output table
  if (nrow(results) > 0) {
    results_dir <- paste0(base_dir, "_", results$species, "/", results$organ)
    results$motif_MEME <- paste0(
      results_dir, "/motif_matrix/", results$motif_MEME)
    results$motif_TRANSFAC <- paste0(
      results_dir, "/motif_matrix/", results$motif_TRANSFAC)
    results$beta_score_matrix <- paste0(
      results_dir, "/beta_score_matrix/", results$beta_score_matrix)
    results$all_peak_file <- paste0(
      results_dir, "/TF_all_peaks/", results$all_peak_file)
    results$peak_with_motif_file <- paste0(
      results_dir, "/TF_peaks_with_motif/", results$peak_with_motif_file)
    results$DNA_methylation_profile <- paste0(
      results_dir, "/DNA_methylation_profile/", results$DNA_methylation_profile)
    results$DNA_methylation_profile_200bp <- paste0(
      results_dir, "/DNA_methylation_profile_200bp/",
      results$DNA_methylation_profile_200bp)
    results$TFBS <- paste0(results_dir, "/TFBS/", results$TFBS)
    results$Ncor_between_MEME_ChIP_and_HOMER <- as.logical(
      results$Ncor_between_MEME_ChIP_and_HOMER)
  }
  return(results)
})

######## S4 class WebDatabase ########
setClass(
  # Set the name for the class
  "WebDatabase",
  # inherits
  contains = "API",
  # # define the slots
  slots = c(
    server = "character"
  ),
  # Set the default values for the slots.
  prototype = list(
    server = "ca"
  )
)

setValidity("WebDatabase", function(object) {
  # check server location
  if (object@server != "sg" && object@server != "ca") {
    "server should be either 'ca' (default) or 'sg'!"
  } else {
    TRUE
  }
})

# helper
WebDatabase <- function(address, server) {
  # make an appropriate API url
  if (missing(address)) {
    if (server == "sg") {
      address <- "https://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/"
    } else {
      address <- "https://methmotif.org/api/table_query/"
    }
  } else if (endsWith(address, suffix = "/index.php")) {
    address <- gsub("index.php", "", address)
    address <- paste0(address, "api/table_query/")
  } else if (endsWith(address, suffix = "/")) {
    address <- paste0(address, "api/table_query/")
  } else {
    address <- paste0(address, "/api/table_query/")
  }

  # set genomeAnnotate files
  gene_name_address <- gsub("api/table_query/", "api/TFregulomeR/genomeAnnotate/", address)
  hg38_gene_name <- paste0(gene_name_address, "hg38_UCSC_to_GeneName.txt")
  hg38_new_gene_name <- paste0(gene_name_address, "hg38_UCSC_to_GeneName_NewVersion.txt")
  mm10_gene_name <- paste0(gene_name_address, "mm10_UCSC_to_GeneName.txt")

  new("WebDatabase", address = address, server = server,
    hg38_gene_name = hg38_gene_name,
    hg38_new_gene_name = hg38_new_gene_name,
    mm10_gene_name = mm10_gene_name)
}

setMethod("show", "WebDatabase", function(object) {
  cat(is(object)[[1]], "\n",
      "  Address: ", object@address, "\n",
      "  Server: ", object@server, "\n",
      sep = ""
  )
})

setMethod("apiRequest", "WebDatabase", function(x, query_index, query_value, id) {
  # prepare query_url
  if (!missing(id)) {
    query_url <- paste0("listTFBS.php?AllTable=F&id=", id)
  } else if (!missing(query_index) && !missing(query_value)) {
    if (sum(query_index) == 0) {
      query_url <- "listTFBS.php?AllTable=T"
    } else {
      query_url <- paste0("listTFBS.php?AllTable=F&",
                          paste0(query_value[query_index == 1], collapse = "&"))
    }
  } else {
    stop("Trying to query database without parameters")
  }

  #parse JSON from API endpoint
  request_content_json <- tryCatch({
    jsonlite::fromJSON(paste0(x@address, query_url))
  },
  error = function(cond) {
    message("There is a warning to connect TFregulomeR API!")
    message("Advice:")
    message("1) Check internet access;")
    message("2) Check dependent package 'jsonlite';")
    message("3) Current TFregulomeR server is implemented in MethMotif database, whose homepage is 'https://bioinfo-csi.nus.edu.sg/methmotif/' or 'https://methmotif.org'. If MethMotif homepage url is no more valid, please Google 'MethMotif', and input the valid MethMotif homepage url using 'TFregulome_url = '.")
    message(paste0("warning: ", cond))
    return(NULL)
  })

  # return the data.frame
  if (is.null(request_content_json)) {
    message("Empty output for the request!")
    return(NULL)
  } else {
    request_content_df <- as.data.frame(request_content_json$TFBS_records)
    # final output and return data.frame
    if (nrow(request_content_df) == 0) {
      if (exists("request_content_json")) {
        message(request_content_json$message)
        return(NULL)
      } else {
        message("No matched records found.")
        return(NULL)
      }
    }
    return(request_content_df)
  }
})
