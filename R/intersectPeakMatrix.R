#' intersectPeakMatrix
#'
#' This function allows you to obtain the pair-wise intersected regions, along with the DNA methylation profiles, between two lists of peak sets, as well as (Meth)Motif logos
#' @param peak_id_x Character of vector, each of which is a unique TFregulomeR ID.
#' @param motif_only_for_id_x Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulomeR ID in peak_id_x.
#' @param user_peak_list_x A list of data.frames, each of which contains user's own bed-format peak regions for peak list x.
#' @param user_peak_x_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_peak_list_x. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulomeR, its TFregulomeR ID should be used here correspondingly.
#' @param peak_id_y Character of vector, each of which is a unique TFregulomeR ID.
#' @param motif_only_for_id_y Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulomeR ID in peak_id_y.
#' @param user_peak_list_y A list of data.frames, each of which contains user's own bed-format peak regions for peak list y.
#' @param user_peak_y_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_peak_list_y. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulomeR, its TFregulomeR ID should be used here correspondingly.
#' @param methylation_profile_in_narrow_region Either TRUE (default) of FALSE. If TRUE, methylation states in 200bp window surrounding peak summits for each intersected peak pair from peak_id_x (peak_id_y) and user_peak_list_x (user_peak_list_y) with TFregulomeR ID.
#' @param external_source a bed-like data.frame files with the fourth column as the score to be profiled in pairwise comparison regions.
#' @param motif_type Motif PFM format, either in MEME by default or TRANSFAC.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @param local_db_path The complete path to the SQLite implementation of TFregulomeR database available at "https://methmotif.org/API_ZIPPED.zip"
#' @return  matrix of IntersectPeakMatrix class objects
#' @keywords intersectPeakMatrix
#' @export
#' @examples
#' peak_id_x <- c("MM1_HSA_K562_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' peak_id_y <- c("MM1_HSA_HepG2_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' intersectPeakMatrix_output <- intersectPeakMatrix(peak_id_x=peak_id_x,
#'                                                   motif_only_for_id_x=TRUE,
#'                                                   peak_id_y=peak_id_y,
#'                                                   motif_only_for_id_y=TRUE)


intersectPeakMatrix <- function(
  peak_id_x,
  motif_only_for_id_x = FALSE,
  user_peak_list_x,
  user_peak_x_id,
  peak_id_y,
  motif_only_for_id_y = FALSE,
  user_peak_list_y,
  user_peak_y_id,
  methylation_profile_in_narrow_region = FALSE,
  external_source,
  motif_type = "MEME",
  server = "ca",
  TFregulome_url,
  local_db_path = NULL
) {
  # check the input argument
  if (missing(peak_id_x) && missing(user_peak_list_x)) {
    stop("No peak list x input. Please input TFregulomeR peaks using TFregulomeR ID(s) by 'peak_id_x = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_peak_list_x = '")
  }
  if (missing(peak_id_y) && missing(user_peak_list_y)) {
    stop("No peak list y input. Please input TFregulomeR peaks using TFregulomeR ID(s) by 'peak_id_y = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_peak_list_y = '")
  }
  if ((!missing(user_peak_list_x) && !is.list(user_peak_list_x)) ||
        (!missing(user_peak_list_y) && !is.list(user_peak_list_y))) {
    stop("The class of input 'user_peak_list_x' and 'user_peak_list_y' should be 'list', a list of bed-like data.frame storing peak regions!")
  }
  if (!is.logical(motif_only_for_id_x) || !is.logical(motif_only_for_id_y)) {
    stop("motif_only_for_id_x and motif_only_for_id_y should be either TRUE or FALSE (default)")
  }
  if (!is.logical(methylation_profile_in_narrow_region)) {
    stop("methylation_profile_in_narrow_region should be either TRUE or FALSE (default)")
  }
  if (!missing(external_source) && !is.data.frame(external_source)) {
    stop("external_source should be data.frame")
  }
  if (motif_type != "MEME" && motif_type != "TRANSFAC") {
    stop("motif_type should be either 'MEME' (default) or 'TRANSFAC'!")
  }

  # if the external source is provided
  external_source_provided <- FALSE
  if (!missing(external_source) && nrow(external_source) > 0) {
    external_source_provided <- TRUE
    external_source_signal <- external_source[, seq(1, 4, 1)]
    colnames(external_source_signal) <- c("chr", "start", "end", "score")
    external_source_signal$id <- paste0(
      "external_source_", rownames(external_source_signal)
    )
    external_source_grange <- GenomicRanges::GRanges(
      external_source_signal$chr,
      IRanges::IRanges(
        external_source_signal$start + 1,
        external_source_signal$end
      ),
      id = external_source_signal$id
    )
  }

  # build api_object
  api_object <- .construct_api(server, TFregulome_url, local_db_path)

  message("TFregulomeR::intersectPeakMatrix() starting ... ...")
  if (methylation_profile_in_narrow_region) {
    message("You chose to profile the methylation levels in 200bp window around peak summits, if there is any peak loaded from TFregulomeR. It will make the program slow. Disable it if you want a speedy analysis and do not care about methylation")
  } else {
    message("You chose NOT to profile the methylation levels in 200bp window around peak summits")
  }

  # loading peak list x
  peak_results_x <- load_peak_list(
    peak_id_x, "x", motif_only_for_id_x, user_peak_list_x,
    user_peak_x_id, motif_type, api_object
  )
  peak_id_x_all <- peak_results_x$peak_id_all
  peak_list_x_all <- peak_results_x$peak_list_all
  is_x_TFregulome <- peak_results_x$is_TFregulome
  # loading peak list y
  peak_results_y <- load_peak_list(
    peak_id_y, "y", motif_only_for_id_y, user_peak_list_y,
    user_peak_y_id, motif_type, api_object
  )
  peak_id_y_all <- peak_results_y$peak_id_all
  peak_list_y_all <- peak_results_y$peak_list_all
  is_y_TFregulome <- peak_results_y$is_TFregulome

  # check if peak list x is empty
  if (length(peak_list_x_all) == 0) {
    message("No input in the peak list x. Function ends here!")
    return(NULL)
  }
  # check if compared peak list is empty
  if (length(peak_list_y_all) == 0) {
    message("No input in the peak list y. Function ends here!")
    return(NULL)
  }

  # start analysing
  intersection_matrix <- list()
  for (i in seq(1, length(peak_list_x_all), 1)) {
    id_x <- peak_id_x_all[i]
    message(paste0("Start analysing list x:", id_x, "... ..."))
    peak_x <- peak_list_x_all[[i]]
    isTFregulome_x <- is_x_TFregulome[i]
    # fetch peak x info
    x_peak_info <- fetch_peak_info(id_x, peak_x, isTFregulome_x, api_object)
    x_elements <- list(id = id_x, peak = peak_x, isTFregulome = isTFregulome_x)
    x_info <- c(x_elements, x_peak_info)

    # start comparing with peak set y
    for (j in seq(1, length(peak_list_y_all), 1)) {
      y_info <- list()
      id_y <- peak_id_y_all[j]
      message(paste0("... ... Start analysing list y:", id_y))
      peak_y <- peak_list_y_all[[j]]
      isTFregulome_y <- is_y_TFregulome[j]
      # fetch peak y info
      y_peak_info <- fetch_peak_info(id_y, peak_y, isTFregulome_y, api_object)
      y_elements <- list(id = id_y, peak = peak_y,
                         isTFregulome = isTFregulome_y)
      y_info <- c(y_elements, y_peak_info)

      # get peak x which intersects with y
      x_intersect_info <- intersect_peak_regions(
        x_info, y_info, external_source_provided, external_source_grange,
        external_source_signal, motif_type, methylation_profile_in_narrow_region
      )
      x_of_y_info <- c(x_info, x_intersect_info)

      # get peak y which intersects with x
      y_intersect_info <- intersect_peak_regions(
        y_info, x_info, external_source_provided, external_source_grange,
        external_source_signal, motif_type, methylation_profile_in_narrow_region
      )
      y_of_x_info <- c(y_info, y_intersect_info)

      #form an IntersectPeakMatrix object
      intersect_id <- paste0(id_x, "_[AND]_", id_y)
      new_IntersectPeakMatrix <- new("IntersectPeakMatrix")
      new_IntersectPeakMatrix <- updateIntersectPeakMatrix(
        theObject = new_IntersectPeakMatrix,
        id = intersect_id,
        id_x = id_x,
        overlap_percentage_x = x_of_y_info$intersect_percentage,
        isxTFregulomeID = isTFregulome_x,
        MethMotif_x = x_of_y_info$MethMotif,
        methylation_profile_x = x_of_y_info$meth_score_distri_target,
        external_signal_x = x_of_y_info$external_signal_in,
        tag_density_x = x_of_y_info$tag_density,
        id_y = id_y,
        overlap_percentage_y = y_of_x_info$intersect_percentage,
        isyTFregulomeID = isTFregulome_y,
        MethMotif_y = y_of_x_info$MethMotif,
        methylation_profile_y = y_of_x_info$meth_score_distri_target,
        external_signal_y = y_of_x_info$external_signal_in,
        tag_density_y = y_of_x_info$tag_density
      )
      intersection_matrix[[intersect_id]] <- new_IntersectPeakMatrix
    }
  }
  intersection_matrix_matrix <- matrix(
    intersection_matrix,
    nrow = length(peak_id_x_all),
    ncol = length(peak_id_y_all),
    byrow = TRUE
  )
  rownames(intersection_matrix_matrix) <- c(peak_id_x_all)
  colnames(intersection_matrix_matrix) <- c(peak_id_y_all)
  return(intersection_matrix_matrix)
}


load_peak_list <- function(
  peak_id, direction, motif_only, user_peak_list,
  user_peak_id, motif_type, api_object
) {
  message(paste0("Loading peak list ", direction, " ... ..."))
  peak_list_all <- list()
  # loading from TFregulomeR server
  TFregulome_peak_id <- c()
  is_TFregulome <- c()
  peak_list_count <- 0
  if (!missing(peak_id) && length(peak_id) > 0) {
    message(paste0("... You have ", length(peak_id)," TFBS(s) requested to be loaded from TFregulomeR server"))
    if (motif_only == TRUE) {
      message(paste0("... You chose to load TF peaks with motif only. Using 'motif_only_for_id_", direction, "' tunes your options"))
    } else {
      message(paste0("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_id_", direction, "' tunes your options"))
    }
    message("... loading TFBS(s) from TFregulomeR now")
    for (i in peak_id) {
      peak_i <- suppressMessages(.loadPeaks(
        id = i,
        includeMotifOnly = motif_only,
        api_object = api_object
      ))
      if (is.null(peak_i)) {
        message(paste0("... ... NO peak file for your id '", i, "'."))
      } else {
        peak_list_count <- peak_list_count + 1
        peak_list_all[[peak_list_count]] <- peak_i
        is_TFregulome <- c(is_TFregulome, TRUE)
        TFregulome_peak_id <- c(TFregulome_peak_id, i)
        message(paste0("... ... peak file loaded successfully for id '", i, "'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulomeR")
  }
  # users' peaks
  if (!missing(user_peak_list) && length(user_peak_list) > 0) {
    message(paste0("... You have ",length(user_peak_list)," customised peak set(s)"))
    if (missing(user_peak_id) || length(user_peak_id) != length(user_peak_list) || length(unique(user_peak_id)) != length(user_peak_list)) {
      message(paste0("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_peak_", direction, "1', 'user_peak_", direction, "2'..."))
      user_peak_id <- paste0("user_peak_", direction, seq(1,length(user_peak_list), 1))
    }
    # new user peak x id in case that any input peak set is empty
    user_peak_id_new <- c()
    for (i in seq(1, length(user_peak_list), 1)) {
      user_peak_i <- user_peak_list[[i]]
      if (nrow(user_peak_i) == 0) {
        message(paste0("... ... Your input peak set '", user_peak_id[i], "' is empty, so SKIP!"))
      } else {
        user_peak_id_new <- c(user_peak_id_new, user_peak_id[i])
        colname_new <- colnames(user_peak_i)
        colname_new[1] <- "chr"
        colname_new[2] <- "start"
        colname_new[3] <- "end"
        # check if the input peak set has fourth column for id
        no_id <- TRUE
        if (length(colname_new) >= 4) {
          if (length(unique(user_peak_i[, 4])) == nrow(user_peak_i)) {
            colname_new[4] <- "id"
            no_id <- FALSE
          }
        }
        colnames(user_peak_i) <- colname_new
        if (no_id) {
          user_peak_i$id <- paste0(
            user_peak_id[i], "_", as.vector(rownames(user_peak_i))
          )
        }
        peak_list_count <- peak_list_count + 1
        peak_list_all[[peak_list_count]] <- user_peak_i
        # test if user input id i match any TFregulomeR ID
        motif_matrix_i <- suppressMessages(.searchMotif(
          id = user_peak_id[i],
          motif_format = motif_type,
          api_object = api_object
        ))
        if (is.null(motif_matrix_i)) {
          is_TFregulome <- c(is_TFregulome, FALSE)
        } else {
          is_TFregulome <- c(is_TFregulome, TRUE)
        }
      }
    }
  } else {
    user_peak_id_new <- c()
  }
  # combine TFregulomeR ID and user ID
  peak_id_all <- c(TFregulome_peak_id, user_peak_id_new)

  # output
  peak_output <- list(
    peak_id_all = peak_id_all, peak_list_all = peak_list_all,
    is_TFregulome = is_TFregulome
  )
  return(peak_output)
}


fetch_peak_info <- function(id, peak, isTFregulome, api_object) {
  if (isTFregulome) {
    # make the request
    request_content_df <- apiRequest(api_object, id = id)
    source_i <- request_content_df[, "source"]
    if (source_i == "MethMotif") {
      isMethMotifID <- TRUE
      motif_seq_path <- request_content_df[1, c("TFBS")]
      meth_file_path <- request_content_df[1, c("DNA_methylation_profile")]
      meth_file_200bp_path <- request_content_df[1, c("DNA_methylation_profile_200bp")]
      WGBS_replicate <- request_content_df[1, c("WGBS_num")]
      peak_info <- list(
        isMethMotifID = isMethMotifID,
        motif_seq_path = motif_seq_path,
        meth_file_path = meth_file_path,
        meth_file_200bp_path = meth_file_200bp_path,
        WGBS_replicate = WGBS_replicate
      )
    } else {
      isMethMotifID <- FALSE
      motif_seq_path <- request_content_df[1, c("TFBS")]
      peak_info <- list(
        isMethMotifID = isMethMotifID,
        motif_seq_path = motif_seq_path
      )
    }

    # if peak x is from TFregulome database, extend peak regions by 100 bp
    peak$start <- peak$start - 99
    peak$end <- peak$end + 100
  } else {
    peak_info <- list(
      isMethMotifID = FALSE
    )
  }

  bed <- GenomicRanges::GRanges(
    peak$chr,
    IRanges::IRanges(peak$start, peak$end),
    id = peak$id
  )
  peak_info[["bed"]] <- bed

  return(peak_info)
}


intersect_peak_regions <- function(x_info, y_info, external_source_provided, external_source_grange, external_source_signal, motif_type, methylation_profile_in_narrow_region) {
  # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
  suppressWarnings(bedx_with_bedy <- subsetByOverlaps(x_info$bed, y_info$bed))
  peakx_with_peaky <- unique(as.data.frame(bedx_with_bedy))
  x_intersect_percentage <- 100 * nrow(peakx_with_peaky) / nrow(x_info$peak)
  MethMotif_x <- new('MethMotif')
  external_signal_in_x <- c(
    "signal_number" = NA, "mean" = NA, "SD" = NA,
    "median" = NA, "quartile_25" = NA, "quartile_75" = NA
  )
  tag_density_x <- c(
    "peak_number" = NA, "mean" = NA, "SD" = NA, "median" = NA,
    "quartile_25" = NA, "quartile_75" = NA
  )
  # collecting all CpG meth scores in the defined methylation profile area
  meth_score_collection_x <- data.frame()
  # methylation distribution only meaningful if we have WGBS and peaks
  is_methProfile_meaningful_x <- FALSE

  # form MethMotif object if the id is TFregulomeR id
  if (x_info$isTFregulome && nrow(peakx_with_peaky) > 0) {
    ############### profile score for the external source #####################
    if (external_source_provided) {
      suppressWarnings(external_signal_of_peakx_with_peaky_grange <- subsetByOverlaps(external_source_grange, bedx_with_bedy))
      external_signal_of_peakx_with_peaky <- unique(as.data.frame(external_signal_of_peakx_with_peaky_grange))
      if (nrow(external_signal_of_peakx_with_peaky) > 0) {
        external_signal_of_peakx_with_peaky_allInfo <- external_source_signal[which(external_source_signal$id %in% external_signal_of_peakx_with_peaky$id), ]
        external_signal_in_x <- c(0, 0, 0, 0, 0, 0)
        names(external_signal_in_x) <- c(
          "signal_number", "mean", "SD", "median", "quartile_25", "quartile_75"
        )
        if (nrow(external_signal_of_peakx_with_peaky_allInfo) > 0) {
          signal_num_x <- nrow(external_signal_of_peakx_with_peaky_allInfo)
          signal_mean_x <- mean(external_signal_of_peakx_with_peaky_allInfo$score)
          signal_sd_x <- sd(external_signal_of_peakx_with_peaky_allInfo$score)
          signal_median_x <- median(external_signal_of_peakx_with_peaky_allInfo$score)
          singal_quartile_x <-quantile(external_signal_of_peakx_with_peaky_allInfo$score)
          singal_quartile_25_x <- as.numeric(singal_quartile_x[2])
          singal_quartile_75_x <- as.numeric(singal_quartile_x[4])
          external_signal_in_x <- c(
            signal_num_x,signal_mean_x, signal_sd_x, signal_median_x,
            singal_quartile_25_x, singal_quartile_75_x
          )
        }
      }
    }
    ############### profile score for the external source #####################

    # tag fold change summary
    peakx_with_peaky_all_Info <- x_info$peak[(x_info$peak$id %in% peakx_with_peaky$id), ]
    if (ncol(peakx_with_peaky_all_Info) <= 4) {
      message(paste0("... ... ... Warning: ", x_info$id, " is derived from TFregulomeR but lack fifth 'read_fold_change' column. So skip read enrichment profiling !!"))
    } else if (!is.numeric(peakx_with_peaky_all_Info[, 5])) {
      message(paste0("... ... ... Warning: ", x_info$id, " is derived from TFregulomeR and fifth column is supposed to be 'tag_old_change' but it's NOT numeric . So skip read enrichment profiling !!"))
    } else {
      tag_x_num <- nrow(peakx_with_peaky_all_Info)
      tag_density_x_median <- median(peakx_with_peaky_all_Info[, 5])
      tag_density_x_mean <- mean(peakx_with_peaky_all_Info[, 5])
      tag_density_x_sd <- sd(peakx_with_peaky_all_Info[, 5])
      tag_density_x_quartile <- quantile(peakx_with_peaky_all_Info[, 5])
      tag_density_x_quartile_25 <- as.numeric(tag_density_x_quartile[2])
      tag_density_x_quartile_75 <- as.numeric(tag_density_x_quartile[4])
      tag_density_x <- c(
        tag_x_num, tag_density_x_mean, tag_density_x_sd, tag_density_x_median,
        tag_density_x_quartile_25,tag_density_x_quartile_75
      )
      names(tag_density_x) <- c(
        "peak_number", "mean", "SD", "median", "quartile_25", "quartile_75"
      )
    }

    # compute motif matrix
    motif_seq_x <- read.delim(x_info$motif_seq_path, sep = "\t", header = FALSE)
    colnames(motif_seq_x) <- c(
      "chr", "start", "end", "strand", "weight", "pvalue", "qvalue", "sequence"
    )
    motif_len_x <- nchar(as.character(motif_seq_x[1, "sequence"]))
    motif_seq_x$id <- paste0(
      x_info$id, "_motif_sequence_", as.vector(rownames(motif_seq_x))
    )
    motif_seq_x_grange <- GenomicRanges::GRanges(
      motif_seq_x$chr,
      IRanges::IRanges(
        motif_seq_x$start + 1,
        motif_seq_x$end
      ),
      id = motif_seq_x$id,
      pvalue = motif_seq_x$pvalue,
      sequence = motif_seq_x$sequence
    )
    suppressWarnings(motif_of_peakx_with_peaky_grange <- subsetByOverlaps(
      motif_seq_x_grange, bedx_with_bedy
    ))
    motif_of_peakx_with_peaky_overlap <- findOverlaps(
      motif_seq_x_grange, bedx_with_bedy, select = "first"
    )
    mcols(motif_seq_x_grange)$peak_id <- mcols(bedx_with_bedy)$id[motif_of_peakx_with_peaky_overlap]
    motif_seq_x_with_peak_id <- as.data.frame(motif_seq_x_grange)
    motif_of_peakx_with_peaky_df <- motif_seq_x_with_peak_id %>%
      dplyr::filter(!is.na(peak_id)) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::distinct(peak_id, .keep_all = TRUE)
    if (nrow(motif_of_peakx_with_peaky_df) > 0) {
      #nPeaks
      motif_matrix_of_peakx_with_y <- formMatrixFromSeq(
        input_sequence = as.vector(motif_of_peakx_with_peaky_df$sequence),
        motif_format = motif_type
      )

      # compute beta score matrix,
      if (x_info$isMethMotifID) {
        # methylation file can be empty
        meth_level_x <- tryCatch(
          read.delim(x_info$meth_file_path, sep = "\t", header = FALSE),
          error = function(e) data.frame()
        )
        # methylation file can be empty
        if (nrow(meth_level_x) == 0) {
          beta_score_matrix_of_peakx_with_y <- formBetaScoreFromSeq(
            input_meth = data.frame(),
            WGBS_replicate = x_info$WGBS_replicate,
            motif_len = motif_len_x
          )
        } else {
          colnames(meth_level_x) <- c(
            "chr", "start", "end", "meth_score", "C_num", "T_num", "seq_chr",
            "seq_start", "seq_end", "strand", "weight", "pvalue", "qvalue",
            "sequence"
          )
          meth_level_x$id <- paste0(
            x_info$id, "_motif_with_CG_", as.vector(rownames(meth_level_x))
          )
          meth_level_x_grange <- GenomicRanges::GRanges(
            meth_level_x$seq_chr,
            IRanges::IRanges(
              meth_level_x$seq_start,
              meth_level_x$seq_end
            ),
            id = meth_level_x$id
          )
          suppressWarnings(meth_level_x_with_y <- unique(as.data.frame(
            subsetByOverlaps(
              meth_level_x_grange, motif_of_peakx_with_peaky_grange
            )
          )))
          meth_level_x_with_y_allInfo <- meth_level_x[which(
            meth_level_x$id %in% meth_level_x_with_y$id
          ), ]
          beta_score_matrix_of_peakx_with_y <- formBetaScoreFromSeq(
            input_meth = meth_level_x_with_y_allInfo,
            WGBS_replicate = x_info$WGBS_replicate,
            motif_len = motif_len_x
          )
        }
      } else {
        beta_score_matrix_of_peakx_with_y <- as.matrix(NA)
      }

      if (motif_type == "TRANSFAC") {
        version_x <- 0
      } else {
        version_x <- 4
      }
      MethMotif_x@MMmotif <- updateMMmotif(
        MethMotif_x@MMmotif,
        motif_format = motif_type,
        version = version_x,
        background = c("A" = 0.25, "C" = 0.25, "G" = 0.25, "T" = 0.25),
        id = paste0(x_info$id, "_overlapped_with_", y_info$id),
        width = motif_len_x,
        nsites = length(motif_of_peakx_with_peaky_grange),
        nPeaks = length(bedx_with_bedy),
        motif_matrix = motif_matrix_of_peakx_with_y
      )
      MethMotif_x@MMBetaScore <- beta_score_matrix_of_peakx_with_y
    }

    # collecting CpG in x peaks that overlap with y peaks
    if (methylation_profile_in_narrow_region) {
      ### if in 200bp around peaks
      if (x_info$isMethMotifID) {
        is_methProfile_meaningful_x <- TRUE
        meth_level_200bp_x <- tryCatch(
          read.delim(x_info$meth_file_200bp_path, sep = "\t", header = FALSE),
          error = function(e) data.frame()
        )
        if (nrow(meth_level_200bp_x) > 0) {
          colnames(meth_level_200bp_x) <- c(
            "chr", "start", "end", "meth_score", "C_num", "T_num"
          )
          meth_level_200bp_x$id <- paste0(
            "200bp_CG_", as.vector(rownames(meth_level_200bp_x))
          )
          meth_level_200bp_x_grange <- GenomicRanges::GRanges(
            meth_level_200bp_x$chr,
            IRanges::IRanges(
              meth_level_200bp_x$start,
              meth_level_200bp_x$end
            ),
            id = meth_level_200bp_x$id
          )
          suppressWarnings(meth_level_in_peakx_200bp <- unique(as.data.frame(
            subsetByOverlaps(meth_level_200bp_x_grange, bedx_with_bedy)
          )))
          meth_level_in_peakx_200bp_allInfo <- unique(
            meth_level_200bp_x[which(
              meth_level_200bp_x$id %in% meth_level_in_peakx_200bp$id
            ), ]
          )
          meth_score_collection_x <- rbind(
            meth_score_collection_x,
            meth_level_in_peakx_200bp_allInfo[, c(
              "chr", "start", "end", "meth_score", "C_num", "T_num"
            )]
          )
        }
      }
    }
  }
  # form methylation score profile - distribution for peak x
  if (is_methProfile_meaningful_x) {
    if (nrow(meth_score_collection_x) > 0) {
      meth_score_distri_target_x <- formBetaScoreDistri(
        input_meth = as.data.frame(meth_score_collection_x$meth_score)
      )
    } else {
      meth_score_distri_target_x <- formBetaScoreDistri(
        input_meth = data.frame()
      )
    }
  } else {
    meth_score_distri_target_x <- matrix()
  }

  intersect_info <- list(
    intersect_percentage = x_intersect_percentage, MethMotif = MethMotif_x,
    meth_score_distri_target = meth_score_distri_target_x,
    external_signal_in = external_signal_in_x, tag_density = tag_density_x
  )
  return(intersect_info)
}


formMatrixFromSeq <- function(input_sequence, motif_format) {
  input_sequence <- as.data.frame(input_sequence)
  input_sequence_matrix <- data.frame(
    do.call(
      "rbind",
      strsplit(as.character(input_sequence$input_sequence), "", fixed = TRUE)
    )
  )

  motif_matrix_TRANSFAC <- matrix(
    rep(-1, 4 * ncol(input_sequence_matrix)),
    ncol = 4
  )
  alphabet <- c("A", "C", "G", "T")
  colnames(motif_matrix_TRANSFAC) <- alphabet

  for (i in seq(1, 4, 1)){
    for (j in seq(1, ncol(input_sequence_matrix), 1)) {
      motif_matrix_TRANSFAC[j, i] <- sum(
        input_sequence_matrix[, j] == alphabet[i]
      )
    }
  }
  motif_matrix_MEME <- motif_matrix_TRANSFAC / nrow(input_sequence_matrix)

  if (motif_format == "MEME") {
    return(motif_matrix_MEME)
  } else {
    return(motif_matrix_TRANSFAC)
  }
}


formBetaScoreFromSeq <- function(input_meth, WGBS_replicate, motif_len) {
  if (nrow(input_meth) == 0) {
    empty_matrix <- TRUE
  } else {
    if (WGBS_replicate == "2") {
      input_meth <- unique(input_meth[, c(
        "chr", "start", "end", "meth_score", "C_num", "T_num", "seq_chr",
        "seq_start", "seq_end", "strand", "sequence"
      )])
      input_meth_d <- input_meth[which(input_meth$strand == "+"), ]
      input_meth_r <- input_meth[which(input_meth$strand == "-"), ]

      input_meth_d$dis <- input_meth_d$start - input_meth_d$seq_start + 1
      input_meth_r$dis <- motif_len - (
        input_meth_r$start - input_meth_r$seq_start + 1
      )

      if (nrow(input_meth_d) == 0 && nrow(input_meth_r) == 0) {
        empty_matrix <- TRUE
      } else if (nrow(input_meth_d) == 0 && nrow(input_meth_r) != 0) {
        input_meth_sub <- input_meth_r[, c("dis", "meth_score")]
        empty_matrix <- FALSE
      }  else if (nrow(input_meth_d) != 0 && nrow(input_meth_r) == 0) {
        input_meth_sub <- input_meth_d[, c("dis", "meth_score")]
        empty_matrix <- FALSE
      } else if(nrow(input_meth_d) != 0 && nrow(input_meth_r) != 0) {
        input_meth_d_sub <- input_meth_d[,c("dis", "meth_score")]
        input_meth_r_sub <- input_meth_r[,c("dis", "meth_score")]
        input_meth_sub <- rbind(input_meth_d_sub, input_meth_r_sub)
        empty_matrix <- FALSE
      }
    } else {
      input_meth <- unique(input_meth[, c(
        "chr", "start", "end", "meth_score", "C_num", "T_num", "seq_chr",
        "seq_start", "seq_end", "strand", "sequence"
      )])
      input_meth_d <- input_meth[which(input_meth$strand == "+"), ]
      input_meth_r <- input_meth[which(input_meth$strand == "-"), ]

      input_meth_d$dis <- input_meth_d$start - input_meth_d$seq_start + 1
      input_meth_r$dis <- motif_len - (
        input_meth_r$start - input_meth_r$seq_start
      )
      # merge read in both strands
      if (nrow(input_meth_d) > 0) {
        for (i in seq(1, nrow(input_meth_d), 1)){
          if (unlist(strsplit(as.character(input_meth_d[i, "sequence"]), split = ""))[as.integer(input_meth_d[i, "dis"])] == "G") {
            input_meth_d[i, "dis"] <- input_meth_d[i, "dis"] - 1
          }
        }
      }
      # merge read in both strands
      if (nrow(input_meth_r) > 0) {
        for (i in seq(1, nrow(input_meth_r), 1)){
          if (unlist(strsplit(as.character(input_meth_r[i, "sequence"]), split = ""))[as.integer(input_meth_r[i, "dis"])] == "G"){
            input_meth_r[i, "dis"] <- input_meth_r[i, "dis"] - 1
          }
        }
      }
      # calculate overall beta score in both strand
      if (nrow(input_meth_d) > 0) {
        input_meth_d$id <- paste(input_meth_d$seq_chr,input_meth_d$seq_start,input_meth_d$seq_end,input_meth_d$dis,sep = "")
        input_meth_d_sub <- data.frame()
        input_meth_d_id_uniq <- unique(input_meth_d$id)
        for (i in input_meth_d_id_uniq){
          input_meth_d_temp <- input_meth_d[which(input_meth_d$id == i), c("C_num", "T_num", "dis")]
          dis_temp <- input_meth_d_temp[1, 3]
          meth_temp <- 100 * sum(input_meth_d_temp[, 1]) / sum(input_meth_d_temp[, c(1, 2)])
          new_add <- data.frame(i, dis_temp, meth_temp)
          input_meth_d_sub <- rbind(input_meth_d_sub, new_add)
        }
      }
      # calculate overall beta score in both strand
      if (nrow(input_meth_r) > 0) {
        input_meth_r$id <- paste(input_meth_r$seq_chr, input_meth_r$seq_start, input_meth_r$seq_end, input_meth_r$dis, sep = "")
        input_meth_r_sub <- data.frame()
        input_meth_r_id_uniq <- unique(input_meth_r$id)
        for (i in input_meth_r_id_uniq){
          input_meth_r_temp <- input_meth_r[which(input_meth_r$id == i), c("C_num", "T_num", "dis")]
          dis_temp <- input_meth_r_temp[1, 3]
          meth_temp <- 100 * sum(input_meth_r_temp[, 1]) / sum(input_meth_r_temp[, c(1, 2)])
          new_add <- data.frame(i, dis_temp, meth_temp)
          input_meth_r_sub <- rbind(input_meth_r_sub, new_add)
        }
      }

      if(nrow(input_meth_d)==0 && nrow(input_meth_r)==0){
        empty_matrix <- TRUE
      } else if(nrow(input_meth_d)==0 && nrow(input_meth_r)!=0){
        input_meth_sub <- input_meth_r_sub[,2:3]
        colnames(input_meth_sub) <- c("dis","meth_score")
        empty_matrix <- FALSE
      } else if(nrow(input_meth_d)!=0 && nrow(input_meth_r)==0){
        input_meth_sub <- input_meth_d_sub[,2:3]
        colnames(input_meth_sub) <- c("dis","meth_score")
        empty_matrix <- FALSE
      } else if(nrow(input_meth_d)!=0 && nrow(input_meth_r)!=0){
        input_meth_sub <- rbind(input_meth_d_sub[,2:3], input_meth_r_sub[,2:3])
        colnames(input_meth_sub) <- c("dis","meth_score")
        empty_matrix <- FALSE
      }

    }
  }

  if(empty_matrix==TRUE){
    plot_matrix <- data.frame(matrix(0,nrow = 3, ncol = motif_len))
    colnames(plot_matrix) <- c(seq(1,motif_len,1))
  } else{
    plot_matrix <- data.frame(c(1,2,3))
    for (i in seq(1,motif_len,1)){
      input_meth_sub_i <- input_meth_sub[which(input_meth_sub$dis == i),]
      unmeth_count <- nrow(input_meth_sub_i[which(input_meth_sub_i$meth_score<10),])
      meth_count <- nrow(input_meth_sub_i[which(input_meth_sub_i$meth>90),])
      inbetween_count <- nrow(input_meth_sub_i)-unmeth_count-meth_count
      new_column <- c(unmeth_count,inbetween_count,meth_count)
      plot_matrix <- data.frame(plot_matrix,new_column)

    }
    plot_matrix <- plot_matrix[,2:ncol(plot_matrix)]
    colnames(plot_matrix) <- c(seq(1,motif_len,1))
  }

  return(as.matrix(plot_matrix))
}
