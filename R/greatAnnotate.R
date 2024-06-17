#' greatAnnotate
#'
#' This function allows you to analyze the gene ontologies of targeting genes by cis-regulatory regions.
#' @param peaks Required. A bed-format genomic regions in data frame.
#' @param assembly The genome assembly of the input regions, currently supporting 'hg19', 'hg38' (default), 'mm9' and 'mm10'.
#' @param return_annotation Either TRUE of FALSE (default). If TRUE, a data.frame containing annotation results will be returned.
#' @param return_html_report Either TRUE of FALSE (default). If TRUE, a dynamic HTML report will be saved.
#' @param pvalue The adjusted p-value which is applied to filter the results, by default 0.01.
#' @param test The statistical test used in GREAT analysis, either 'binomial' (default) or 'hypergeometric'.
#' @param great_rule Equivalent to the rGREAT input 'rule', 'basalPlusExt' (default, basal plus extension), 'twoClosest' (two nearest genes), or 'oneClosest' (single nearest gene).
#' @param great_adv_upstream Equivalent to the rGREAT input 'adv_upstream' (upstream extension, kb). Only applicable when 'great_rule' is  'basalPlusExt', by default 5.0 (kb).
#' @param great_adv_downstream Equivalent to the rGREAT input 'adv_downstream' (downstream extension, kb). Only applicable when 'great_rule' is  'basalPlusExt', by default 1.0 (kb).
#' @param great_adv_span Equivalent to the rGREAT input 'adv_span' (maximal distal region). Only applicable when 'great_rule' is  'basalPlusExt', by default 1000.0 (kb).
#' @param great_adv_twoDistance Equivalent to the rGREAT input 'adv_twoDistance' (region range to be considered). Only applicable when 'great_rule' is  'twoClosest', by default 1000.0 (kb).
#' @param great_adv_oneDistance Equivalent to the rGREAT input 'adv_oneDistance' (region range to be considered). Only applicable when 'great_rule' is  'oneClosest', by default 1000.0 (kb).
#' @param request_interval The minimal gap time between two requests using greatAnnotate, by default 60 (s).
#' @param great_version Equivalent to the rGREAT input 'version', by default 4.0.4
#' @return  a data.frame, or an HTML report depending on the options.
#' @keywords greatAnnotate
#' @export
#' @examples
#' require(rGREAT)
#' require(rbokeh)
#' K562_CEBPB_regions <- loadPeaks(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_regions_annotation <- greatAnnotate(
#'   peaks = K562_CEBPB_regions[1:500, ],
#'   return_annotation = TRUE,
#'   return_html_report = TRUE
#' )

greatAnnotate <- function(
  peaks, assembly = "hg38",
  return_annotation = FALSE,
  return_html_report = FALSE,
  pvalue = 0.01,
  test = "binomial",
  great_rule = "basalPlusExt",
  great_adv_upstream = 5.0,
  great_adv_downstream = 1.0,
  great_adv_span = 1000.0,
  great_adv_twoDistance = 1000.0,
  great_adv_oneDistance = 1000.0,
  request_interval = 60,
  great_version = "4.0.4",
  cores = 1,
  local_version = FALSE
) {
  # check input arguments
  if (missing(peaks)) {
    stop("please provide peak regions using 'peaks ='!")
  }
  if (!is.data.frame(peaks)) {
    stop("The 'peaks' should be a BED-format data.frame!")
  }
  if (!(assembly %in% c("hg38", "hg19", "mm10", "mm9"))) {
    stop("Currently greatAnnotate only supports hg19, hg38, mm9 and mm10.")
  }
  if (!is.logical(return_annotation)) {
    stop("'return_annotation' should be either TRUE (T) or FALSE (F, default)")
  }
  if (!is.logical(return_html_report)) {
    stop("'return_html_report' should be either TRUE (T) or FALSE (F, default)")
  }
  if (!is.numeric(pvalue)) {
    stop("'pvalue' should be a numeric!")
  }
  if (!(test %in% c("binomial", "hypergeometric"))) {
    stop("'test' should be either 'binomial' (default) or 'hypergeometric'!")
  }
  if (!(great_rule %in% c("basalPlusExt", "twoClosest", "oneClosest"))) {
    stop("According to 'rGREAT', 'great_rule' shoule be 'basalPlusExt' (default), 'twoClosest', or 'oneClosest'!")
  }
  if (is.double(cores)) {
    cores <- as.integer(cores)
  } else if (!is.integer(cores)) {
    stop("'cores' should be an integer no greater than the core count of your CPU!")
  }
  if (!is.logical(local_version)) {
    stop("'local_version' should be either TRUE (T) or FALSE (F, default)")
  }
  # check loaded package
  if (!("rGREAT" %in% (.packages()))) {
    stop("GREAT R package 'rGREAT' (>=1.16.1) is NOT loaded yet!")
  }
  # if (return_html_report && !("rbokeh" %in% (.packages()))) {
  #   stop("A dynamic html report requires 'rbokeh' pakcage (>=0.5.0). 'rbokeh' is NOT loaded yet!")
  # }

  #message
  message("Start greatAnnotate ...")
  if (return_annotation) {
    message("... ... You chose to return annotated results in a data.frame.")
  } else {
    message("... ... You chose NOT to return annotated results in a data.frame.")
  }
  if (return_html_report) {
    message("... ... You chose to return an HTML report.")
  } else {
    message("... ... You chose NOT to return an HTML report.")
  }
  if (!return_annotation && ! return_html_report) {
    message("... ... You chose no action. EXIT!!")
    return(NULL)
  }

  peaks <- peaks[, c(1, 2, 3)]
  colnames(peaks) <- c("chr", "start", "end")
  peaks$id <- paste0("greatAnnotate_peak_", as.vector(rownames(peaks)))
  # great analysis
  message("... ... start GREAT analysis")
  if (local_version) {
    peaks.gr <- GenomicRanges::GRanges(
      peaks$chr,
      IRanges::IRanges(peaks$start, peaks$end),
      id = peaks$id
    )
    if (assembly %in% c("hg38", "hg19")) {
      gap_chromo <- paste0("chr", c(1:22, "X", "Y"))
    } else if (assembly %in% c("mm10", "mm9")) {
      gap_chromo <- paste0("chr", c(1:19, "X", "Y"))
    }
    gap <- rGREAT::getGapFromUCSC(assembly, gap_chromo)
    great_gene_sets <- c("GO:BP", "GO:CC", "GO:MF")
    # multiple by 1000 to convert from Kb to base pairs
    # local rgreat takes base pairs while the submission takes Kb
    jobs <- lapply(
      great_gene_sets, local_great, peaks = peaks.gr,
      assembly = assembly, great_rule = great_rule,
      great_adv_upstream = great_adv_upstream * 1000,
      great_adv_downstream = great_adv_downstream * 1000,
      great_adv_span = great_adv_span * 1000,
      great_adv_twoDistance = great_adv_twoDistance * 1000,
      great_adv_oneDistance = great_adv_oneDistance * 1000,
      cores = cores, gap = gap
    )
  } else {
    if (great_rule == "basalPlusExt") {
      jobs <- suppressMessages(rGREAT::submitGreatJob(
        gr = peaks, species = assembly, rule = great_rule,
        adv_upstream = great_adv_upstream,
        adv_downstream = great_adv_downstream, adv_span = great_adv_span,
        request_interval = request_interval, version = great_version
      ))
    } else if (great_rule == "twoClosest") {
      jobs <- suppressMessages(rGREAT::submitGreatJob(
        gr = peaks, species = assembly, rule = great_rule,
        adv_twoDistance = great_adv_twoDistance,
        request_interval = request_interval, version = great_version
      ))
    } else {
      jobs <- suppressMessages(rGREAT::submitGreatJob(
        gr = peaks, species = assembly, rule = great_rule,
        adv_oneDistance = great_adv_oneDistance,
        request_interval = request_interval, version = great_version
      ))
    }
  }

  message("... ... getting enrichment tables")
  if (local_version) {
    get_enrichment <- function(job) {
      suppressMessages(rGREAT::getEnrichmentTables(job))
    }
    tb <- lapply(jobs, get_enrichment)
    names(tb) <- c(
      "Biological Process", "Cellular Component", "Molecular Function"
    )
  } else {
    tb <- suppressMessages(rGREAT::getEnrichmentTables(jobs))
  }
  handle_enrichment_table <- function(gene_set, test, pvalue, local_version) {
    go_tb <- tb[[gene_set]]
    if (local_version) {
      if (test == "binomial") {
        go_tb <- go_tb %>%
          dplyr::select(id, description, observed_gene_hits, p_adjust)
      } else {
        go_tb <- go_tb %>%
          dplyr::select(id, description, observed_gene_hits, p_adjust_hyper) %>%
          dplyr::rename(p_adjust = p_adjust_hyper)
      }
    } else {
      if (test == "binomial") {
        go_tb <- go_tb %>%
          dplyr::select(ID, name, Hyper_Observed_Gene_Hits, Binom_Adjp_BH) %>%
          dplyr::rename(p_adjust = Binom_Adjp_BH)
      } else {
        go_tb <- go_tb %>%
          dplyr::select(ID, name, Hyper_Observed_Gene_Hits, Hyper_Adjp_BH) %>%
          dplyr::rename(p_adjust = Hyper_Adjp_BH)
      }
    }
    go_tb_sorted_pass <- go_tb %>%
      dplyr::filter(p_adjust < pvalue) %>%
      dplyr::arrange(p_adjust) %>%
      dplyr::mutate(go_id = gene_set) %>%
      dplyr::select(go_id, dplyr::everything())
    go_tb_sorted_pass
  }

  # added to align the names between local and submitted rgreat
  if (!local_version) {
    names(tb) <- substring(names(tb), first = 4)
  }
  filtered_tb <- lapply(
    names(tb), handle_enrichment_table, test = test, pvalue = pvalue,
    local_version = local_version
  )

  all <- do.call(rbind, filtered_tb)
  if (nrow(all) == 0) {
    message("... ... No result for the current settings!")
    return(NULL)
  } else {
    if (return_html_report) {
      formatted_time <- format(Sys.time(), "%Y-%m-%d_%H:%M:%S")
      report_prefix <- paste0("greatAnnotate_result-", formatted_time)
      report_name <- paste0(report_prefix, ".html")

      html_output <- formHTMLoutput(all, report_prefix)

      htmltools::save_html(html_output, report_name)
      message(paste0(
        "... ... An html report has been generated as '", report_name, "'!"
      ))
    }
    if (return_annotation) {
      colnames(all) <- c(
        "category", "ID", "name", "number_of_targeting_genes", "adjusted_pvalue"
      )
      message(
        "... ... The annotation results have been returned in a data.frame!"
      )
      return(all)
    }
  }
}

formHTMLoutput <- function(all, report_prefix) {
  colnames(all) <- c(
    "go_id", "id", "description", "observed_gene_hits", "p_adjust"
  )
  # required to get crosstalk between scatter plot and table
  shared_all <- crosstalk::SharedData$new(all)

  plot <- ggplot2::ggplot(
    shared_all,
    ggplot2::aes(
      x = -log10(p_adjust),
      y = observed_gene_hits,
      color = go_id,
      label1 = description,
      label2 = p_adjust,
      label3 = observed_gene_hits
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      # title = "GREAT analysis",
      x = "-log10(adjusted p-value)",
      y = "Number of genes",
      color = "GO Gene set"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  # plotly offers the interactive functionality on top of a ggplot
  plotly_plot <- plotly::ggplotly(
    plot,
    tooltip = c("label1", "label2", "label3"),
    height = 600, width = 600
  )

  filter_bar <- crosstalk::filter_checkbox(
    "go_id", "Go Gene set", shared_all, ~go_id, inline = TRUE
  )

  dt_table <- DT::datatable(
    shared_all,
    width = "100%",
    colnames = c(
      "Gene set", "ID", "Description", "Observed gene hits", "Adjusted p-value"
    ),
    rownames = FALSE,
    extensions = "Buttons",
    options = list(
      dom = "Bfrtip",
      pageLength = 25,
      buttons = list(
        "copy",
        list(
          extend = "collection",
          buttons = list(
            list(extend = "csv", title = report_prefix),
            list(extend = "excel", title = report_prefix)
          ),
          text = "Download"
        )
      )
    )
  )

  base_css <- "display: flex; justify-content: center; align-items: center;"

  html_page <- htmltools::tagList(
    htmltools::tags$head(htmltools::tags$title("Crosstalk Example")),
    htmltools::tags$body(
      htmltools::tags$div(
        style = paste0("flex-direction: column; ", base_css),
        htmltools::tags$h2("TFregulomeR - GreatAnnotate Result"),
        htmltools::tags$div(
          style = base_css,
          htmltools::tags$h3("Figure - Gene ontology analysis")
        ),
        htmltools::tags$div(style = base_css, plotly_plot),
        htmltools::tags$div(
          style = paste0("width: 75%; margin-top: 20px; ", base_css),
          filter_bar
        ),
        htmltools::tags$div(
          style = base_css,
          htmltools::tags$h3(
            "Table - Gene ontology analysis of targeting genes"
          )
        ),
        htmltools::tags$div(style = "width: 75%;", dt_table)
      )
    )
  )

  return(html_page)
}

local_great <- function(
  gene_set, peaks, assembly, great_rule, great_adv_upstream,
  great_adv_downstream, great_adv_span, great_adv_twoDistance,
  great_adv_oneDistance, cores, gap
) {
  if (great_rule == "basalPlusExt") {
    result <- rGREAT::great(
      gr = peaks, gene_sets = gene_set,
      tss_source = paste0("txdb:", assembly),
      mode = great_rule,
      basal_upstream = great_adv_upstream,
      basal_downstream = great_adv_downstream,
      extension = great_adv_span,
      cores = cores,
      exclude = gap
    )
  } else if (great_rule == "twoClosest") {
    result <- rGREAT::great(
      gr = peaks, gene_sets = gene_set,
      tss_source = paste0("txdb:", assembly),
      mode = great_rule,
      extension = great_adv_twoDistance,
      cores = cores,
      exclude = gap
    )
  } else {
    result <- rGREAT::great(
      gr = peaks, gene_sets = gene_set,
      tss_source = paste0("txdb:", assembly),
      mode = great_rule,
      extension = great_adv_oneDistance,
      cores = cores,
      exclude = gap
    )
  }
  n_gene_total <- length(unique(result@extended_tss$gene_id))
  result@table$p_value_hyper <- exp(phyper(
    result@table$observed_gene_hits - 1,
    result@table$gene_set_size,
    n_gene_total - result@table$gene_set_size,
    result@n_gene_gr,
    lower.tail = FALSE,
    log.p = TRUE
  ))
  result@table$p_value <- exp(pbinom(
    result@table$observed_region_hits - 1,
    result@n_total,
    result@table$genome_fraction,
    lower.tail = FALSE,
    log.p = TRUE
  ))
  return(result)
}
