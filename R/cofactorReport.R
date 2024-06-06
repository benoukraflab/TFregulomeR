#' cofactorReport
#'
#' This function allows you to get a PDF report of top cofactors along with DNA methylation for a TF.
#' @param intersectPeakMatrix Output of function 'intersectPeakMatrix()'.
#' @param top_num Number of highest co-binding factors to report for each TF (up to 10). By default the number is 10.
#' @param cobinding_threshold Only the co-factors with co-binding percentages more than this threshold value will be reported. By default the threshold is 0.05.
#' @return A PDF file
#' @keywords cofactorReport
#' @export
#' @examples
#' peak_id_x <- c("MM1_HSA_K562_CEBPB")
#' peak_id_y <- c("MM1_HSA_K562_CEBPD", "MM1_HSA_K562_ATF4")
#' intersect_output <- intersectPeakMatrix(peak_id_x=peak_id_x,
#'                                         motif_only_for_id_x=TRUE,
#'                                         peak_id_y=peak_id_y)
#' cofactorReport(intersectPeakMatrix = intersect_output)

cofactorReport <- function(
  intersectPeakMatrix,
  top_num = 10,
  cobinding_threshold = 0.05
) {
  # check input arguments
  if (missing(intersectPeakMatrix)) {
    stop("Please provide output of 'intersectPeakMatrix()' using 'intersectPeakMatrix ='!")
  }
  # check the validity of input intersectPeakMatrix
  if (class(intersectPeakMatrix[1, 1][[1]])[1] != "IntersectPeakMatrix") {
    stop("The input 'intersectPeakMatrix' is not valid. Please use the output of function 'intersectPeakMatrix()'")
  }
  # check other input arguments
  if (top_num <= 0) {
    stop("Invalid input for 'top_num'. It should be > 0 but <=10")
  }
  if (top_num > 10) {
    stop("This function only reports no more than 10 cofactors for each TF.")
  }
  if (cobinding_threshold <= 0 || cobinding_threshold > 1) {
    stop("'cobinding_threshold' should be >0 but <=1.")
  }

  #### start reporting ####
  message("Start cofactorReport ...")
  message(
    paste0("... The maximum number of cofactors to be reported is ", top_num)
  )
  message(paste0(
    "... The minimum percent of co-binding peaks for a cofactor is ",
    cobinding_threshold * 100, "%"
  ))
  message("... Each peak set derived from TFregulomeR compendium in 'peak list x' will be reported in an individual PDF file")

  for (i in seq(1, nrow(intersectPeakMatrix), 1)) {
    intersectPeakMatrix_i <- intersectPeakMatrix[i, , drop = FALSE]
    id_i <- rownames(intersectPeakMatrix_i)
    # judge whether peak set id_i is derived from TFregulomeR compendium
    is_from_TFregulomeR <- intersectPeakMatrix_i[1, 1][[1]]@isxTFregulomeID
    if (!is_from_TFregulomeR) {
      message(paste0(
        "... ... peak id '", id_i,
        "' in 'peak list x' is NOT a TFregulomeR ID. SKIP!!"
      ))
      next
    }
    message(paste0("... ... Start reporting peak id '", id_i, "' ..."))
    suppressMessages(intersectPeakMatrix_res_i <- intersectPeakMatrixResult(
      intersectPeakMatrix_i,
      return_intersection_matrix = TRUE,
      angle_of_matrix = "x"
    ))
    intersect_matrix_i <- intersectPeakMatrix_res_i$intersection_matrix
    intersect_matrix_t_i <- as.data.frame(t(intersect_matrix_i))

    intersect_matrix_filter_i <- intersect_matrix_t_i %>%
      dplyr::filter(!!as.name(id_i) >= cobinding_threshold * 100)

    if (nrow(intersect_matrix_filter_i) == 0) {
      message(paste0("... ... ... The number of cofactors passing 'cobinding_threshold' for peak id '", id_i, "' is 0. SKIP!!"))
      next
    }

    intersect_matrix_order_i <- intersect_matrix_filter_i %>%
      dplyr::arrange(!!as.name(id_i))
    if (nrow(intersect_matrix_filter_i) > top_num) {
      message(paste0("... ... ... The number of cofactors passing 'cobinding_threshold' for peak id '", id_i, "' is ", nrow(intersect_matrix_filter_i), ". Only top ", top_num, " will be selected."))
      intersect_matrix_len <- nrow(intersect_matrix_filter_i)
      intersect_matrix_order_i <- intersect_matrix_order_i %>%
        dplyr::slice((intersect_matrix_len - top_num + 1):intersect_matrix_len)
    }
    intersect_matrix_order_i <- intersect_matrix_order_i %>%
      dplyr::mutate(x = id_i) %>%
      tibble::rownames_to_column(var = "y") %>%
      dplyr::rename(value = !!as.name(id_i)) %>%
      dplyr::relocate(x, y, value)
    #### cobinding heatmap ####
    intersect_matrix_heatmap_i <- intersect_matrix_order_i %>%
      dplyr::mutate(
        new_y = paste(
          sapply(strsplit(y, "_"), function(x) tail(x, 1)),
          rev(rownames(intersect_matrix_order_i)),
          sep = "-"
        )
      ) %>%
      dplyr::mutate(
        new_x = sapply(strsplit(x, "_"), function(x) tail(x, 1))
      ) %>%
      dplyr::mutate(
        new_y = factor(new_y, levels = as.character(new_y))
      )
    cobinding_ylabel <- as.character(
      intersect_matrix_heatmap_i$new_y[
        rev(seq(1, nrow(intersect_matrix_heatmap_i), 1))
      ]
    )
    cobinding_ylabel_new <- paste0(
      as.character(intersect_matrix_heatmap_i$new_x),
      "\n+\n", cobinding_ylabel
    )
    colors_cobinding <- colorRampPalette(
      c("white", "#D46A6A", "#801515", "#550000")
    )(11)
    cobinding_gradientn <- ggplot2::scale_fill_gradientn(
      colours = colors_cobinding,
      breaks = c(seq(0, 100, length = 11)),
      limits = c(0, 100)
    )
    base_theme <- ggplot2::theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.position = "none"
    )
    p1 <- ggplot2::ggplot(
      intersect_matrix_heatmap_i,
      ggplot2::aes(
        x = new_x,
        y = new_y,
        fill = value
      )
    ) +
      ggplot2::geom_tile() +
      cobinding_gradientn +
      base_theme +
      ggplot2::theme(
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 5)
      ) +
      ggplot2::scale_y_discrete(
        breaks = cobinding_ylabel,
        labels = cobinding_ylabel_new
      )
    #### cobinding heatmap legend ####
    legend1_matrix <- as.data.frame(matrix(nrow = 11, ncol = 2))
    colnames(legend1_matrix) <- c("cobinding", "value")
    legend1_matrix[, 1] <- "x"
    legend1_matrix[, 2] <- seq(0, 100, 10)
    p_legend1 <- ggplot2::ggplot(
      legend1_matrix,
      ggplot2::aes(
        x = value,
        y = cobinding,
        fill = value
      )
    ) +
      ggplot2::geom_tile() +
      cobinding_gradientn +
      base_theme +
      ggplot2::theme(
        axis.text = element_text(size = 10)
      ) +
      ggplot2::scale_y_discrete(breaks = c("x"), labels = c("co-binding(%)"))
    #### motifs ####
    empty_background <- ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
    motif_plots <- list()
    count <- 1
    for (id_y in rev(intersect_matrix_heatmap_i$y)) {
      motif <- t(
        intersectPeakMatrix_i[id_i, id_y][[1]]@MethMotif_x@MMmotif@motif_matrix
      )
      top_margin <- 10
      bottom_margin <- 0
      if (count == length(intersect_matrix_heatmap_i$y)) {
        top_margin <- 0
        bottom_margin <- 10
      }

      motif_plot <- ggplot2::ggplot() +
        suppressWarnings(ggseqlogo::geom_logo(data = motif, method = "bits")) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(
            t = top_margin, r = 0, b = bottom_margin, l = 0, unit = "pt"
          )
        ) +
        empty_background

      motif_plots[[count]] <- motif_plot
      count <- count + 1
    }
    if (nrow(intersect_matrix_heatmap_i) == 1) {
      p2 <- gridExtra::arrangeGrob(
        grobs = motif_plots,
        layout_matrix = as.matrix(c(NA, 1, 2))
      )
    } else if (nrow(intersect_matrix_heatmap_i) == 2) {
      p2 <- gridExtra::arrangeGrob(
        grobs = motif_plots,
        layout_matrix = as.matrix(c(NA, 1, NA, NA, 2, 3))
      )
    } else if (nrow(intersect_matrix_heatmap_i) == 3) {
      p2 <- gridExtra::arrangeGrob(
        grobs = motif_plots,
        layout_matrix = as.matrix(
          c(NA, 1, 1, 1, NA,
            NA, 2, 2, 2, NA,
            NA, 3, 3, 3, 4)
        )
      )
    } else {
      p2 <- gridExtra::arrangeGrob(
        grobs = motif_plots,
        nrow = nrow(intersect_matrix_heatmap_i)
      )
    }

    #### methylation within motif heatmap ####
    meth_test <- intersectPeakMatrix_i[1, 1][[1]]@MethMotif_x@MMBetaScore[1, 1]
    if (!is.na(meth_test)) {
      methylation_matrix_inside <- intersect_matrix_heatmap_i %>%
        dplyr::select(new_y)

      extract_methylation_inside <- function(id_y, mm_object, id_i) {
        meth_matrix <- mm_object[id_i, id_y][[1]]@MethMotif_x@MMBetaScore
        if (sum(meth_matrix) > 0) {
          un_meth <- 100 * sum(meth_matrix[1, ]) / sum(meth_matrix)
          between <- 100 * sum(meth_matrix[2, ]) / sum(meth_matrix)
          meth <- 100 * sum(meth_matrix[3, ]) / sum(meth_matrix)
          meth_output <- c(un_meth, between, meth)
        } else {
          meth_output <- c(0, 0, 0)
        }
        meth_output
      }
      meth_inside_info <- lapply(
        intersect_matrix_order_i$y, extract_methylation_inside,
        mm_object = intersectPeakMatrix_i, id_i = id_i
      )
      meth_inside_info <- as.data.frame(do.call(rbind, meth_inside_info))
      meth_columns <- c("unmeth", "in-between", "meth")
      colnames(meth_inside_info) <- meth_columns

      methylation_matrix_inside <- cbind(
        methylation_matrix_inside, meth_inside_info
      )
      meth_levels <- as.character(meth_columns)
      methylation_matrix_inside <- methylation_matrix_inside %>%
        # rounds to the first decimal point
        tidyr::pivot_longer(
          meth_columns,
          names_to = "meth_interval",
          values_to = "value"
        ) %>%
        dplyr::mutate(
          meth_interval = factor(meth_interval, levels = meth_levels)
        )

      colors_meth <- colorRampPalette(
        c("white", "#7887AB", "#4F628E", "#2E4172", "#061539")
      )(11)
      meth_gradientn <- ggplot2::scale_fill_gradientn(
        colours = colors_meth,
        breaks = c(seq(0, 100, length = 11)),
        limits = c(0, 100)
      )
      p3 <- ggplot2::ggplot(
        methylation_matrix_inside,
        ggplot2::aes(
          x = meth_interval,
          y = new_y,
          fill = value
        )
      ) +
        ggplot2::geom_tile() +
        meth_gradientn +
        base_theme +
        ggplot2::theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6)
        )
    } else {
      message(paste(
        "... ... ... No methylation information within motif for id '",
        id_i, "', since it is not from MethMotif."
      ))
      p3 <- ggplot2::ggplot() + empty_background
    }

    #### methylation surrounding heatmap ####
    if (!is.na(intersectPeakMatrix_i[1, 1][[1]]@methylation_profile_x[1, 1])) {
      methylation_matrix <- intersect_matrix_heatmap_i %>%
        dplyr::select(new_y)

      extract_methylation <- function(id_y, mm_object, id_i) {
        meth_matrix <- mm_object[id_i, id_y][[1]]@methylation_profile_x
        if (sum(meth_matrix) > 0) {
          un_meth <- 100 * meth_matrix[1] / sum(meth_matrix)
          between <- 100 * sum(meth_matrix[2:9]) / sum(meth_matrix)
          meth <- 100 * meth_matrix[10] / sum(meth_matrix)
          meth_output <- c(un_meth, between, meth)
        } else {
          meth_output <- c(0, 0, 0)
        }
        meth_output
      }
      meth_info <- lapply(
        intersect_matrix_order_i$y, extract_methylation,
        mm_object = intersectPeakMatrix_i, id_i = id_i
      )
      meth_info <- as.data.frame(do.call(rbind, meth_info))
      colnames(meth_info) <- meth_columns

      methylation_matrix <- cbind(methylation_matrix, meth_info)
      methylation_matrix <- methylation_matrix %>%
        # rounds to the first decimal point
        tidyr::pivot_longer(
          meth_columns,
          names_to = "meth_interval",
          values_to = "value"
        ) %>%
        dplyr::mutate(
          meth_interval = factor(meth_interval, levels = meth_levels)
        )
      p4 <- ggplot2::ggplot(
        methylation_matrix,
        ggplot2::aes(
          x = meth_interval,
          y = new_y,
          fill = value
        )
      ) +
        ggplot2::geom_tile() +
        meth_gradientn +
        base_theme +
        ggplot2::theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6)
        )
    } else {
      message(paste("... ... ... No methylation information in 200bp peak regions for id '", id_i, "', since it is not from MethMotif or 'methylation_profile_in_narrow_region=FALSE' during intersectPeakMatrix()."))
      p4 <- ggplot2::ggplot() + empty_background
    }

    #### methylation legend ####
    p_legend2 <- ggplot2::ggplot(
      legend1_matrix,
      ggplot2::aes(
        x = value,
        y = cobinding,
        fill = value
      )
    ) +
      ggplot2::geom_tile() +
      meth_gradientn +
      base_theme +
      ggplot2::theme(
        axis.text = element_text(size = 10)
      ) +
      ggplot2::scale_y_discrete(breaks = c("x"), labels = c("CG percent(%)"))

    #### read enrichment score ####
    extract_fold_change <- function(quartile) {
      suppressMessages(tag_res_i <- intersectPeakMatrixResult(
        intersectPeakMatrix_i,
        return_tag_density = TRUE,
        tag_density_value = quartile
      ))
      tag_i <- tag_res_i$tag_density_matrix
      tag_i_t <- as.data.frame(t(tag_i))
      tag_i_t_filter <- tag_i_t[
        intersect_matrix_heatmap_i$y,
        unique(intersect_matrix_heatmap_i$x),
        drop = FALSE
      ]
      colnames(tag_i_t_filter) <- "y"
      tag_i_t_filter$x <- seq(1, nrow(tag_i_t_filter), 1)
      tag_i_t_filter
    }
    tag_cases <- c("median", "quartile_25", "quartile_75")
    tag_results_i <- lapply(tag_cases, extract_fold_change)
    names(tag_results_i) <- c("median", "q1", "q3")

    tag_max_value <- max(c(
      tag_results_i$median$y,
      tag_results_i$q1$y,
      tag_results_i$q3$y
    ))
    tag_max_value_int <- (as.integer(tag_max_value / 10) + 1) * 10
    tag_x_max_list <- c(1, 0.7, 0.5, 0.5, 0.4, 0.4, 0.3, 0.3, 0, 0)
    tag_x_min <- 1 - tag_x_max_list[nrow(tag_results_i$median)]
    tag_x_max <- nrow(tag_results_i$median) + tag_x_max_list[
      nrow(tag_results_i$median)
    ]

    tag_quartile_total1 <- rbind(tag_results_i$q1, tag_results_i$q3)
    colnames(tag_results_i$q1) <- c("y1", "x1")
    colnames(tag_results_i$q3) <- c("y3", "x3")
    tag_quartile_total2 <- cbind(tag_results_i$q1, tag_results_i$q3)

    p_tag <- ggplot2::ggplot(
      tag_results_i$median,
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_point(shape = 19, size = 3) +
      ggplot2::ylim(0, tag_max_value_int) +
      ggplot2::xlim(tag_x_min, tag_x_max) +
      ggplot2::theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6)
      ) +
      ggplot2::geom_point(
        data = tag_quartile_total1,
        shape = "|",
        size = 3,
        mapping = ggplot2::aes(x = x, y = y)
      ) +
      ggplot2::geom_segment(
        data = tag_quartile_total2,
        mapping = ggplot2::aes(
          x = x1, y = y1,
          xend = x3, yend = y3
        )
      ) +
      ggplot2::labs(y = "Read enrichment score") +
      ggplot2::coord_flip()


    #### Final arrangement ####
    text1 <- grid::textGrob("cofactors")
    text2 <- grid::textGrob("motif")
    text3 <- grid::textGrob("methylation in motif")
    text4 <- grid::textGrob("methylation in\n 200bp regions")
    text5 <- grid::textGrob("read enrichment")

    pdf_i_name <- paste0(id_i, "_cofactor_report.pdf")
    pdf(pdf_i_name)
    gridExtra::grid.arrange(
      text1, text2, text3, text4, text5,
      p1, p2, p3, p4, p_tag, p_legend1, p_legend2,
      layout_matrix = rbind(
        c(1, 2, 2, 3, 3, 4, 4, 5, 5),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(6, 7, 7, 8, 8, 9, 9, 10, 10),
        c(NA, 11, 11, 11, NA, 12, 12, 12, NA)
      )
    )
    dev.off()
    message(paste0(
      "... ... ... Cofactor report for id '", id_i,
      "' has been saved as ", pdf_i_name
    ))
  }
}
