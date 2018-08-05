# Main functions to create the scatter plot

#' Calculate the rate of transcription and translation.
#'
#' Integrates the count from multiple sequencing experiments to produce the
#' rates of transcription and translation in contrast to the egg developmental
#' stage.
#'
#' @param stage Developmental stage; \code{string}.
#' @param totalCounts Read counts from total RNA sequencing experiment.
#' @param boundCounts Read counts form polysome profiling bound fraction
#' sequencing.
#' @param unboundCounts Read counts form polysome profiling unbound fraction
#' sequencing.
#' @param totalFPKM FPKM values from total RNA sequencing experiment.
#' @param cluster Clustering information for RNA fractions.
#' @param cutoff FPKM cutoff value.
#'
#' @return A data frame with rates and clustering information
#'
#' @seealso Materials and Methods, Winata CL, Development, 11 December 2017
#' \url{http://dev.biologists.org/content/early/2017/12/08/dev.159566}
#'
#' @export
getRTvsRB_FPKM <- function(stage,
                           totalCounts,
                           boundCounts,
                           unboundCounts,
                           totalFPKM,
                           cluster,
                           cutoff = 1) {
    # print stage
    cat("Stage:", stage, "\n")
    # set the pattern
    pat <- sprintf('(Egg)|(%s)', stage)
    # get the FPKM for the indicated stage
    totalStageVSEgg_FPKM <-
        totalFPKM[, grep(pat, colnames(totalFPKM))]
    # subset the genes, based on the FPKM higher than the indicated cutoff value
    totalStageVSEgg_FPKM_Cut <-
        totalStageVSEgg_FPKM[apply(
            totalStageVSEgg_FPKM,
            MARGIN = 1,
            FUN = function(x) {
                all(x > cutoff)
            }
        ),]

    # print the number of genes after cutoff
    cat(sprintf(
        'Number of genes with FPKM > %i in total fractions: %i\n',
        cutoff,
        dim(totalStageVSEgg_FPKM_Cut)[1]
    ))
    # subset the bound counts for the indicated stage with the genes names
    # from the total subset, and column names for the indicated stage
    boundStageVSEgg_Cut <-
        boundCounts[rownames(totalStageVSEgg_FPKM_Cut),
                    grep(pat, colnames(boundCounts))]
    # get rid of genes with 0 counts in bound fraction
    boundStageVSEgg_Cut <-
        boundStageVSEgg_Cut[!apply(
            boundStageVSEgg_Cut,
            MARGIN = 1,
            FUN = function(x) {
                any(x == 0)
            }
        ),]
    cat(sprintf(
        'Without genes with 0 counts in bound fraction: %i\n',
        dim(boundStageVSEgg_Cut)[1]
    ))
    # subset the genes from the unbound counts
    unboundStageVSEgg_Cut <-
        unboundCounts[rownames(boundStageVSEgg_Cut),
                      grep(pat, colnames(unboundCounts))]
    # get rid of those with 0 counts in unbound fraction
    unboundStageVSEgg_Cut <-
        unboundStageVSEgg_Cut[!apply(
            unboundStageVSEgg_Cut,
            MARGIN = 1,
            FUN = function(x) {
                any(x == 0)
            }
        ),]
    # subset those genes from the final total counts and bound counts
    totalStageVSEgg_Cut <-
        totalCounts[rownames(unboundStageVSEgg_Cut),
                    grep(pat, colnames(totalCounts))]
    boundStageVSEgg_Cut <-
        boundStageVSEgg_Cut[rownames(unboundStageVSEgg_Cut),]
    # print the number of genes
    cat(sprintf(
        'Without genes with 0 in bound and unbound fraction: %i\n',
        dim(totalStageVSEgg_Cut)[1]
    ))
    # calculate and return RB and RT ratios
    RT <-
        totalStageVSEgg_Cut[, grep(stage, colnames(totalStageVSEgg_Cut))] /
        totalStageVSEgg_Cut[, grep('Egg', colnames(totalStageVSEgg_Cut))]
    RB <-
        boundStageVSEgg_Cut[, grep(stage, colnames(boundStageVSEgg_Cut))] /
        boundStageVSEgg_Cut[, grep('Egg', colnames(boundStageVSEgg_Cut))]
    RU <-
        unboundStageVSEgg_Cut[, grep(stage, colnames(unboundStageVSEgg_Cut))] /
        unboundStageVSEgg_Cut[, grep('Egg', colnames(unboundStageVSEgg_Cut))]
    RBvsRTvsRU <- cbind(RB, RT, RU)
    # change the gene names so they are compatible with clustering gene names
    rownames(RBvsRTvsRU) <- sub('(ENSDARG\\d+).\\d+\\b',
                                '\\1',
                                rownames(RBvsRTvsRU),
                                perl = TRUE)
    RBvsRTvsRU <- as.data.frame(RBvsRTvsRU)
    # merge the ratio and clustering tables
    RBvsRTvsRU_clust <- merge(
        RBvsRTvsRU,
        cluster,
        by.x = 0,
        by.y = 1,
        all.x = TRUE
    )
    # print the number of unclustered genes
    cat(sprintf(
        'Number of unclustered genes: %i\n',
        length(RBvsRTvsRU_clust$cluster[is.na(RBvsRTvsRU_clust$cluster)])
    ))
    # label the unclustered genes
    levels(RBvsRTvsRU_clust$cluster) <-
        c(levels(RBvsRTvsRU_clust$cluster),
          'unclustered')
    RBvsRTvsRU_clust[is.na(RBvsRTvsRU_clust)] <- 'unclustered'
    RBvsRTvsRU_sanity <-
        cbind(RBvsRTvsRU_clust,
              sanity = apply(RBvsRTvsRU_clust,
                             MARGIN = 1,
                             function(x) {
                                 ((as.numeric(x[2]) <= as.numeric(x[3])) &
                                      (as.numeric(x[3]) <= as.numeric(x[4]))) |
                                     ((as.numeric(x[2]) > as.numeric(x[3])) &
                                          (as.numeric(x[3]) > as.numeric(x[4])))
                             }))
    cat(sprintf(
        'Number of genes passing sanity test: %i\n',
        length(RBvsRTvsRU_sanity$sanity[RBvsRTvsRU_sanity$sanity])
    ))
    # return the table
    return(RBvsRTvsRU_sanity)
}


#' Prepare the data for plotting
#'
#' Filter maternal genes and subset the rows to match the limits of the graph
#'
#' @param data A data frame for operations
#' @param xlim,ylim X and Y axis limits.
#' @param title Dataset name
#'
#' @return A processed data frame
#'
#' @keywords internal
prepare_data <- function(data, xlim, ylim, title) {
    # Filter only maternal genes
    data_m <-
        data %>% dplyr::filter(.data$cluster == "maternal CPA" |
                                   .data$cluster == "maternal NoCPA")
    # Arrange the data for plotting and subset the rows for those matching the
    # graph limits
    data_lim <- data_m %>% arrange(desc(.data$cluster)) %>%
        mutate(lRB = log2(.data$RB)) %>% mutate(lRT = log2(.data$RT)) %>%
        dplyr::filter(.data$lRB > xlim[1] & .data$lRB < xlim[2]) %>%
        dplyr::filter(.data$lRT > ylim[1] & .data$lRT < ylim[2])
    # Print how many genes were ommited by the graph scale
    cat(sprintf(
            '%s\n%i genes were filtered according to the graph limits criteria.\n\n',
            title,
            nrow(data_m) - nrow(data_lim)
        )
    )
    return(data_lim)
}


#' Plot the scatterplot
#'
#' Uses ggplot2 and ggExtra to produce the scatterplot with density plots on the
#' margins.
#'
#' @param data A data frame for plotting.
#' @param xlim,ylim X and Y axis limits.
#' @param title Plot title.
#' @param type Type of marginal density plot passed to ggMarginal function.
#' @param bins Number of histogram bins.
#' @param colour Colour values for different clusters.
#'
#' @return A ggplot2 plot object.
#'
#' @export
#' @importFrom dplyr filter arrange desc mutate
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom rlang .data
plotLogMat <- function(data,
                       xlim = c(-5, 5),
                       ylim = c(-5, 5),
                       title = NULL,
                       type = "histogram",
                       bins = 100,
                       colour = c("#ef8a62", "#67a9cf")) {
    data_lim <- prepare_data(data, xlim, ylim, title)
    p <- ggplot() +
        geom_point(
            data = data_lim,
            mapping = aes_string(x = 'lRB', y = 'lRT'),
            alpha = 0
        ) +
        ggthemes::geom_rangeframe(data = data.frame(x = xlim, y = ylim),
                                  aes_string(x = 'x', y = 'y')) +
        ggthemes::theme_tufte() +
        theme(text = element_text(size = 18)) +
        xlim(xlim) +
        ylim(ylim) +
        geom_polygon(aes(
            x = c(
                min(xlim[2], ylim[2]),
                max(xlim[1], ylim[1]),
                max(xlim[1], ylim[1]),
                xlim[2],
                xlim[2]
            ),
            y = c(
                min(xlim[2], ylim[2]),
                max(xlim[1], ylim[1]),
                ylim[1],
                ylim[1],
                min(xlim[2], ylim[2])
            )
        ),
        fill = 'grey20',
        alpha = 0.1) +
        geom_point(data = data_lim,
                   mapping = aes_string(x = 'lRB',
                                        y = 'lRT',
                                        colour = 'cluster')) +
        scale_color_manual(values = colour) +
        geom_vline(
            xintercept = 0,
            colour = 'black',
            linetype = 5,
            size = 0.5,
            alpha = 1
        ) +
        geom_hline(
            yintercept = 0,
            colour = 'black',
            linetype = 5,
            size = 0.5,
            alpha = 1
        ) +
        ggtitle(title) +
        xlab('log2(\u2206TL)') +
        ylab('log2(\u2206TS)') +
        theme(legend.position = 'bottom')
    if (type == "histogram") {
        return(
            ggExtra::ggMarginal(
                p,
                data = data_lim,
                type = type,
                bins = bins,
                fill = 'transparent'
            )
        )
    } else {
        return(ggExtra::ggMarginal(p,
                                   data = data_lim,
                                   type = type))
    }
}
