# Append to a ggplot2 object to add background transparency
transparent_background <- theme(#legend.position='bottom',
                                legend.background = element_rect(fill='transparent', colour=NA),
                                legend.key = element_rect(fill='transparent', colour=NA),
                                legend.box.background = element_rect(fill='transparent', colour=NA),
                                panel.background = element_rect(fill='transparent', colour=NA),
                                plot.background = element_rect(fill='transparent', colour=NA))

# Function to obtain default ggplot colour scheme
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Sharing legend from here: https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
# Which took it from https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(plt) {
  tmp <- ggplot_gtable(ggplot_build(plt))
  g_legend_from_gtable(tmp)
}

g_legend_from_gtable <- function(gtable) {
  leg <- which(sapply(gtable$grobs, function(x) x$name) == "guide-box")
  gtable$grobs[[leg]]
}


heatmap_mutation_extended <- function(df, genes, group_var, q=0.05, prob_var=NULL,
                                      idcol = 'ID',
                                      bar_vars=NULL, plot_all_genes=FALSE, colour_schemes=NULL,
                                      gene_text_size=6, group_text_size=8,
                                      q_axis_text_size=6, q_axis_label_text_size=10,
                                      legend_title_size=10, legend_text_size=8, legend_margin=-25,
                                      y_order=c('dendrogram', 'group'),
                                      bar_position=c('top', 'bottom')) {
    y_order <- match.arg(y_order)
    bar_position <- match.arg(bar_position)
    if (!idcol %in% colnames(df)) {
        message("Error: df must contain `idcol`.")
    }

    if (is.null(colour_schemes)) {
        colour_schemes <- list()
    }

    if (!group_var %in% colnames(df)) {
        message("Error: df must contain group_var.")
    }

    if (!all(sapply(genes, function(x) x %in% colnames(df)))) {
        message("Error: df must contain all columns in genes.")
    }

    if (is.null(prob_var)) {
        prob_var <- 'foo'
        df[[prob_var]] <- 0
    }

    df <- df %>%
        dplyr::rename(cluster = group_var, prob = prob_var, ID=idcol)
    cluster_levels <- levels(df$cluster)
    nclust <- length(cluster_levels)

    if (is.null(bar_vars)) {
        df_long <- df %>%
            dplyr::select(ID, cluster, prob, dplyr::one_of(genes)) %>%
            tidyr::gather(gene, is_mut, -ID, -cluster, -prob)
    } else {
        df_long <- df %>%
            dplyr::select(ID, cluster, prob, dplyr::one_of(genes), dplyr::one_of(bar_vars)) %>%
            tidyr::gather(gene, is_mut, -ID, -cluster, -prob, -dplyr::one_of(bar_vars))
    }

    # Add factor level for all clusters + no mutation (which will be white)
    df_long <- df_long %>%
        dplyr::mutate(mut_cluster=factor(ifelse(is_mut == 1, as.character(cluster), 'none'),
                                         levels=c('none', cluster_levels)))


    chisq_results <- as.data.frame(stats::setNames(lapply(cluster_levels, function(clust) {
        df <- df_long %>%
            dplyr::mutate(in_group = as.factor(as.numeric(cluster == clust)))

        sapply(genes, function(g) {
            df %>%
                dplyr::filter(gene == g) %>%
                droplevels() %>%
                dplyr::summarise(stats::chisq.test(is_mut, in_group)$p.value) %>%
                unlist()
        })

    }), cluster_levels))

    sig_long <- chisq_results %>%
        dplyr::mutate(gene = genes) %>%
        tidyr::gather(cluster, raw_p, -gene) %>%
        dplyr::mutate(cluster = factor(gsub("\\.", "/", cluster), levels=cluster_levels),
                      adj_p = stats::p.adjust(raw_p, method='BH'),
                      is_sig = adj_p< q) %>%
        dplyr::left_join(df_long %>% dplyr::group_by(gene) %>% dplyr::summarise(overall_prop = mean(is_mut == 1)), by='gene') %>%
        dplyr::left_join(df_long %>% dplyr::group_by(cluster, gene) %>% dplyr::summarise(cluster_prop = mean(is_mut == 1)), by=c('cluster', 'gene')) %>%
        dplyr::mutate(is_increase = cluster_prop > overall_prop,
               to_plot = is_sig & is_increase)

    # Filter to significant genes that have increase over general population frequency that will be plotted
    if (! plot_all_genes) {
        sig_long <- sig_long %>%
            dplyr::filter(to_plot)
    }

    scale_x_reordered <- function(..., sep = "___") {
        reg <- paste0(sep, ".+$")
        ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
    }

    if (nrow(sig_long) == 0) {
        message("No significantly enriched mutations detected.")
        return()
    }

    # Grab the genes to use for heatmap
    genes_to_plot <- sig_long %>%
                    dplyr::distinct(gene) %>%
                    dplyr::pull(gene)

    # And order them for both plots
    if (y_order == 'dendrogram') {
        cor_dist_muts <- stats::as.dist(philentropy::distance(t(df[genes_to_plot]), method='dice'))
        cor_dist_muts[is.na(cor_dist_muts)] <- 0
        gene_clust <- stats::hclust(cor_dist_muts, method='ward.D2')
        gene_order <- genes_to_plot[gene_clust$order]
    } else if (y_order == 'group') {
        gene_order <- sig_long %>%
                        dplyr::group_by(gene) %>%     # If have a gene sig for 2 clusters, then order by most sig
                        dplyr::top_n(1, -adj_p) %>%
                        dplyr::arrange(dplyr::desc(cluster), dplyr::desc(adj_p)) %>%
                        dplyr::pull(gene)
    } else {
        stop(paste0("Unknown y_order '", y_order, "'. Please select from 'dendrogram' and 'group'."))
    }

    df_long <- df_long %>%
        dplyr::filter(gene %in% genes_to_plot) %>%
        dplyr::mutate(gene = factor(clean_gene_names(as.character(gene)), levels=clean_gene_names(gene_order)))

    sig_long <- sig_long %>%
                    dplyr::mutate(gene = factor(clean_gene_names(as.character(gene)), levels=clean_gene_names(gene_order)))

    # Identify colour scheme to use for this many clusters
    if (is.null(colour_schemes[[group_var]])) {
        colours_main <- gg_color_hue(nclust)
    } else {
        colours_main <- colour_schemes[[group_var]]
    }

    # Turn into plot
    plt_q <- sig_long %>%
        dplyr::mutate(val = -log10(adj_p)) %>%
        ggplot2::ggplot(ggplot2::aes(x=gene,  y=val, fill=cluster)) +
            ggplot2::geom_col(position='dodge', width=0.5) +
            ggplot2::labs(y=latex2exp::TeX("$-\\log_{10}(q)"), x="") +
            scale_x_reordered() +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.major.y = ggplot2::element_blank(),
                legend.position = 'bottom',
                axis.text.x = ggplot2::element_text(size=q_axis_text_size),
                axis.title.x = ggplot2::element_text(size=q_axis_label_text_size),
                panel.border = ggplot2::element_blank()
            ) +
            ggplot2::scale_y_continuous(expand=c(0, 0)) +
            ggplot2::scale_fill_manual("", drop=FALSE, guide=F,
                                       values=colours_main) +
            transparent_background

    # Order patients by their probability of being in their assigned cluster
    if (!is.null(prob_var)) {
        pat_order_raw <- df_long %>%
            dplyr::distinct(ID, cluster, prob) %>%
            dplyr::mutate(pat_order = -prob) %>%
            dplyr::select(ID, cluster, pat_order)
    } else {
        pat_order_raw <- dplyr::bind_rows(lapply(stats::setNames(cluster_levels, cluster_levels), function(in_clust) {
            clust_df <- df_long %>%
                dplyr::filter(cluster == in_clust) %>%
                dplyr::select(ID, gene, is_mut) %>%
                tidyr::spread(gene, is_mut)

            if (length(unique(clust_df$ID)) > 2) {
                pat_dist <- stats::as.dist(philentropy::distance(clust_df[, -1], method='dice'))
                pat_dist[is.na(pat_dist)] <- 0
                pat_clust <- stats::hclust(pat_dist, method='ward.D2')
                data.frame(ID=clust_df$ID,
                           pat_order=pat_clust$order)
            } else {
                data.frame(ID=clust_df$ID,
                           pat_order=1)
            }

        }), .id='cluster')
    }

    pat_order_raw$cluster <- factor(pat_order_raw$cluster, levels=cluster_levels)
    patient_order <- pat_order_raw %>%
        dplyr::arrange(cluster, pat_order) %>%
        dplyr::pull(ID)
    df_long$ID <- factor(df_long$ID, levels=patient_order)

    plt_heatmap <- ggplot2::ggplot() +
            ggplot2::geom_tile(ggplot2::aes(x=ID, y=gene, fill=mut_cluster), data=df_long) +
            ggplot2::scale_fill_manual(guide=F, values=c('white', colours_main), drop=F) +
            ggplot2::scale_y_discrete(position="right") +
            ggplot2::scale_alpha_continuous(guide=F) +
            ggplot2::facet_grid(cols=ggplot2::vars(cluster), space='free_x', scales='free_x') +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_text(size=gene_text_size),
                           strip.text.x = ggplot2::element_text(size=group_text_size),
                           axis.ticks = ggplot2::element_blank(),
                           axis.title = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank(),
                           panel.spacing.x = ggplot2::unit(0, "lines"),
                           axis.line = ggplot2::element_blank()) +
            transparent_background

    # Add in facet background colours, removing white
    plt_heatmap_coloured <- colour_facets(plt_heatmap, colours_main)

    # Combine heatmap and q plot together
    comb_plot <- combine_heatmap_q_plots(plt_heatmap_coloured, plt_q)

    # Then deal with the top horizontal bar
    bar_grobs <- lapply(bar_vars, create_top_bar, df_long, colour_schemes, legend_title_size, legend_text_size, legend_margin)
    for (grb in bar_grobs) {
        comb_plot <- add_topbar_grob(comb_plot, grb, bar_position)
    }

    comb_plot
}

# Cleans gene names by removing the trailing _S indicative of somatic hypermutation and replaces
# my notation of _OR_del with _HD instead.
# x: A character vector of gene names to be replaced.
clean_gene_names <- function(x) {
    gsub("_S$", "", gsub("_OR_del$", "_HD", x))
}

# Change facet colours of a plot, taken from
# https://github.com/tidyverse/ggplot2/issues/2096#issuecomment-389825118
colour_facets <- function(plt, colours) {
    bld <- ggplot2::ggplot_build(plt)
    g <- ggplot2::ggplot_gtable(bld)
    strip_right <- which(grepl('strip-t-', g$layout$name))

    if (length(colours) < length(strip_right)) {
        stop("Error: have got fewer colours than facets")
    }

    for (i in seq_along(strip_right)) {
        strp <- strip_right[i]
        j <- which(grepl('rect', g$grobs[[strp]]$grobs[[1]]$childrenOrder))
        g$grobs[[strp]]$grobs[[1]]$children[[j]]$gp$fill <- colours[i]
    }
    g
}

# Combines a heatmap and a q plot into a single grob
# heatmap: A heatmap, either as ggplot object or grob
# qplot: The q-plot, either as ggplot object or grob
combine_heatmap_q_plots <- function(heatmap, qplot) {
    if (! 'gtable' %in% class(heatmap)) {
        heatmap <- ggplot2::ggplotGrob(heatmap)
    }
    if (! 'gtable' %in% class(qplot)) {
        qplot <- ggplot2::ggplotGrob(qplot)
    }

    # Add empty column onto heatmap where can put these panels
    comb <- gtable::gtable_add_cols(heatmap, ggplot2::unit(0.2, "npc"))

    # Add panel, axis, and xlabels
    comb <- add_panel(comb, qplot, "panel", l=ncol(comb), r=ncol(comb))
    comb <- add_panel(comb, qplot, "axis-b", set_height=TRUE, l=ncol(comb), r=ncol(comb))
    comb <- add_panel(comb, qplot, "xlab-b", set_height=TRUE, l=ncol(comb), r=ncol(comb))
    comb
}

# Adds various components from a right-hand plot into a the
# furtherest column of a left-hand plot.
# Assumes 't' and 'b' components are always the same which is a fair assumption
# left: Grob for left plot
# right: Grob for right plot
# component: String detailing which component to merge over, i.e. 'panel'
# set_height_of_left: Whether to set the height of this component on the left
#   hand plot to that of the same component on the right. Useful when adding
#   a component that wasn't used on left but was on right, and so if used default
#   left plot's sizing then it would be too small.
add_panel <- function(old, new, component, set_height=FALSE, set_width=FALSE,
                      l=NULL, r=NULL, t=NULL, b=NULL) {

    # Grab grid values of both components
    old_layer <- which(grepl(component, old$layout$name))
    new_layer <- which(grepl(component, new$layout$name))
    old_grid <- lapply(c('l'='l', 'r'='r', 't'='t', 'b'='b'),
                       function(x) old$layout[old_layer, x])
    new_grid <- lapply(c('l'='l', 'r'='r', 't'='t', 'b'='b'),
                       function(x) new$layout[new_layer, x])

    if (is.null(t)) {
        t <- max(old_grid$t)
    }
    if (is.null(b)) {
        b <- min(old_grid$b)
    }
    if (is.null(l)) {
        l <- min(old_grid$l)
    }
    if (is.null(r)) {
        r <- max(old_grid$r)
    }

    # Also set height/width of old component to accomodate new component to new
    if (set_height) {
        old$heights[t] <- new$heights[max(new_grid$t)]
    }
    if (set_width) {
        old$widths[r] <- new$widths[max(new_grid$r)]
    }

    gtable::gtable_add_grob(old, new$grobs[[new_layer]],
                            t=t, b=b, l=l, r=r)
}

# Creates top bar for combined mutation + qval plot, showing belongingness
# to an additional group
# var: string representing column in data frame to represent a group belonging bar
# df: data frame in long format, must have ID and var in it
# colour_schemes: List with colour schemes, should be an entry under 'var' name if provided
#    else default to ggplot2 hue based method. See heatmap_mutation_extended for precise
#    specification.
create_top_bar <- function(var, df, colour_schemes, legend_title_size=10, legend_text_size=8, legend_margin = -25) {
    ngroups <- length(unique(df[[var]]))
    if (is.null(colour_schemes[[var]])) {
        cols <- gg_color_hue(ngroups)
    } else {
        cols <- colour_schemes[[var]]
    }

    plt <- df %>%
                dplyr::rename(bar = var) %>%
                dplyr::distinct(ID, bar, cluster) %>%
            ggplot2::ggplot(ggplot2::aes(x=ID, y=1, fill=bar)) +
                ggplot2::geom_tile() +
                ggplot2::theme_bw() +
                ggplot2::scale_fill_manual(var, values=cols) +
                ggplot2::facet_grid(cols=ggplot2::vars(cluster), space='free_x', scales='free_x') +
                ggplot2::theme(legend.position="right",
                               legend.direction="horizontal",
                               legend.justification="left",
                               panel.grid = ggplot2::element_blank(),
                               panel.border = ggplot2::element_blank(),
                               panel.spacing.x = ggplot2::unit(0, "lines"),
                               axis.title = ggplot2::element_blank(),
                               axis.text = ggplot2::element_blank(),
                               axis.ticks = ggplot2::element_blank(),
                               legend.text = ggplot2::element_text(size=legend_text_size),
                               legend.title = ggplot2::element_text(size=legend_title_size),
                               legend.key.size = ggplot2::unit("0.02", "npc"),
                               legend.margin = ggplot2::margin(l=legend_margin)) +
                transparent_background
    ggplot2::ggplotGrob(plt)
}
