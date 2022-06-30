# This produces a Volcano plot with statistically significant p < 0.01 genes
# colored red (up-regulated) and blue (down-regulated). The top 20 most
# statistically significant and top 20 most highly regulated genes
# with an absolute average log2 fold change above 0.25 are labeled with
# their gene symbol.

# Input data is a dataframe from the Seurat FindMarkers output (with extra
# column for regulation). Meant to be iterated over the identity variable,
# which will subset the dataframe from the cluster column. Works best if the
# DGE is done without a threshold, logfc.threshold = 0, then non-
# significantly regulated genes, yet still expressed, appear as grey dots.

VolcanoPlot <- function(df, identity, title) {

  # subset each cluster, calc -log10(p-adj)
  deg <- df[df$cluster == identity, ]
  deg$neglog10p <- -(log10(deg$p_val_adj))

  # filter significantly upregulated genes
  deg_up <- deg[deg$regulation == "Up", ]
  deg_up_sig <- deg_up[deg_up$p_val_adj < 0.01, ]
  deg_up_sig_fc <- deg_up_sig[deg_up_sig$avg_log2FC > 0.25, ]
  deg_up_sig_fc_top20sig <- deg_up_sig_fc %>% slice_min(p_val_adj, n = 20)
  deg_up_sig_fc_top20fc <- deg_up_sig_fc %>% slice_max(avg_log2FC, n = 20)
  deg_up_tops <- rbind(deg_up_sig_fc_top20sig, deg_up_sig_fc_top20fc)
  deg_up_tops <- distinct(deg_up_tops)

  # filter significantly downregulated genes
  deg_down <- deg[deg$regulation == "Down", ]
  deg_down_sig <- deg_down[deg_down$p_val_adj < 0.01, ]
  deg_down_sig_fc <- deg_down_sig[deg_down_sig$avg_log2FC < -0.25, ]
  deg_down_sig_fc_top20sig <- deg_down_sig_fc %>% slice_min(p_val_adj, n = 20)
  deg_down_sig_fc_top20fc <- deg_down_sig_fc %>% slice_min(avg_log2FC, n = 20)
  deg_down_tops <- rbind(deg_down_sig_fc_top20sig, deg_down_sig_fc_top20fc)
  deg_down_tops <- distinct(deg_down_tops)

  # volcano plot
  vp <- ggplot(
    deg,
    aes(
      x = avg_log2FC,
      y = neglog10p
    )
  )
  vc <- vp + geom_point(colour = "grey") +
    geom_point(
      data = deg_up_sig,
      colour = "#eb2a0e"
    ) +
    geom_point(
      data = deg_down_sig,
      colour = "#2664ad"
    ) +
    geom_text_repel(
      data = deg_up_tops,
      aes(label = gene),
      colour = "#eb2a0e",
      fontface = "italic",
      size = 3,
      force = 2,
      nudge_x = max(abs(deg$avg_log2FC)),
      nudge_y = median(abs(deg$neglog10p)),
      direction = "both",
      segment.size = 0.1,
      segment.alpha = 0.3,
      segment.curvature = 0,
      max.overlaps = Inf
    ) +
    geom_text_repel(
      data = deg_down_tops,
      aes(label = gene),
      colour = "#2664ad",
      fontface = "italic",
      size = 3,
      force = 2,
      nudge_x = -(max(abs(deg$avg_log2FC))),
      nudge_y = median(abs(deg$neglog10p)),
      direction = "both",
      segment.size = 0.1,
      segment.alpha = 0.3,
      segment.curvature = 0,
      max.overlaps = Inf
    ) +
    xlab(expression("avg log"[2] * "(FC)")) +
    # ylab("-log10 (p-adj)") +
    ylab(expression("-log"[10] * "(p-adj)")) +
    xlim(
      -max(abs(deg$avg_log2FC)),
      max(abs(deg$avg_log2FC))
    ) +
    geom_hline(
      yintercept = -log10(0.01),
      linetype = "dashed",
      color = "grey",
      size = 0.1
    ) +
    # geom_vline(xintercept = -0.25,
    #            linetype="dashed",
    #            color = "grey",
    #            size = 0.1) +
    # geom_vline(xintercept = 0.25,
    #            linetype="dashed",
    #            color = "grey",
    #            size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14)
    ) +
    ggtitle(title)
  print(vc)
}
