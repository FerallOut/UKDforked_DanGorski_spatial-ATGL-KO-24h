ORA_down <- function (df, identity, detected_genes, title) {
  
  # rename variable
  i <- identity
  
  # subset each cluster, calc -log10(p-adj)
  deg <- df[df$cluster == i,]
  deg$neglog10p <- -(log10(deg$p_val_adj))
  
  # filter significant
  deg_up <- deg[deg$regulation == "Up",]
  deg_up_sig <- deg_up[deg_up$p_val_adj < 0.01,]
  deg_up_sig_fc <- deg_up_sig[deg_up_sig$avg_log2FC > 0.25,]
  deg_down <- deg[deg$regulation == "Down",]
  deg_down_sig <- deg_down[deg_down$p_val_adj < 0.01,]
  deg_down_sig_fc <- deg_down_sig[deg_down_sig$avg_log2FC < -0.25,]
  
  # ORA downregulated
  if  (length(deg_down_sig_fc$gene) > 5) {
    ORA_down <- enrichGO(gene =  deg_down_sig_fc$gene,
                         universe = detected_genes,
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         keyType = "SYMBOL",
                         pAdjustMethod = "bonferroni", 
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = F)
    
    # breakdown top GOBP terms
    ORA_down_sub <- ORA_down
    row.names(ORA_down_sub@result) <- NULL
    top_GO_down <- head(ORA_down_sub@result[,c("Description", "geneID")], n=5)
    colnames(top_GO_down) <- c("GO term", "Genes") 

    # plot ORA downregulated
  if (min(ORA_down@result$p.adjust) < 0.05){
    ora_down <- dotplot(ORA_down,
                        showCategory = 15,
                        label_format = 30,
                        orderBy = "x") +
      labs(title = title,
           subtitle = paste0("Niche: ", i)) +
      theme(legend.position = "right",
            legend.box = "vertical")
      print(ora_down)
      
      # create table of top GOBP terms with genes
      top_GO_down %>% 
        kbl() %>% 
        kable_styling(bootstrap_options = c("striped", 
                                            "hover", 
                                            "condensed", 
                                            "responsive"),
                      font_size = 12)
   }
  } else {
    print(ggplot() +
            theme_void() +
            geom_text(aes(0,0, label = "No significant downregulation")))
  }
}
