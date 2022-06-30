# load packages
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(Seurat)){
  install.packages("Seurat") # v4.0.1
  library(Seurat)
}
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(patchwork)){
  install.packages("patchwork")
  library(patchwork)
}

# load functions and data
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/hearts.Rdata")
deg <- read.csv("results/differential-gene-expression/deg_clusters.csv")

# ui
ui <- fluidPage(
  verticalLayout(
    titlePanel("Bottermann et al. spatial-Atgl-KO-24h"),
    wellPanel(textInput("gene",
                        label = "Gene (e.g. Col1a1)"),
              splitLayout(cellWidths = c("50%", "50%"), 
                          plotOutput("Dim"),
                          plotOutput("Feature_split")),
              splitLayout(cellWidths = c("50%", "50%"), 
                          plotOutput("Violin_split"),
                          dataTableOutput("table"))),
    helpText("daniel.gorski@uni-duesseldorf.de")))

# server
server <- function(input, output) {
  output$Dim <- renderPlot({
    a <- SpatialDimPlot(hearts,
                        images = "Control",
                        pt.size.factor = 1.6
    ) +
      labs(fill = "Cluster") +
      ggtitle("Control") +
      theme(
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        title = element_text(size = 16)
      ) +
      guides(fill = guide_legend(override.aes = list(size = 5))) + NoLegend()
    b <- SpatialDimPlot(hearts,
                        images = "KO",
                        pt.size.factor = 1.6
    ) +
      labs(fill = "Cluster") +
      ggtitle("KO") +
      theme(
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3, "point"),
        title = element_text(size = 16)
      ) +
      guides(fill = guide_legend(override.aes = list(size = 5)))
    print(a + b + plot_layout(guides = "collect"))
  })
  output$Feature_split <- renderPlot({
    SpatialFeaturePlotScaled(
      object = hearts,
      group = "genotype",
      group.1 = "Control",
      group.2 = "KO",
      feature_of_interest = input$gene,
      from.meta.data = FALSE,
      group.1.title = "Control",
      group.2.title = "KO")
  })
  output$Violin_split <- renderPlot({
    p <- VlnPlot(hearts,
                 features = input$gene,
                 split.by = "genotype",
                 assay = "Spatial",
                 split.plot = TRUE
    ) +
      theme(
        plot.title = element_text(face = "italic", size = 20),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "right"
      )
    print(p)
  })
  output$table <- renderDataTable(
    deg %>% filter(gene == input$gene) %>% select(gene,
                                                  cluster,
                                                  regulation,
                                                  avg_log2FC,
                                                  pct.1,
                                                  pct.2,
                                                  p_val_adj),
    options = list(paging = FALSE, searching = FALSE)
  )
}

# app----
runApp(list(ui = ui, server = server), launch.browser = TRUE)
