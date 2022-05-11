#!/bin/Rscript

library(shiny)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(tidyr)
library(ggtree)
library(aplot)
library(ggfortify)
library(ggplotify)
library(easyGgplot2)
library(ComplexHeatmap)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)

######defining functions
box2.func <- function(df, df2) {
  plot <- data.frame()
  for(i in 1:length(colnames(df))){
    expression <- df[,i]
    gene <- rownames(df)
    plottmp <- data.frame(expression=expression, 
                          gene=gene, 
                          sample=c(colnames(df)[i]), 
                          inducer=unlist(strsplit(colnames(df)[i], "_"))[3], 
                          container=unlist(strsplit(colnames(df)[i], "_"))[1],
                          type = c("all_genes_in_cells"))
    plot <- rbind(plot, plottmp)
  }
  plot2 <- data.frame()
  for(i in 1:length(colnames(df2))){
    expression <- df2[,i]
    gene <- rownames(df2)
    plottmp <- data.frame(expression=expression, 
                          gene=gene, 
                          sample=c(colnames(df2)[i]), 
                          inducer=unlist(strsplit(colnames(df2)[i], "_"))[3], 
                          container=unlist(strsplit(colnames(df2)[i], "_"))[1],
                          type = c("gene_on_circuits"))
    plot2 <- rbind(plot2, plottmp)
  }
  plot <- rbind(plot, plot2)
  return(plot)
}
pca2d.func <- function(mat) {
  pcamat <- prcomp(t(na.omit(mat)), scale. = T)
  pctmp <- as.data.frame(pcamat$x[,1:2])
  pctmp$inducer <- matchtreatment$rep[match(rownames(pctmp), matchtreatment$V2)]
  pctmp$container <- matchtreatment$container[match(rownames(pctmp), matchtreatment$V2)]
  ggplot(pctmp, aes(x=PC1, y=PC2, color=inducer, shape=container)) +
    geom_point(size = 5) + 
    stat_ellipse()
}
heat.func <- function(df) {
  plot <- data.frame()
  for(i in 1:length(colnames(df))){
    expression <- df[,i]
    gene <- rownames(df)
    plottmp <- data.frame(expression=expression, 
                          gene=gene, 
                          sample=c(colnames(df)[i]), 
                          inducer=unlist(strsplit(colnames(df)[i], "_"))[3], 
                          container=unlist(strsplit(colnames(df)[i], "_"))[1],
                          type = c("all"))
    plot <- rbind(plot, plottmp)
  }
  phc <- hclust(dist(as.matrix(t(df)))) %>% ggtree() + layout_dendrogram()
  
  pp <- ggplot(plot, aes(x=sample, y=gene, fill=expression)) + geom_tile() +
    theme_minimal()+
    scale_fill_viridis_c() +
    scale_y_discrete(position="right")+
    xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_text(angle = 90))
  
  p <- pp %>% insert_top(phc, height=0.1) 
  
  as.ggplot(p)
}
comp.heat.func <- function(df, annotab) {
  colanno <- HeatmapAnnotation(inducer = annotab$rep)
  ba <- columnAnnotation(dist = anno_boxplot(as.matrix(df)), height = unit(2, "cm"))
  ha <- rowAnnotation(dist = anno_density(as.matrix(df)), width = unit(2, "cm"))
  p <- Heatmap(df, 
               column_split=annotab$container, 
               border = TRUE, 
               top_annotation = colanno, 
               right_annotation = ha, 
               bottom_annotation = ba,
               width = unit(10, "cm"), height = unit(10, "cm"))
  as.ggplot(p)
}
level2.func <- function(df) {
  plot <- data.frame(row.names = rownames(df))
  for(i in 1:length(colnames(df))){
    expression <- df[,i]
    gene <- rownames(df)
    tmpvec <- df[,i]
    tmp <- df[,i]
    tmp[tmp < mean(tmpvec) - 0.5 * sd(tmpvec)] <- -1
    tmp[tmp <= mean(tmpvec) + 0.5 * sd(tmpvec) & tmp >= mean(tmpvec) - 0.5 * sd(tmpvec)] <- 0
    tmp[tmp > mean(tmpvec) + 0.5 * sd(tmpvec)] <- 1
    plottmp <- data.frame(expression=tmp)
    colnames(plottmp) <- colnames(df)[i]
    plot <- cbind(plot, plottmp)
  }
  return(plot)
  
}
level.cir.func <- function(df, skew) {
  plot <- data.frame(row.names = rownames(df))
  for(i in 1:length(colnames(df))){
    expression <- df[,i]
    gene <- rownames(df)
    tmpvec <- df[,i]
    tmp <- df[,i]
    tmp[tmp < mean(tmpvec) - skew * sd(tmpvec)] <- -1
    tmp[tmp <= mean(tmpvec) + skew * sd(tmpvec) & tmp >= mean(tmpvec) - skew * sd(tmpvec)] <- 0
    tmp[tmp > mean(tmpvec) + skew * sd(tmpvec)] <- 1
    plottmp <- data.frame(expression=tmp)
    colnames(plottmp) <- colnames(df)[i]
    plot <- cbind(plot, plottmp)
  }
  return(plot)
  
}
level4.func <- function(df, df2, skew, skew2) {
  plot <- data.frame(row.names = rownames(df))
  for(i in 1:length(colnames(df))){
    expression <- df[,i]
    gene <- rownames(df)
    tmpvec <- df2[,match(colnames(df)[i], colnames(df2))]
    tmp <- df[,i]
    print(sort.default(tmpvec)[round(skew *length(tmpvec))])
    print(sort.default(tmpvec)[round(skew2 *length(tmpvec))])
    tmp[tmp < sort.default(tmpvec)[round(skew2 *length(tmpvec))]] <- -1
    tmp[tmp <= sort.default(tmpvec)[round(skew *length(tmpvec))] & tmp >= sort.default(tmpvec)[round(skew2 *length(tmpvec))] ] <- 0
    tmp[tmp > sort.default(tmpvec)[round(skew *length(tmpvec))] ] <- 1
    plottmp <- data.frame(expression=tmp)
    colnames(plottmp) <- colnames(df)[i]
    plot <- cbind(plot, plottmp)
  }
  return(plot)
  
}
findgenes.func <- function(list, df) {
  outtab <- data.frame()
  for(i in list) {
    outtmp <- df[grep(i, rownames(df)),]
    outtab <- rbind(outtab, outtmp)
  }
  outtab <- unique.data.frame(outtab)
  return(outtab)
}
input1.mani.func <- function(vec1, vec2) {
  if(vec1 == "none") {
    exptab <- assay(keepgenes)
  }
  if(vec1 == "log2") {
    exptab <- countreads
  }
  if(vec1 == "rld") {
    exptab <- assay(rld_keepgenes)
  }
  if(vec1 == "fpm") {
    exptab <- fpm_keepgenes
  }
  samp <- subset(matchtreatment, matchtreatment$container %in% vec2)$V2
  exptab <- exptab[,match(samp, colnames(exptab))]
  return(exptab)
}
input2.mani.func <- function(vec1, vec2, df1) {
  if(vec1 == "circuitgenes") {
    targettab <- df1[match(circuitgenes, rownames(df1)),]
  }
  if(vec1 == "sensors") {
    targettab <- df1[match(sensors, rownames(df1)),]
  }
  if(vec1 == "gates") {
    targettab <- df1[match(gates, rownames(df1)),]
  }
  if(vec1 == "selfdefined") {
    genelist <- trimws(unlist(strsplit(as.character(vec2), ",")))
    targettab <- findgenes.func(genelist, df1)
  }
  return(targettab)
}

######loading data
#load(".RData")
circuitgenes <- rownames(countreads[-grep("locus", rownames(countreads)),])
sensors <- c("LacI",  "TetR",  "KanR", "AraC" )
gates <- c("AmtR",  "LitR",  "BM3R1", "SrpR",  "PhlF",  "YFP")

######shiny parts


menu <- controlbarMenu(
  id = "controlbarMenu",
  controlbarItem(
    "Tab 1",
    "Welcome to tab 1"
  ),
  controlbarItem(
    "Tab 2",
    checkboxGroupInput("selectk", 
                       h4("please select sample types from:"), 
                       choices = unique(matchtreatment$container),
                       selected = unique(matchtreatment$container)
    ),
    selectInput("norm", 
                h4("normalisation and transformation methods"), 
                choices =  c("none", 
                             "log2", 
                             "rld", 
                             "fpm"),
                selected = "log2"
    ),
    selectInput("part", 
                h4("circuit parts to visualise"), 
                choices =  c("circuitgenes", 
                             "sensors", 
                             "gates", 
                             "selfdefined"),
                selected = "circuitgenes"
    ),
    textInput("text", h4("type in other genes or gene families to visualise"), 
              value = "ara,caiT")
  ),
  controlbarItem(
    "Tab 3",
    "Welcome to tab 3"
  ),
)


shinyApp(
  ui = dashboardPage(
    skin = "midnight",
    header = dashboardHeader(
      title = "Circuit Navigator",
      dropdownMenu(type = "messages",
                   messageItem(
                     from = "Developer",
                     message = "Qian Li at Zhejiang Lab",
                     time = "2022"
                     
                   ),
                   messageItem(
                     from = "Support",
                     message = "https://github.com/lynceuslq",
                     icon = icon("life-ring")
                   )
                   
      )
    ),
    body = dashboardBody(
      tabItems(
        tabItem(tabName = "menu_1", "Content 1"),
        tabItem(
          tabName = "menu_2",
          fluidRow(
            box(title = "Dataset information",
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "tabset1", height = "250px",
                dataTableOutput('tab1') 
            ), 
            box(
              title = "PCA plot", 
              solidHeader = TRUE, 
              height = "500px",
              collapsible = TRUE,
              plotOutput("plot2", height ='500px') 
            )
            
          ),
          fluidRow(
            box(
              title = "distribution plot", 
              closable = TRUE, 
              width = 12,
              height = "800px",
              solidHeader = FALSE, 
              collapsible = TRUE,
              actionButton("update", "plot types"),
              sidebar = boxSidebar(
                id = "mycardsidebar",
                width = 30,
                selectInput("plott", 
                            h4("change distribution plot types"), 
                            choices =  c("box", "violin"),
                            selected = "violin"
                )
              ),
              plotOutput("distPlot")
            )
          ),
          
          fluidRow(
            box(
              title = "plot target gene expression", 
              closable = TRUE, 
              width = 12,
              height = "1000px",
              solidHeader = FALSE, 
              collapsible = TRUE,
              actionButton("update2", "methods for categorisation"),
              sidebar = boxSidebar(
                id = "mycardsidebar2",
                width = 30,
                selectInput("lev", 
                            h4("ways to categorise gene expression"), 
                            choices =  c("within circuits", 
                                         "comparison against all genes", 
                                         "comparison against house-keeping genes"),
                            selected = "within circuits"
                ),
                sliderInput("csk", 
                            "gene expression beyond the folds of standard deviation from mean:", 
                            value = 1, min = 0, max = 2, step= 0.05),
                
                sliderInput("upbound", 
                            "gene expression beyond the fraction of all genes:", 
                            value = 0.95, min = 0.8, max = 1),
                
                sliderInput("downbound", 
                            "gene expression below the fraction of all genes:", 
                            value = 0.9, min = 0.5, max = 1)
              ),
              #plotOutput("heatPlot")
              div(style='height:500px;overflow-y: scroll;',
                  uiOutput("plot"))
              
            )
            
          )
        ),
        tabItem(tabName = "menu_3", "Content 3")
      )

      
    ),
    
    sidebar = dashboardSidebar(
      minified = FALSE, 
      collapsed = FALSE,
      sidebarMenu(
        id = "sidebarMenu",
        lapply(1:3, function(i) {
          menuItem(
            sprintf("Page %s", i), 
            tabName = sprintf("menu_%s", i), 
            icon = icon("circle")
          )
        })
      )
    ),
    
    controlbar = dashboardControlbar(
      id = "controlbar",
      menu
    ),
    
    
  ),
  server = function(input, output, session) {
    
    observeEvent(input$sidebarMenu, {
      idx <- strsplit(input$sidebarMenu, "_")[[1]][2]
      if (idx == 2) {
        updateControlbar("controlbar")
      }
      updateControlbarMenu("controlbarMenu", selected = idx)
    })
    
    observeEvent(input$controlbarMenu, {
      if (input$controlbarMenu == "Tab 2") updateBoxSidebar("boxSidebar")
    })
    
    
    output$tab1 <- renderDataTable({
      matchtreatment
      
    },options = list(pageLength = 10))
    
    output$plot2 <- renderPlot({
      dat <- input1.mani.func(input$norm, input$selectk)
      show <- input2.mani.func(input$part, input$text, dat)
      p <- pca2d.func(dat)
      p
    })
    
    observe(print(input$mycardsidebar))
    
    output$distPlot <- renderPlot({
      dat <- input1.mani.func(input$norm, input$selectk)
      show <- input2.mani.func(input$part, input$text, dat)
      if(input$plott == "violin") {
        p <-ggplot(box2.func(dat, show), aes(x=sample, y=expression, fill=type)) +
          geom_violin(position=position_dodge(1)) +
          #geom_boxplot(width=0.1) +
          #facet_wrap(~type) + 
          theme(axis.text.x = element_text(angle = 90))
      }
      if(input$plott == "box") {
        p <-ggplot(box2.func(dat, show), aes(x=sample, y=expression, fill=type)) +
          geom_boxplot(position=position_dodge(1)) +
          #geom_boxplot(width=0.1) +
          #facet_wrap(~type) + 
          theme(axis.text.x = element_text(angle = 90))
        
      }
      p
    })
    
    observeEvent(input$update, {
      updateBoxSidebar("mycardsidebar")
    })
    
    
    observe(print(input$mycardsidebar2))
    
    output$plot <- renderUI({
      output$plot3 <- renderPlot({
        dat <- input1.mani.func(input$norm, input$selectk)
        show <- input2.mani.func(input$part, input$text, dat)
        p1 <- comp.heat.func(show, matchtreatment) + theme(text = element_text(size = 20)) 
        
        if(input$lev == "within circuits") {
          p2 <- comp.heat.func(level.cir.func(show, as.numeric(input$csk)), matchtreatment) + theme(text = element_text(size = 20)) 
        }
        
        if(input$lev == "comparison against all genes") {
          p2 <- comp.heat.func(level4.func(show, 
                                           dat, 
                                           as.numeric(input$upbound), 
                                           as.numeric(input$downbound)), 
                               matchtreatment
          ) + theme(text = element_text(size = 20)) 
        }
        
        p <- p1 + p2
        
        p
      })
      plotOutput('plot3', height ='800px')  
    })
    
    
    observeEvent(input$update2, {
      updateBoxSidebar("mycardsidebar2")
    })
    
  }
)

