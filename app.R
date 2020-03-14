# This is a Shiny web application. 
# http://shiny.rstudio.com/

# packages used in this app
library(shiny)
library(ggplot2)
library(gplots)
library(DESeq2)
library(RColorBrewer)
library(shinythemes)
library(pheatmap)
library(reshape2)
#source("DT")
library(DT)
source("mydds.R")
source("cmcdistance.R")


# Define UI for application 
ui <- fluidPage(
  
  theme = shinytheme("cerulean"),
  # Application title
  titlePanel("Shiny-DEG"),
  helpText(h4("a web application to analyze and visualize differentially expressed genes in RNA-seq")),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_file_type",h3("Use example data or upload your own data"),c("Example Data"="examplecounts","Upload Data"="upload"),selected = "examplecounts"),
      
      # Input: Select a file ----
      conditionalPanel(condition = "input.data_file_type=='upload'",
                       radioButtons("file_type",h4(em("please choose data type")),c("MATRIX"="matrix",'CSV'="csv"),selected = "matrix"),
                       
                       # Input: Select a file ----
                       conditionalPanel(condition = "input.file_type=='matrix'",
                                        fileInput("file1", "Input matrix Data",
                                                  multiple = TRUE)
                       ),
                       conditionalPanel(condition = "input.file_type=='csv'",
                                        fileInput("file3", "Input csv Data",
                                                  multiple = TRUE,
                                                  accept = c("text/csv","text/comma-separated-values,text/plain",".csv")
                                        )
                       )
      ),
      br(),
      actionButton("goButton", "Submit !"),
      tags$hr(),
      h3("DEG Analysis"),
      # Input: Select a Group ----
      #choose condition and replicates
      #choose experimental design, default replicates are 3, so we do not need replicates for now
      selectInput("Group_factor",h5(em("Please select experimental design")),list("Single factor"="type1","Multi-factor"="type2"),selected = "type1"),
      conditionalPanel(condition = "input.Group_factor=='type2'",
                       selectInput("multi_factor","Please select multi-factor design",list("factor1","factor2"),selected = "factor1")),
                       
      sliderInput("slider1", h5("FDR"),min = 0.0000000001, max = 0.05, value = 0.05),
      sliderInput("slider2", h5("log2Foldchange"),min = 0.5, max = 5.0, value = 2,step = 0.5),
      
      # Horizontal line ----
      tags$hr(),
      h3("DEG Visualization"),
      h4("Heatmap Figure"),
      #Z-score Choice
      selectInput("score", label=("Z-score Choice"),list("by matrix","by column")),
      #Distance Choice
      selectInput("distance", label=("Distance Choice"),list("Euclidean"="euclidean","Pearson correlation distance"="correlation","Manhattan"="manhattan")),
      #Method Choice
      selectInput("method", label=("Method Choice"),list("Average"="average","Complete"="complete","Median"="median","Single"="single","Centroid"="centroid","Ward.D"="ward.D2")),
      #Heatmap title
      textInput("text", "Figure Title",value = "DEG"),
      #dispaly gene name
      checkboxInput("checkbox1", "Show Gene Name", value = FALSE),
      #dispaly gene cluster
      checkboxInput("checkbox2", "Show Gene Cluster", value = FALSE),
      #heatmap color
      selectInput("color",label = ("Color Choice"),list("RdYlBu","NvWiFr")),
      
      # Horizontal line ----
      tags$hr(),
      h4("PCA Figure"),
      #PCA title
      textInput("text2", "Figure Title", 
                value = "PCA"),
      # dispaly PCA legend
      checkboxInput("checkbox4", "Show Legend", value = T),
      
      tags$hr(),
      h3("Download Tables and Figures"),
      #  download DEG Table
      h4("Download DEG Tables"),
      #Table Format Choice
      radioButtons("checkGroup2", "table Format Choice",
                   choices = list("csv" = 1, "txt" = 2),
                   selected = 1),
      downloadButton("downloadCsv", "Download DEG Table"),
      h4("Download Figures"),
      #Figure Format Choice
      radioButtons("checkGroup", "Figure Format Choice",
                   choices = list( "JPEG" = 2,"PDF"=3),
                   selected = 2),
      
      #Download Heatmap Figure
      downloadButton("downloadFigure", "Download Heatmap Figure"),
      
      #Download PCA Figure
      downloadButton("downloadFigure1", "Download PCA Figure"),
      
      #Download Boxplot Figure
      downloadButton("downloadFigure2", "Download Boxplot Figure"),
      
      #Download volcanoplot Figure
      downloadButton("downloadFigure3", "Download Volcano Plot Figure")
    ),
    # Show a plot of the generated distribution
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Instruction",
                           strong(h4(" Shiny-DEG is a web-based platform to help you analyze RNA-seq data and plot high quality figures.")),
                           br(),
                           em(h4("The Shiny-DEG allows users to visualize differentially expressed genes (DEG) starting with count data.")),
                           h4(em("Explore the app's features with the example data set pre-loaded.")),
                           em(h4("Upload your genes Expression data first,then submit your data.")),
                           br(),
                           strong(h4("Data Requirments")),
                           h5("1. Data must be uploaded as a matrix or CSV file"),
                           h5("2. File must be the raw counts,not normalized data,e.g.FPKM,TPKM,TPM"),
                           h5("3. File must have a header row."),
                           h5("4. First column must be gene identifiers."),
                           br(),
                           em(h4("Example Data format")),
                           h5("Each row denotes a gene, each column denotes a sample."),
                           em(h5("single factor data format:")),
                           img(src = "example3.png", height = 245, width = 700),
                            br(),
                           em(h5("multi-factor data format:")),
                           img(src = "example5.png", height = 210, width = 700),
                           br(),
                            h4(" The Shiny-DEG workflow."),
                           img(src = "Workflow.png", height = 515, width = 700),
                           tags$hr()
                  ),
                  
                  
                  tabPanel("InputData",
                           #tags$hr(),
                           br(),
                           dataTableOutput('countdataDT')),
                  tabPanel("DEG",
                           #tags$hr(),
                           br(),
                           dataTableOutput('countdataDT2')),
                  tabPanel("Boxplot", 
                           #tags$hr(),
                           br(),
                           plotOutput("boxplot"),
                           br(),
                           plotOutput("density"),
                           br(),
                           h4("Summary:"),
                           verbatimTextOutput("analysis1")
                           
                           ),
                  tabPanel("Volcano plot", 
                           #tags$hr(),
                           br(),
                           plotOutput("plot1", height = 600),
                           br(),br(),
                           plotOutput("scatterplot", height = 600,
                                      dblclick = "scatterplot_dblclick",
                                      brush = brushOpts(
                                        id = "scatterplot_brush",
                                        resetOnNew = TRUE
                                      )
                           )
                    ),
                  tabPanel("Heatmap", 
                           br(),
                           br(),
                           plotOutput("plot")
                           ),
                  tabPanel("PCA", 
                           br(),
                           plotOutput("plot2")),
                  tabPanel("Help", 
                           br(),
                           h4("Source code can be found on github:",a("https://github.com/344968067/shiny-DEG")),
                           br(),
                           ####specify each column
                           h4("DEG Table"),
                           h5("Column A provide gene name."),
                           h5("Column C and column G provide Fold Changes and FDR,respectively."),
                           h5(" We use both log2FC and FDR to filter DEG."),
                           img(src = "example4.png", height = 210, width = 600),
                           br(),
                           h4("DEG Visualization"),
                           img(src = "example.png", height = 300, width = 250),
                           img(src = "example2.png", height = 300, width = 225),
                           img(src = "PCA1.png", height = 250, width = 250),
                           img(src = "volcanoplot.png", height = 300, width = 300),
                           img(src = "Boxplot.png", height = 200, width = 400),
                           br(),
                          h4("For more questions or suggestions, please contact us: Dr. Sufang Wang,", a("email: sufangwang@nwpu.edu.cn")),
                           tags$hr()
                           )
                          )
                        )
                      )
                  )

# Define server
server <- function(input, output) {
  #Use Example file or upload your own data
  dataInput<-reactive({
    input$goButton
    isolate({
    validate(
      need((input$data_file_type=="examplecounts")|((!is.null(input$file1))|(!is.null(input$file3))),
           message = "Please select a file")
    )
    if(input$data_file_type=="examplecounts"){
      inFile<-read.delim("Data/example.matrix",
                         header = TRUE,
                         sep = "\t",
                         quote = "\t",dec = ".",
                         fill = TRUE,
                         stringsAsFactors = F,
                         row.names = 1)
    }else {
      if(input$file_type=="matrix"){
        req(input$file1)
        inFile<-read.delim(input$file1$datapath,
                           header = TRUE,
                           sep = "\t",
                           quote = "\t",dec = ".",### when upload matrix data saved by ourself ,but raw matrix,be careful of dec
                           fill = TRUE,
                           stringsAsFactors = F,
                           row.names = 1)
      }else{      
        req(input$file3)
        inFile <- read.csv(input$file3$datapath,
                           header = TRUE,
                           sep = ",",row.names = 1)
      }
    }
    return(inFile)
    })
  })
  #ensure condition by group choice
  datasetInputcondition <- reactive({
    input$goButton
    isolate({
    inFile <- dataInput()
    if (is.null(inFile))
      return(NULL)
    if(input$data_file_type=="examplecounts"){
      condition=c(rep("ctrl",3),rep("exp",3))
    }else{
      if (dim(inFile)[1] != 1) {
        if (input$Group_factor=="type1"){
           condition=c(rep("ctrl",3),rep("exp",3))
        }else{
          return(NULL)
        }
      }
     }
      return(condition)
    })
  })
  
  # screen gene and get DESeq result by our defined function (mydds)
  datasetInput <- reactive({
    inFile <- dataInput()
    condition<-datasetInputcondition()
    if (input$Group_factor=="type1"){
      mydds(inFile,condition)
    }else if(input$Group_factor=="type2"){ 
      mycountData <- round(inFile,digits=0)
      #multi-factor design
      factor1 = c( rep("a1",6),rep("a2",6))
      factor2=c(rep("c1",3),rep("c2",3),rep("c1",3),rep("c2",3))
      data1Design <- data.frame(row.names = colnames( mycountData ),factor1 =as.factor(factor1),factor2=as.factor(factor2))
      mycolData <- data1Design
      mydds <- DESeqDataSetFromMatrix(countData = mycountData, colData = mycolData,design = ~factor1+factor2)
      mydds <- mydds[ rowSums(counts(mydds)) > 100, ]
      mydds <- DESeq(mydds)
      if(input$multi_factor=="factor1"){
        myres2 <- results(mydds,contrast=c("factor1","a1","a2"))
      }else{
        myres2 <- results(mydds,contrast=c("factor2","c1","c2"))
      }
      return(myres2)
    }
  })
  
  #screen DEG data by self-define FDR and log2FC,and sort data
  datasetInput3 <- reactive({
    # self-define FDR and log2FC
    selfFDR <- input$slider1
    selflog2FC <- input$slider2
    resSig <- subset(datasetInput(), (datasetInput()$padj < selfFDR & abs(datasetInput()$log2FoldChange) >= selflog2FC))
    #sort the Data
    resSig<-resSig[order(resSig$log2FoldChange,decreasing = FALSE),]
    return(resSig)
  })
  output$analysis1<-renderPrint({
    alldata<- dataInput()
    colnames( alldata ) <- sub("(*).genes.results","\\1",colnames(alldata))
    summary(alldata)
  }) 
  
  #Display chosen data(upload,previous or example data) by tables
  output$countdataDT <- renderDataTable({
    alldata <-dataInput()
    colnames( alldata ) <- sub("(*).genes.results","\\1",colnames(alldata))
    alldata<-as.data.frame(alldata,keep.rownames=TRUE)
  })
  
  #Display whole DEG data by table
  output$countdataDT2 <- renderDataTable({
    if (is.null(datasetInput())){
      return(NULL)
    }
    #head(myres)
    resSig<- datasetInput3()
    resSig<-as.data.frame(resSig,keep.rownames=TRUE)
  })
  #Density plot
  output$density <- renderPlot({
    alldata <- dataInput()
    colnames( alldata ) <- sub("(*).genes.results","\\1",colnames(alldata))
    alldata<- melt(alldata,measure.vars = colnames(alldata),variable.name = "Sample",value.name = "Counts")
    alldata$Sample=as.factor(alldata$Sample)
    ggplot(alldata,aes(x=log10(Counts+1),fill=Sample))+geom_density()+ggtitle("Density plot")+theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))
  })
  #Boxplot
  output$boxplot <- renderPlot({
    alldata <- dataInput()
    colnames( alldata ) <- sub("(*).genes.results","\\1",colnames(alldata))
    #alldata<- melt(alldata,measure.vars = colnames(alldata),variable.name = "Sample",value.name = "Counts")
    #alldata$Sample=as.factor(alldata$Sample)
    #ggplot(alldata, aes(x=Sample, y=log10(Counts+1))) + geom_boxplot(aes(fill=Sample))
    boxplot(log10(alldata+1),col=rainbow(9),pch=20, main="Boxplot", cex=1.0, xlab="Group", ylab="log[10](Counts+1)")
   })
  
  #volcano plot--Single zoomable plot (on below)
  ranges <- reactiveValues(x = NULL, y = NULL)
  output$scatterplot <- renderPlot({
    myres<-datasetInput()
    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
    degTotal<-myres
    topT <- as.data.frame(degTotal)
    # self-define FDR and log2FC
    selfFDR <- input$slider1
    selflog2FC <- input$slider2
    #plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~Fold~Change), ylab=bquote(~-log[10]~FDR),xlim=c(-13,13))
    ggplot(topT, aes(log2FoldChange, -log10(padj),color=log10(topT$baseMean))) + scale_color_gradient(low="green", high="red")+ggtitle("Zoomable Volcano Plot")+theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
      geom_point() +xlim(-12,12)+geom_vline(xintercept=c(-selflog2FC,0,selflog2FC), linetype="dotted")+ geom_hline(aes(yintercept=-log10(max(topT$pvalue[topT$padj<selfFDR], na.rm=TRUE))),linetype="dashed")+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+labs(color="log10(baseMean)")
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$scatterplot_dblclick, {
    brush <- input$scatterplot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    }else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  #volcano plot
  output$plot1 <- renderPlot({
    myres<-datasetInput()
    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
    degTotal<-myres
    topT <- as.data.frame(degTotal)
    # self-define FDR and log2FC
    selfFDR <- input$slider1
    selflog2FC <- input$slider2
    #Adjusted P values (FDR Q values)
    with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~Fold~Change), ylab=bquote(~-log[10]~FDR),xlim=c(-11,11)))
    with(subset(topT, padj<=selfFDR & log2FoldChange>=selflog2FC), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
    with(subset(topT, padj<=selfFDR & log2FoldChange<=(-selflog2FC)), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))

    #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
    abline(v=0, col="black", lty=3, lwd=1.0)
    abline(v=-selflog2FC, col="black", lty=4, lwd=2.0)
    abline(v=selflog2FC, col="black", lty=4, lwd=2.0)
    abline(h=-log10(max(topT$pvalue[topT$padj<selfFDR], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
    legend("topright", legend=c("Up","Down","Normal"),title = "Significance",col=c("red","green","black"), pch=16, xpd=T, cex=0.5,horiz=F)
    
  })
  #correlation plot-sample
  output$correlation<-renderPlot({
    all_data<-dataInput()
    # new matrix colnames,which match raw data
    if(input$data_file_type=="examplecounts"){
      names<-colnames(dataInput())
    }else{
      if(input$data_file_type!="examplecounts"){
          names <-c("Ctrl1","Ctrl2","Ctrl3","Exp1","Exp2","Exp3")
       }
     }
    colnames(all_data)<-names
    matrix<-cor(all_data[1:dim(all_data)[2]])
    pheatmap(matrix, cellwidth = 60,treeheight_col = 50,treeheight_row = 50,display_numbers = T,number_color = "black",cellheight=36,main = "Sample Correlation Plot")
  })
  #heatmap
  output$plot <- renderPlot({
    resSig<-datasetInput3()
    deg<-row.names(resSig)
    sig <- matrix (0,nc=dim(dataInput())[2],nr=length(deg))
    
    # new matrix colnames,which match raw data
    if(input$data_file_type=="examplecounts"){
      names <-c("Ctrl1","Ctrl2","Ctrl3","Exp1","Exp2","Exp3")
    }else{
      if (input$Group_factor=="type1"){
        names<-colnames(dataInput())
      }else if (input$Group_factor=="type2"){
          names <-c("F1_L1","F1_L1","F1_L1","F1_L2","F1_L2","F1_L2","F2_L1","F2_L1","F2_L1","F2_L2","F2_L1","F2_L2")
        }
    }
    sig <- as.data.frame(sig)
    rownames(sig) <- deg
    colnames(sig) <- names
    for(i in 1:length(deg)){
      sig[i,] <- (dataInput())[which(row.names(dataInput()) == deg[i] ), ]
    }
    # log transform then normalize
    sig_matrix <- data.matrix(sig)
    sig_matrix <- log10(sig_matrix+1)
    file <- sig_matrix
    #normalized data either by column or by matrix
    if (input$score == "by column"){
      activity.mean <- apply(file, 2,mean, na.rm=T)
      activity.sd <- apply(file, 2,sd,na.rm=T)
      zscore.mat <- sweep(file, 2, activity.mean, "-")
      zscore.mat <- sweep(zscore.mat, 2, activity.sd, "/")
    } else if (input$score == "by matrix") {
      activity.mean <- mean(file, na.rm=T)
      activity.sd <- sd(file, na.rm=T)
      zscore.mat <- sweep(file, 1, activity.mean, "-")
      zscore.mat <- sweep(zscore.mat, 1, activity.sd, "/")
    }
    dist.choice <- input$distance
    if(input$color=="RdYlBu"){
    rc <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(deg))
    }else{
    rc<-colorRampPalette(c("navy", "white", "firebrick3"))(length(deg))
    }
    if (input$checkbox1==FALSE){
      labRow=NA
    } else{
      labRow=NULL  
    }
    if (input$checkbox2==FALSE){
      Rowv=NA
    }else{
      Rowv=NULL
    }
    #heatmap(zscore.mat,main=input$text,col=rc, Rowv = Rowv,RowSideColors = rc,labRow = labRow)
    pheatmap(zscore.mat,color =rc, treeheight_row = 60,treeheight_col = 60,main=input$text,cluster_rows = input$checkbox2,cluster_cols= T,show_rownames= input$checkbox1,show_colnames= T, border=F,clustering_distance_row=input$distance,clustering_method=input$method,scale="row", cellwidth = 60,fontsize = 10)
  })

  #PCA Figure
  output$plot2 <- renderPlot({
    dataT <- t(dataInput())
    dataT2 <-  dataT[, colSums(dataT != 0) > 0.1]
    dataT3 <- log10(dataT2+1)
    dataPCA <- prcomp(dataT3)
    dataPCA3 <- dataPCA
    pc1var <- (summary(dataPCA)$importance[2])*100
    pc1var <-paste0(as.character(pc1var),"%")
    pc1var <-paste("PC1 explains", pc1var,"variance")
    pc2var <- (summary(dataPCA)$importance[5])*100
    pc2var <-paste0(as.character(pc2var),"%")
    pc2var <-paste("PC2 explains", pc2var,"variance")
    if(input$data_file_type=="examplecounts"){
      refClass <- c(1,1,1,2,2,2)
    }else{
      if (input$Group_factor=="type1"){
        refClass <- c(1,1,1,2,2,2)
      }else if (input$Group_factor=="type2"){
        refClass <- c(1,1,1,2,2,2,3,3,3,4,4,4)      
      }
    }
    refClass <- factor(refClass)
    par(pin=c(3.6,4)) 
    if(input$Group_factor=="type1"){
      plot(dataPCA$x[,1:2],col =refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
      if(input$checkbox4==T){
        legend(35,9, legend=c("Ctrl","Exp"),col=c(1,2), pch=16, xpd=T, cex=0.7,horiz=F)
      }
    }else if(input$Group_factor=="type2"){
      plot(dataPCA$x[,1:2],col = refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
      if(input$checkbox4==T){
        legend(35,9, legend=c("F1-L1","F1-L2","F2-L1","F2-L2"),col=c(1,2,3,4), pch=16, xpd=T, cex=0.7,horiz=F)
      }
    }
  })
  output$instructionspdf <- downloadHandler(filename="Instructions.pdf",
                                            content=function(file){
                                              file.copy("instructions/Instructions.pdf",file)
                                            })
  
  # Downloadable  csv of DEG dataset ----
  output$downloadCsv <- downloadHandler(
    filename = function() {
      if(input$checkGroup2==1){
        paste("DEG", ".csv", sep = "")
      }else{
        paste("DEG", ".txt", sep = "")
      }
    },
    content = function(file) {
      resSig<-datasetInput3()
      #write DEG data to csv
      if(input$checkGroup2==1){
        write.csv(resSig, file)
      }else {
        write.table(resSig, file) 
      }
    })
  
  # Download heatmap Figure ----
  output$downloadFigure <- downloadHandler(
    filename = function(){
      if(input$checkGroup==2){
        paste('Heatmap.jpeg')
      }else{
        paste('Heatmap.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==2){
        jpeg(file,width=4000,height=2300,res=500)
      }else{
        pdf(file)
      }
      resSig<-datasetInput3()
      deg<-row.names(resSig)
      sig <- matrix (0,nc=dim(dataInput())[2],nr=length(deg))
      if(input$data_file_type=="examplecounts"){
        names <-c("Ctrl1","Ctrl2","Ctrl3","Exp1","Exp2","Exp3")

      }else{
        if (input$Group_factor=="type1"){
          names<-colnames(dataInput())          
        }else if (input$Group_factor=="type2"){
          names <-c("F1_L1","F1_L1","F1_L1","F1_L2","F1_L2","F1_L2","F2_L1","F2_L1","F2_L1","F2_L2","F2_L1","F2_L2")
        }
      }
      sig <- as.data.frame(sig)
      rownames(sig) <- deg
      colnames(sig) <- names
      for(i in 1:length(deg)){
        sig[i,] <- (dataInput())[which(row.names(dataInput()) == deg[i] ), ]
      }
      # log transform then normalize
      sig_matrix <- data.matrix(sig)
      sig_matrix <- log10(sig_matrix+1)
      file <- sig_matrix
      #normalized data either by column or by matrix
      if (input$score == "by column"){
        activity.mean <- apply(file, 2,mean, na.rm=T)
        activity.sd <- apply(file, 2,sd,na.rm=T)
        zscore.mat <- sweep(file, 2, activity.mean, "-")
        zscore.mat <- sweep(zscore.mat, 2, activity.sd, "/")
      } else if (input$score == "by matrix") {
        activity.mean <- mean(file, na.rm=T)
        activity.sd <- sd(file, na.rm=T)
        zscore.mat <- sweep(file, 1, activity.mean, "-")
        zscore.mat <- sweep(zscore.mat, 1, activity.sd, "/")
      }
      
      if(input$color=="RdYlBu"){
        rc <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(deg))
      }else{
        rc<-colorRampPalette(c("navy", "white", "firebrick3"))(length(deg))
      }
      if (input$checkbox1==FALSE)
      {
        labRow=NA
      }else{
        labRow=NULL
      }
      if (input$checkbox2==FALSE){
        Rowv=NA
      }else{
        Rowv=NULL
      }
      pheatmap(zscore.mat,color =rc, treeheight_row = 60,treeheight_col = 60,main=input$text,cluster_rows = input$checkbox2,cluster_cols= T,show_rownames= input$checkbox1,show_colnames= T, border=F,clustering_distance_row=input$distance,clustering_method=input$method,scale="row", cellwidth = 60,fontsize = 10)
      dev.off()
    })
  
  #download PCA Figure---
  output$downloadFigure1 <- downloadHandler(
    filename = function(){
      if(input$checkGroup==2){
        paste('PCA.jpeg')
      }else{
        paste('PCA.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      dataT <- t(dataInput())
      dataT2 <-  dataT[, colSums(dataT != 0) > 0.1]
      dataT3 <- log10(dataT2+1)
      dataPCA <- prcomp(dataT3)
      dataPCA3 <- dataPCA
      pc1var <- (summary(dataPCA)$importance[2])*100
      pc1var <-paste0(as.character(pc1var),"%")
      pc1var <-paste("PC1 explains", pc1var,"variance")
      pc2var <- (summary(dataPCA)$importance[5])*100
      pc2var <-paste0(as.character(pc2var),"%")
      pc2var <-paste("PC2 explains", pc2var,"variance")
      if(input$data_file_type=="examplecounts"){
        refClass <- c(1,1,1,2,2,2)
      }else{
        if (input$Group_factor=="type1"){
          refClass <- c(1,1,1,2,2,2)
        }else if (input$Group_factor=="type2"){
          refClass <- c(1,1,1,2,2,2,3,3,3,4,4,4)      
        }
      }
      refClass <- factor(refClass)
      par(pin=c(4.5,4.2))
      if(input$Group_factor=="type1"){
        plot(dataPCA$x[,1:2],col =refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
        if(input$checkbox4==T){
          legend(35,9, legend=c("Ctrl","Exp"),col=c(1,2), pch=16, xpd=T, cex=0.7,horiz=F)
        }
      }else if(input$Group_factor=="type2"){
        plot(dataPCA$x[,1:2],col = refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
        if(input$checkbox4==T){
          legend(35,9, legend=c("F1-L1","F1-L2","F2-L1","F2-L2"),col=c(1,2,3,4), pch=16, xpd=T, cex=0.7,horiz=F)
        }
      }
      
      dev.off() 
    })
  #download Boxplot Figure---
  output$downloadFigure2 <- downloadHandler(
    filename = function(){
      if(input$checkGroup==2){
        paste('Boxplot.jpeg')
      }else{
        paste('Boxplot.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      alldata <- dataInput()
      boxplot(log10(alldata+1),col=rainbow(9),pch=20, main="Boxplot", cex=1.0, xlab="Group", ylab="log[10](Counts+1)")
      dev.off() 
    })
  
  #download Volcano Figure---
  output$downloadFigure3 <- downloadHandler(
    filename = function(){
      if(input$checkGroup==2){
        paste('Volcano.jpeg')
      }else{
        paste('Volcano.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      myres<-datasetInput()
      par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
      degTotal<-myres
      topT <- as.data.frame(degTotal)
      # self-define FDR and log2FC
      selfFDR <- input$slider1
      selflog2FC <- input$slider2
      #Adjusted P values (FDR Q values)
      with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~Fold~Change), ylab=bquote(~-log[10]~FDR),xlim=c(-11,11)))
      with(subset(topT, padj<=selfFDR & log2FoldChange>=selflog2FC), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
      with(subset(topT, padj<=selfFDR & log2FoldChange<=(-selflog2FC)), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))
      #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
      abline(v=0, col="black", lty=3, lwd=1.0)
      abline(v=-selflog2FC, col="black", lty=4, lwd=2.0)
      abline(v=selflog2FC, col="black", lty=4, lwd=2.0)
      abline(h=-log10(max(topT$pvalue[topT$padj<selfFDR], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
      legend("topright", legend=c("Up","Down","Normal"),title = "Significance",col=c("red","green","black"), pch=16, xpd=T, cex=0.5,horiz=F)
      dev.off() 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

