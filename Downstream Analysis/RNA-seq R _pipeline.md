# RNA-seq Analysis Using R: DESeq2 Workflow
***
>#### This notebook demonstrates a step-by-step workflow for analyzing RNA-seq data using DESeq2. The analysis includes preprocessing the expression matrix, handling metadata, identifying differentially expressed genes (DEGs), and Plotting important figures such as volcano plots, Heatmap, and PCA.
## The Workflow goes as follow:
![Downstream v2](https://github.com/user-attachments/assets/95bff580-b610-44a5-9b61-c0ce52823371)
***
## 1- Loading libraries
```{r}
 setwd("E:/1.Fresh Grad/02_EgComBio2023/NGS RNA-seq")
    library(readr)
    library(org.Hs.eg.db)
    library(dplyr)
    library(ggplot2)
    library(ggfortify)             
    library(DESeq2)
    library(gridExtra)
    library(ggrepel)
    library(ComplexHeatmap)
    library(EnsDb.Hsapiens.v86)
    library(GOplot)
```
***
## 2- Loading the data and metadata
```{r}
    merged_counts <- read_delim("merged_counts.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    merged_counts <- merged_counts[,1:29]
    Metadata <- read_csv("Metadata.csv")
```
***

## 3- Gene Annotation
#### Convert gene IDs to gene symbols using EnsDb.Hsapiens.v86
 ```{r}
   Ens_genename <- as.character(merged_counts$Transcript_ID)
    # Remove version numbers (e.g., ".1") from Ensembl transcript IDs
    Ens_genename <- sub("\\..*", "", Ens_genename)
    gene_symbol <- mapIds(
      EnsDb.Hsapiens.v86,
      keys = Ens_genename,
      keytype = "TXID",  # Ensembl transcript ID
      column = "SYMBOL", # Gene symbol
      multiVals = "first"
    )
    # Create a data frame with gene IDs and corresponding gene symbols
    gene_df <- data.frame(
        Transcript_ID = names(gene_symbol),
        Gene_symbol = as.vector(gene_symbol),
        stringsAsFactors = FALSE)
    gene_df <- na.omit(gene_df)             #remove NA genes (not found)
    mtx <- as.matrix(merged_counts)
    mtx[, 1] <- gsub("\\..*", "", mtx[, 1])             #remove ensembl gene version
    
    
    # Merge the gene symbols with the raw counts data
    mtx <- merge(mtx, gene_df, by="Transcript_ID", all.x=TRUE)
    mtx <- mtx %>% dplyr::select(Gene_symbol, everything(), -Transcript_ID)
```
***

## 4- Preprocessing
###    (a) metadata
```{r Preprocesing I}
    meta <- Metadata %>%
        dplyr::select(Run,sample_Name, source_name) %>%
        rename(sample_name = sample_Name,  Condition = source_name)
    
    meta <- meta %>%
        mutate(Condition = case_when(
            Condition == "Triple-Negative Breast Cancer" ~ "TNBC",
            Condition == "normal breast epithelial" ~ "Control",
            TRUE ~ Condition))
  


```
###    (b) Preprocess expression matrix
#### Clean the data (NAs, Duplicates, Zero variance genes)
#### A- Handling Duplicates
```{r}
    mtx.agg <- mtx %>% dplyr::select(everything(),-Gene_symbol)
    mtx.agg <-apply(mtx.agg, 2, as.integer)
    mtx.agg = aggregate(mtx.agg, by = list(mtx$Gene_symbol), FUN = mean)
    
    gene_names=unlist(mtx.agg$Group.1)
    mtx.agg <- mtx.agg[,-1]
    mtx.agg <-apply(mtx.agg, 2, as.integer)
    rownames(mtx.agg) <- gene_names
```
#### B- Removing Zero variance genes
```{r}
   varrow <- apply(mtx.agg, 1, var, na.rm=TRUE)
    cons_var <- (varrow == 0 | is.na(varrow))
    dim(mtx.agg)
    mtx.agg.var <- mtx.agg[!cons_var,]
    dim(mtx.agg.var)
```
***
## 5- Exploratory Data Analysis: 
### A- Box plot and Histogram
```{r Exploration}
    boxplot(log2(merged_counts[,-1]+1), main="Exploratory Box plot", ylab="log2_count", las=2)
    boxplot(merged_counts[,-1], main="Exploratory Box plot", ylab="Counts", las=2)
    hist(as.matrix(log2(merged_counts[,-1]+1)), main="Exploratory Histogram", xlab="Sample", ylab="log2_count", breaks=50)

    par(mfrow=c(1,2))
    hist(mtx.agg.var, main ="Histogram before 1ry analysis")    # Preprocessed data histogram
    boxplot(mtx.agg.var, main = "boxplot of Raw data", col = seq(1:ncol(mtx.agg.var)), las = 2)
```
***
### B- PCA
```{r PCA}
    pca_result <- prcomp(t(mtx.agg.var), scale. = TRUE)
    pca_data <-pca_result$x
    pca_result$x
    autoplot(pca_result, data = meta, colour = 'Condition',frame = T,label = T, label.size = 3,shape="Condition",main="All Data PCA")+ theme_classic()+ 
    theme(plot.title =element_text(
        size = 20,           # Increase title size
        face = "bold",       # Make title bold
        color = "darkblue",  # Change title color
        hjust = 0.5        ))
```
***
## 6- Differential expression analysis (DEseq2)
#### Prepare the data and Filter out low-count genes
```{r DEseq2}

    ## remove the low counts data(recomended by deseq2 (asmaa))
    Exp <- mtx.agg.var
    dim(Exp)                                #20654    28
    Exp=Exp[which(rowSums(Exp) > 100),]
    dim(Exp)                                #11656    28
    # Ensure column names in expression data match metadata
    all(colnames(Exp) %in% meta$Run)
    # Ensure column names in expression data are in the same order as metadata
    all(colnames(exp) == meta$Run)
    
    

    
    
```
#### Create DESeqDataSet object and Run DESeq2 analysis
```{r Follow DEseq2}
    dds = DESeqDataSetFromMatrix( countData = Exp, colData = meta , design = ~ Condition)
    dds.run = DESeq(dds)
    ### direct results or specifying the contrast (to make a res object based on two specific conditions/treatment)
    res=results(dds.run, contrast = c("Condition","TNBC","Control") )
    # remove nulls (recomended by deseq2 (asmaa))
    res=res[complete.cases(res), ] 
    summary(res)
    plotMA(res)
```
#### Extract results
```{r}
#Extract DEMs
    res.df=as.data.frame(res)
    res.dems=res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange) >2,]
    res.dems=res.dems[order(res.dems$log2FoldChange), ]
    
    dems=rownames(res.dems)
    write.table(dems, file="dems_tnbc.txt", quote = F, col.names = T, row.names = F)
```
***
## 7- Data normalization & Scaling (VST)
```{r Normalization}
   #### get the normalized and logged transformed values of all exp data
    ntd=normTransform(dds)
    Exp.norm= assay(ntd)
    
    ### get the normalized expression levels of the dems ###################
    dems.exp=Exp.norm[dems,]

    ###  Scaling the Data ###################
    Exp.norm.scaled=scale(Exp.norm,center = TRUE, scale = TRUE)
```
***
### 8- Visualize Normalization
```{r visualize}
# Adjust margins to prevent truncation
    par(mfrow = c(1,2), mar = c(6, 5, 5, 2))  # (bottom, left, top, right)
# Boxplot after normalization
    boxplot(Exp.norm, col = "red",
            border = "black",
            horizontal = FALSE,
            notch = FALSE,
            main = "After Normalization",
            xlab = "Samples",
            ylab = "Expression Levels",
            cex.axis = 0.8,  # Adjust axis text size
            cex.lab = 1,     # Adjust label size
            cex.main = 1.2)  # Adjust title size
    
# Boxplot after normalization & scaling
    boxplot(Exp.norm.scaled, col = "violet",
            border = "black",
            horizontal = FALSE,
            notch = FALSE,
            main = "After Normalization & Scaling",
            xlab = "Samples",
            ylab = "Expression Levels",
            cex.axis = 0.8,
            cex.lab = 1,
            cex.main = 1.2)
    
# Reset graphical parameters after plotting
    par(mfrow = c(1,1))

    dfExp<- stack(as.data.frame(Exp.norm))
    dfExp2<- stack(as.data.frame(Exp.norm.scaled))
    plot1 <- ggplot(dfExp, aes(x = values, fill = ind)) +
        geom_density(alpha = 0.5) +
        theme_classic()+
        ggtitle("After Normalization")+
        theme(plot.title = element_text(
            size = 20,           # Increase title size
            face = "bold",       # Make title bold
            color = "darkblue",  # Change title color
            hjust = 0.5        )) # Center the title (0.5 for center)
    
    plot2 <- ggplot(dfExp2, aes(x = values, fill = ind)) +
        geom_density(alpha = 0.5) +
        ggtitle("After Normalization & Scaling", )+
        theme_classic()+
        theme(plot.title = element_text(
            size = 20,           # Increase title size
            face = "bold",       # Make title bold
            color = "darkblue",  # Change title color
            hjust = 0.5        )) # Center the title (0.5 for center)
    grid.arrange(plot1, plot2, ncol = 2)

```
***
## 9- Downstream plots

###    (b) PCA on DEGs
```{r PCA again}
     options(repr.plot.width=10,repr.plot.height=8)
    #par(mfrow=c(1,2))
    pca.2 <- prcomp(t((dems.exp)))
    autoplot(pca.2, data = meta, colour = 'Condition',frame = T,label = T, label.size = 4,shape="Condition",main="DEGs PCA")+ 
        stat_ellipse(aes(color = Condition), level = 0.95, type = "norm", linetype = 1, size = 1)+
        theme_classic() +
        theme(plot.title =element_text(size = 20, face = "bold",color = "darkblue",hjust = 0.5))
```
###    (b) Volcano plot
```{r Volcano plot}
res_df <- res.df
    res_df$gene_symbol <- rownames(res_df)
    theme_set(theme_classic(base_size = 20) +
                  theme(
                      axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                      axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                      plot.title = element_text(hjust = 0.5)))
    res_df$diffexpressed <- "NA"
    res_df$diffexpressed[res_df$log2FoldChange > 2 & res_df$padj < 0.05] <- "UP"
    res_df$diffexpressed[res_df$log2FoldChange < -2 & res_df$padj < 0.05] <- "DOWN"
    #head(res_df[order(res_df$padj) & res_df$diffexpressed == 'DOWN', ])
    #for down miRNAs
    res_df$delabel <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$log2FoldChange), "gene_symbol"],5), res_df$gene_symbol, NA)
    #for UP miRNAs
    res_df$delabel2 <- ifelse(res_df$gene_symbol %in% head(res_df[order(res_df$log2FoldChange,decreasing=T), "gene_symbol"], 5), res_df$gene_symbol, NA)
    
    ggplot(res_df, aes(x = log2FoldChange, y = -log(padj), col = res_df$diffexpressed),label=delabel) +
        geom_point(alpha = 0.4) +
        scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
        labs(
            title = "Volcano Plot",
            x = "log Fold change",
            y = "-log10 adjusted Pvalue"
        ) +
        geom_vline(xintercept = c(2, -2), col="black",linetype = "dashed") +
        geom_hline(yintercept = -log(0.05), col="black",linetype = "dashed")+
        geom_point(size = 2)+
        scale_color_manual(values = c("#00AFBB", "grey", "limegreen"), # to set the colours of our variable
                           labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
        coord_cartesian(ylim = c(0, 400), xlim = c(-14, 14)) + # since some genes can have minuslog10padj of inf, we set these limits
        labs(color = 'Regulation', #legend_title, 
             x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adj.p-value")) + 
        scale_x_continuous(breaks = seq(-14, 14, 2))+ # to customise the breaks in the x axis
        ggtitle('TNBC vs healthy Controls Volcano Plot')+  # Plot title 
        geom_text_repel(label=res_df$delabel,max.overlaps = Inf)+ # To show all labels 
        geom_text_repel(label=res_df$delabel2,max.overlaps = Inf)
```
###    (a) Heatmap for DEGs
```{r Heatmap}
   meta.sorted=meta[order(meta$Run), ]
    sam_condition <- meta.sorted$Condition
    column_annot <- HeatmapAnnotation(
        Condition = sam_condition,
        col = list(Condition = c("TNBC" = "purple", "Control" = "violet")))
    
    # Generate the heatmap
    Heatmap(
        matrix = dems.exp,
        top_annotation = column_annot,
        row_title = "Top DEGs",
        column_title = "Samples",
        cluster_rows = TRUE,
        cluster_columns = F,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 2),
        column_names_rot = 45,
        column_names_centered = TRUE)  
```
## 10- Enrichment Analysis
### A- Enrichment Analysis
```{r Enrichment Analysis}

DAVID_BP <- read_delim("DAVID_BP.txt", delim = "\t",escape_double = FALSE, trim_ws = TRUE)
    DAVID_CC <- read_delim("DAVID_CC.txt", delim = "\t",escape_double = FALSE, trim_ws = TRUE)
    DAVID_MF <- read_delim("DAVID_MF.txt", delim = "\t",escape_double = FALSE, trim_ws = TRUE)
    DAVID <- rbind(DAVID_BP,DAVID_CC,DAVID_MF)
    ID <- t(data.frame(strsplit(DAVID$Term,"~")))
    DAVID <- data.frame(ID,DAVID) 
    DAVID <- DAVID %>% dplyr::select(ID=X1,Term=X2,Category,adj_pval=FDR,everything(),-Term)
    Res.dems <- res.dems 
    Res.dems$ID <- rownames(res.dems)
    Res.dems <- Res.dems %>% dplyr::select(ID,logFC=log2FoldChange,P.Value=pvalue,adj.P.Val=padj,everything())
    
    circ <- circle_dat(DAVID, Res.dems)
```
### B- Enrichment Analysis Plots
```{r Enrichment Analysis Plots}
## Generate a simple barplot
    GOBar(subset(circ, category == 'GOTERM_BP_DIRECT'))
    GOBar(circ, display = 'multiple', title = 'Z-score Coloured Barplot', zsc.col = c('red', 'white', 'navyblue'))
# Generate the bubble plot with a label threshold of 3
    GOBubble(circ, labels = 3,ID=F,table.col=T)
# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
    GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), bg.col = F,display = 'multiple', labels = 3,ID=F)   
    #GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)
# Reduce redundant terms with a gene overlap >= 0.75...
    reduced_circ <- reduce_overlap(circ, overlap = 0.75)
    GOBubble(reduced_circ, labels = 2.8,ID=F)

#Generate a circular visualization of the results of gene- annotation enrichment analysis
    GOCircle(circ,zsc.col = c('violet', 'white', 'cyan'),lfc.col=c('red','green'))    
    
```
