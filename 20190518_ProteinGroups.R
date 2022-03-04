setwd("K:/Collaborations/Matthew_Kraushar/20190517_Ebp1_pSILAC_AHA/output")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(gplots)
library(data.table)







                    #### proteinGroups table load and preparation ####

PG <- fread("../txt_RQ_allFiles/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PG_noRQ <- fread("../txt_noRQ_allFiles/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)




PG$Gene.names <- sapply(strsplit(PG$Gene.names, ";"), "[", 1)
PG$Majority.protein.IDs <- sapply(strsplit(PG$Majority.protein.IDs, ";"), "[", 1)

PG_noRQ$Gene.names <- sapply(strsplit(PG_noRQ$Gene.names, ";"), "[", 1)
PG_noRQ$Majority.protein.IDs <- sapply(strsplit(PG_noRQ$Majority.protein.IDs, ";"), "[", 1)



#### Make sure both tables have same length
PG <- subset(PG, Majority.protein.IDs %in% PG_noRQ$Majority.protein.IDs)
PG_noRQ <- subset(PG_noRQ, Majority.protein.IDs %in% PG$Majority.protein.IDs)





intensities = c("iBAQ.L.AHA_For",
                "iBAQ.M.AHA_For",
                "iBAQ.H.AHA_For",
                "iBAQ.L.AHA_Rev",
                "iBAQ.M.AHA_Rev",
                "iBAQ.H.AHA_Rev",
                "iBAQ.L.For",
                "iBAQ.M.For",
                "iBAQ.H.For",
                "iBAQ.L.Rev",
                "iBAQ.M.Rev",
                "iBAQ.H.Rev")


ratios=c("Ratio.H.M.normalized.AHA_For",
         "Ratio.H.M.normalized.AHA_Rev",
         "Ratio.H.M.normalized.For",
         "Ratio.H.M.normalized.Rev")






#### Define Requantified ratios and intensities ####

### Ratios
requantified.ratios <- c()
unscrupulous.ratios <- c()
for (i in ratios) {
  PG[paste0("Requantified.", i)] <- ifelse(is.na(PG_noRQ[i]) &
                                             !is.na(PG[i]), T, F)
  
  requantified.ratios <- c(requantified.ratios, paste0("Requantified.", i))
  unscrupulous.ratios <- c(unscrupulous.ratios, paste0("Unscrupulous.", i))
  rm(i)
}



### intensities
requantified.intensities <- c()
for (i in intensities) {
  PG[paste0("Requantified.", i)] <- ifelse(PG_noRQ[i] == 0 &
                                             PG[i] != 0, T, F)
  
  requantified.intensities <- c(requantified.intensities, paste0("Requantified.", i))
  rm(i)
}









#### Unscrupulous Requantification ####



# Select Unscrupulous Requantified ratios (ratios between two channels when both are requantified)
#For
PG$Unscrupulous.Ratio.H.M.normalized.AHA_For <- ifelse((PG$Requantified.iBAQ.H.AHA_For & PG$Requantified.iBAQ.M.AHA_For) & !is.na(PG$Ratio.H.M.AHA_For), T, F)

PG$Unscrupulous.Ratio.H.M.normalized.For <- ifelse((PG$Requantified.iBAQ.H.For & PG$Requantified.iBAQ.M.For) & !is.na(PG$Ratio.H.M.For), T, F)

PG$Unscrupulous.Ratio.H.M.normalized.AHA_Rev <- ifelse((PG$Requantified.iBAQ.H.AHA_Rev & PG$Requantified.iBAQ.M.AHA_Rev) & !is.na(PG$Ratio.H.M.AHA_Rev), T, F)

PG$Unscrupulous.Ratio.H.M.normalized.Rev <- ifelse((PG$Requantified.iBAQ.H.Rev & PG$Requantified.iBAQ.M.Rev) & !is.na(PG$Ratio.H.M.Rev), T, F)












# Delet Unscrupulous ratios
for (i in 1:length(ratios)) {
  PG[,ratios[i]] <- ifelse(PG[,unscrupulous.ratios[i]], NA, PG[,ratios[i]])
  
  rm(i)
}





#### Filters ####

#PG <- PG_noRQ
#filter out contaminants, Rev and only identified by site
PG <- subset(PG, Reverse != "+")
PG <- subset(PG, Potential.contaminant != "+")
PG <- subset(PG, Only.identified.by.site != "+")




#transform Intensity and ratios to Log2 (or 10 if you prefer)
PG[c(intensities, ratios)] = log2(PG[c(intensities, ratios)])
# change Inf and NaN values for NA
is.na(PG[c(intensities, ratios)]) <- sapply(PG[c(intensities, ratios)], is.infinite)
is.na(PG[c(intensities, ratios)]) <- sapply(PG[c(intensities, ratios)], is.nan)



# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(intensities)) {
  cat(intensities[i])
  cat("\t")
  cat("\t")
  cat(nrow(PG[!is.na(PG[intensities[i]]),]))
  cat("\t")
  cat(mean(PG[,intensities[i]], na.rm = T))
#  cat("\t")
#  cat(sum(PG[,requantified.intensities[i]], na.rm = T))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have ratios values not NA) in each L-H group pair
for (i in 1:length(ratios)) {
  #group
  cat(ratios[i])
  cat("\t")
  cat("\t")
  #quantified ratios
  cat(nrow(PG[!is.na(PG[ratios[i]]),]))
  cat("\t")
  #ratio mean
  cat(mean(PG[,ratios[i]], na.rm = T))
  cat("\t")
  #Requantified ratios
  cat(sum(!is.na(PG[ratios[i]]) & PG[, requantified.ratios[i]], na.rm = T))
  cat("\t")
  #Percentage of requantified ratios
  cat(100*(sum(!is.na(PG[ratios[i]]) & PG[, requantified.ratios[i]], na.rm = T)) / nrow(PG[!is.na(PG[ratios[i]]),]))
  cat("\n")
  rm(i)
}


# how many proteins were identified (have intensities values not NA) in each L-H group pair
for (i in 1:length(unscrupulous.ratios)) {
  cat(unscrupulous.ratios[i])
  cat("\t")
  cat("\t")
  cat(sum(PG[,unscrupulous.ratios[i]]))
  cat("\t")
  cat(100*(sum(PG[,unscrupulous.ratios[i]])/sum(!PG[,unscrupulous.ratios[i]])))
  cat("\n")
  rm(i)
}





rm(PG_noRQ, requantified.intensities, requantified.ratios, unscrupulous.ratios)




#### Define PG_summ as working table ####

PG_summ <- subset(PG, select = c("Majority.protein.IDs",
                                 "Gene.names", "Peptide.IDs",
                                 intensities, ratios))

#PG_summ$Gene.names <- sapply(strsplit(PG_summ$Gene.names, ";"), "[", 1)
PG_summ$Majority.protein.IDs <- sapply(strsplit(PG_summ$Majority.protein.IDs, ";"), "[", 1)



#### Pa2g4 was quantified with Requantify and is wrong, based on Western Blot
PG_summ$Ratio.H.M.normalized.For[PG_summ$Gene.names == "Pa2g4"] <- NA
PG_summ$iBAQ.H.For[PG_summ$Gene.names == "Pa2g4"] <- NA



#### Invert Reverse ratios ####
PG_summ$Ratio.H.M.normalized.AHA_Rev <- -PG_summ$Ratio.H.M.normalized.AHA_Rev
PG_summ$Ratio.H.M.normalized.Rev <- -PG_summ$Ratio.H.M.normalized.Rev







#### Write PG_summ ####
write.table(PG_summ, file = "PG_summ.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)








## boxplot ####
MyBoxplot <- function(df, variables,
                      AxisName_y = ("Log2"), LabelNames = variables,
                      limits_y_min = round(min(df[variables], na.rm = T))-1,
                      limits_y_max = round(max(df[variables], na.rm = T))+1,
                      limits_breaks = round((limits_y_max+limits_y_min)/15),
                      ptn = NULL) {
  
  all <- melt(df,
              id.vars = c("Gene.names"), measure.vars=variables)
  
  
  
  bplot <- ggplot(all, aes(variable,value)) +
    
    geom_violin(fill = rgb(0,0,0,.25), na.rm = T) +
    
    geom_boxplot(fill = rgb(0,0,0,0), width = 0.75, na.rm = T, notch = T) +
    
    ylab(AxisName_y) +
    coord_flip(ylim = c(limits_y_min, limits_y_max)) +
    scale_y_continuous(breaks = seq(limits_y_min, limits_y_max, limits_breaks)) +
    scale_x_discrete(labels = LabelNames) +
    
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20),
          axis.title.x = element_text(face="bold",
                                      size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=0.4,
                                      size=20))
  
  if (!is.null(ptn)) {
    
    all2 <- subset(all, Gene.names %in% ptn)
    
    bplot <- bplot +
      
      geom_point(data = all2,
                 aes(x = variable, col = Gene.names),
                 size = 4) +
      
      geom_line(data = all2,
                aes(x = variable, col = Gene.names, group = Gene.names),
                size = 2, linetype = 1)
    
    
  }
  
  
  
  
  return(bplot)
  
}


#select genes you want to see patterns
ptn = "Pa2g4"


plot <- MyBoxplot(PG_summ, intensities, AxisName_y = "Intensities (Log2 iBAQ)",
                  LabelNames = intensities,
                  ptn = ptn)


png(paste("Boxplot_intensities.png", sep = ""), width = 750, height = 500, pointsize = 25)
grid.arrange(plot,
             nrow = 1)
dev.off()
rm(plot, ptn)










#### Scatter Plots ####


MyScatterPlot <- function(df, X, Y, KnockdownName, threshold = 1, axis = 4.5) {
  plot <- ggplot(data = df, aes(x = df[,X], y = df[,Y])) +
    
    
    #ggtitle(label = KnockdownName) +
    ylab(paste0("si", KnockdownName, "/", "siControl (Log2_FC)")) +
    xlab(paste0("si", KnockdownName, "/", "siControl (Log2_FC)")) +
    coord_cartesian(xlim = c(-axis, axis),ylim = c(-axis, axis)) +
    scale_x_continuous(breaks=seq(-axis, axis, 1)) +
    scale_y_continuous(breaks=seq(-axis, axis, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.x = element_text(face="bold",
                                      size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.x  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20),
          axis.title.y = element_text(face="bold",
                                      size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20)) +
    
    annotate(geom = "segment",
             x = -threshold, xend = -2*axis,
             y = -threshold, yend = -threshold,
             colour = rgb(.6,0,.85,1),
             linetype = "dashed",
             size = 1.5,
             lineend = "round") +
    
    annotate(geom = "segment",
             x = -threshold, xend = -threshold,
             y = -threshold, yend = -2*axis,
             colour = rgb(.6,0,.85,1),
             linetype = "dashed",
             size = 1.5,
             lineend = "round") +
    
    
    annotate(geom = "segment",
             x = threshold, xend = 2*axis,
             y = threshold, yend = threshold,
             colour = rgb(.9,.7,0,1),
             linetype = "dashed",
             size = 1.5,
             lineend = "round") +
    
    annotate(geom = "segment",
             x = threshold, xend = threshold,
             y = threshold, yend = 2*axis,
             colour = rgb(.9,.7,0,1),
             linetype = "dashed",
             size = 1.5,
             lineend = "round") +
    
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    
    
    geom_point(na.rm = T,
               shape = 21,
               stroke = 0,
               size = 4,
               fill = densCols(x = df[,X],
                               y = df[,Y],
                               colramp = colorRampPalette(c(rgb(.5,.5,.5),rgb(0,0,0)))),
               color = rgb(0,0,0,0)) +
    
    annotate(geom = "text", label = nrow(subset(df, !is.na(df[,X]) & !is.na(df[,Y]))),
             x = -axis,
             y = axis,
             colour = rgb(.25,.25,.25),
             size = 5) +
    
    geom_point(data = subset(df, Gene.names == KnockdownName),
               aes(x = subset(df, Gene.names == KnockdownName)[,X],
                   y = subset(df, Gene.names == KnockdownName)[,Y]),
               na.rm = T,
               shape = 21,
               stroke = 0,
               size = 4,
               fill = rgb(1,0,0,1),
               color = rgb(0,0,0,0)) +
    
    geom_text_repel(data = subset(df, Gene.names == KnockdownName),
                    aes(x = subset(df, Gene.names == KnockdownName)[,X],
                        y = subset(df, Gene.names == KnockdownName)[,Y],
                        label = Gene.names),
                    size = 5,
                    fontface = "bold",
                    colour = rgb(1,0,0,1),
                    na.rm = T) +
    
    
    ### Increase translation
    geom_point(data = subset(df, df[,X] > threshold & df[,Y] > threshold),
               aes(x = subset(df, df[,X] > threshold & df[,Y] > threshold)[,X],
                   y = subset(df, df[,X] > threshold & df[,Y] > threshold)[,Y]),
               na.rm = T,
               shape = 21,
               stroke = 2,
               size = 3.5,
               fill = rgb(0,0,0,0),
               color = rgb(.9,.7,0,1)) +
    
    annotate(geom = "text", label = nrow(subset(df, df[,X] > log2(1.25) & df[,Y] > log2(1.25) &
                                           (!is.na(df[,X]) & !is.na(df[,Y])))),
             x = -axis+1,
             y = axis-.5,
             colour = rgb(.9,.7,0),
             size = 5) +
    
    annotate(geom = "text", label = nrow(subset(df, df[,X] > log2(1.5) & df[,Y] > log2(1.5) &
                                                  (!is.na(df[,X]) & !is.na(df[,Y])))),
             x = -axis+.5,
             y = axis-.5,
             colour = rgb(.9,.7,0),
             size = 5) +
    
    annotate(geom = "text", label = nrow(subset(df, df[,X] > log2(2) & df[,Y] > log2(2) &
                                                  (!is.na(df[,X]) & !is.na(df[,Y])))),
             x = -axis,
             y = axis-.5,
             colour = rgb(.9,.7,0),
             size = 5) +
    
    geom_text_repel(data = subset(df, df[,X] > threshold & df[,Y] > threshold),
                    aes(x = subset(df, df[,X] > threshold & df[,Y] > threshold)[,X],
                        y = subset(df, df[,X] > threshold & df[,Y] > threshold)[,Y],
                        label = Gene.names),
                    size = 5,
                    fontface = "bold",
                    colour = rgb(0,0,0,1),
                    na.rm = T) +
    
    
    
    
    ### Decrease translation
    geom_point(data = subset(df, df[,X] < -threshold & df[,Y] < -threshold),
               aes(x = subset(df, df[,X] < -threshold & df[,Y] < -threshold)[,X],
                   y = subset(df, df[,X] < -threshold & df[,Y] < -threshold)[,Y]),
               na.rm = T,
               shape = 21,
               stroke = 2,
               size = 3.5,
               fill = rgb(0,0,0,0),
               color = rgb(.6,0,.85,1)) +
    
    annotate(geom = "text", label = nrow(subset(df, df[,X] < -log2(1.25) & df[,Y] < -log2(1.25)  &
                                                  (!is.na(df[,X]) & !is.na(df[,Y])))),
             x = -axis+1,
             y = axis-1,
             colour = rgb(.6,0,.85),
             size = 5) +
    
    annotate(geom = "text", label = nrow(subset(df, df[,X] < -log2(1.5) & df[,Y] < -log2(1.5)  &
                                                  (!is.na(df[,X]) & !is.na(df[,Y])))),
             x = -axis+.5,
             y = axis-1,
             colour = rgb(.6,0,.85),
             size = 5) +
    
    annotate(geom = "text", label = nrow(subset(df, df[,X] < -log2(2) & df[,Y] < -log2(2)  &
                                                  (!is.na(df[,X]) & !is.na(df[,Y])))),
             x = -axis,
             y = axis-1,
             colour = rgb(.6,0,.85),
             size = 5) +
    
    geom_text_repel(data = subset(df, df[,X] < -threshold & df[,Y] < -threshold),
                    aes(x = subset(df, df[,X] < -threshold & df[,Y] < -threshold)[,X],
                        y = subset(df, df[,X] < -threshold & df[,Y] < -threshold)[,Y],
                        label = Gene.names),
                    size = 5,
                    fontface = "bold",
                    colour = rgb(0,0,0,1),
                    na.rm = T)
  
  
  
  
  
  return(plot)
}





pSILAC_AHA <- MyScatterPlot(df = PG_summ, X = ratios[1], Y = ratios[2], KnockdownName = "Pa2g4", threshold = 1, axis = 3) +
  ggtitle("pSILAC AHA")

pSILAC <- MyScatterPlot(df = PG_summ, X = ratios[3], Y = ratios[4], KnockdownName = "Pa2g4", threshold = 1, axis = 3) +
  ggtitle("pSILAC")




png("Ratios_plot.png", width = 1000, height = 500, pointsize = 25)
grid.arrange(pSILAC_AHA, pSILAC, ncol=2)
dev.off()


pdf("Ratios_plot_AHA.pdf")
grid.arrange(pSILAC_AHA)
dev.off()

pdf("Ratios_plot.pdf")
grid.arrange(pSILAC)
dev.off()

rm(pSILAC_AHA, pSILAC)



#Number proteins in AHA experiment
table(!is.na(PG_summ$Ratio.H.M.normalized.AHA_For) & !is.na(PG_summ$Ratio.H.M.normalized.AHA_Rev))

#Number proteins in pSILAC experiment
table(!is.na(PG_summ$Ratio.H.M.normalized.For) & !is.na(PG_summ$Ratio.H.M.normalized.Rev))



#### Venn diagrams ####
# set parameters
threshold <- log2(1.25)

#colect gene names
Total_increase_AHA_pSILAC <- unique(subset(PG_summ, (PG_summ[,"Ratio.H.M.normalized.AHA_For"] > threshold & 
                                                       PG_summ[,"Ratio.H.M.normalized.AHA_Rev"] > threshold) &
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.AHA_For"]) & 
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.AHA_Rev"]))$Gene.names)

Total_increase_pSILAC <- unique(subset(PG_summ, (PG_summ[,"Ratio.H.M.normalized.For"] > threshold & 
                                                       PG_summ[,"Ratio.H.M.normalized.Rev"] > threshold) &
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.For"]) & 
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.Rev"]))$Gene.names)


Total_decrease_AHA_pSILAC <- unique(subset(PG_summ, (PG_summ[,"Ratio.H.M.normalized.AHA_For"] < -threshold & 
                                                       PG_summ[,"Ratio.H.M.normalized.AHA_Rev"] < -threshold) &
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.AHA_For"]) & 
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.AHA_Rev"]))$Gene.names)

Total_decrease_pSILAC <- unique(subset(PG_summ, (PG_summ[,"Ratio.H.M.normalized.For"] < -threshold & 
                                                       PG_summ[,"Ratio.H.M.normalized.Rev"] < -threshold) &
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.For"]) & 
                                             !is.na(PG_summ[,"Ratio.H.M.normalized.Rev"]))$Gene.names)


MyVennDiagram <- function(A, B, C = NULL) {
  require(eulerr)
  
  
  if (is.null(C)) {
    A_B = length(intersect(A, B))
    len_A = length(A) - A_B
    len_B = length(B) - A_B
    
    venn_plot <- euler(c(A = len_A, B = len_B, "A&B" = A_B))
  } else {
    
    A_B_noC = length(intersect(A, B)) - length(intersect(A,intersect(B, C)))
    A_C_noB = length(intersect(A, C)) - length(intersect(A,intersect(B, C)))
    B_C_noA = length(intersect(B, C)) - length(intersect(A,intersect(B, C)))
    A_B_C = length(intersect(A,intersect(B, C)))
    
    len_A = length(A) - A_B_C - A_B_noC - A_C_noB
    len_B = length(B) - A_B_C - A_B_noC - B_C_noA
    len_C = length(C) - A_B_C - A_C_noB - B_C_noA
    
    venn_plot <- euler(c(A = len_A, B = len_B, C = len_C,
                         "A&B" = A_B_noC, "A&C" = A_C_noB, "B&C" = B_C_noA,
                         "A&B&C" = A_B_C))
  }
  
  return(venn_plot)
  
}

pdf(file = paste0("Venn_plot_increase_FC", 2**threshold, ".pdf"), height = 5, width = 5, useDingbats = F)
plot(MyVennDiagram(A = Total_increase_AHA_pSILAC,
                   B = Total_increase_pSILAC),
     legend = TRUE, quantities= TRUE)
dev.off()

pdf(file = paste0("Venn_plot_decrease_FC", 2**threshold, ".pdf"), height = 5, width = 5, useDingbats = F)
plot(MyVennDiagram(A = Total_decrease_AHA_pSILAC,
                   B = Total_decrease_pSILAC),
     legend = TRUE, quantities= TRUE)
dev.off()









#### Writing tables ####

# All quantified proteins in either experiments
foo <- subset(PG_summ, (!is.na(PG_summ$Ratio.H.M.normalized.AHA_For) |
                                       !is.na(PG_summ$Ratio.H.M.normalized.AHA_Rev)) |
                             (!is.na(PG_summ$Ratio.H.M.normalized.For) |
                                !is.na(PG_summ$Ratio.H.M.normalized.Rev)))
write.table(foo$Majority.protein.IDs, file = "allQuant.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = F)




# Proteins quantified in For and Rev
foo <- subset(PG_summ, (!is.na(PG_summ$Ratio.H.M.normalized.AHA_For) &
                                       !is.na(PG_summ$Ratio.H.M.normalized.AHA_Rev)) |
                             (!is.na(PG_summ$Ratio.H.M.normalized.For) &
                                !is.na(PG_summ$Ratio.H.M.normalized.Rev)))
write.table(foo$Majority.protein.IDs, file = "allQuant_duplicate.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = F)



# proteins increasing either in AHA or normal pSILAC
foo_increase <- subset(PG_summ, (PG_summ$Ratio.H.M.normalized.AHA_For > 1 &
                          PG_summ$Ratio.H.M.normalized.AHA_Rev < -1) |
                (PG_summ$Ratio.H.M.normalized.For > 1 &
                   PG_summ$Ratio.H.M.normalized.Rev < -1))
write.table(foo_increase$Majority.protein.IDs, file = "allIncrease.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = F)


# proteins decreasing either in AHA or normal pSILAC
foo_decrease <- subset(PG_summ, (PG_summ$Ratio.H.M.normalized.AHA_For < -1 &
                                   PG_summ$Ratio.H.M.normalized.AHA_Rev > 1) |
                         (PG_summ$Ratio.H.M.normalized.For < -1 &
                            PG_summ$Ratio.H.M.normalized.Rev > 1))
write.table(foo_decrease$Majority.protein.IDs, file = "allDecrease.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = F)


write.table(c(foo_decrease$Majority.protein.IDs, foo_increase$Majority.protein.IDs), file = "allChange.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = F)



rm(foo)


#### GO annotation ####
#Use Uniprot to add annotations to genes
Uniprot <- fread("K:/Datasets/20190122_Uniprot_MM.tab", sep = "\t", stringsAsFactors = F, data.table = F)
#Uniprot_Reviewed <- subset(Uniprot, Status == "reviewed")

PG_summ <- merge(PG_summ, Uniprot[,c("Entry",
                                              "Gene ontology (molecular function)",
                                              "Gene ontology (biological process)",
                                              "Gene ontology (cellular component)")],
                 by.x = "Majority.protein.IDs", by.y = "Entry",
                 all.x = T, all.y = F)

rm(Uniprot)



GO_terms = c("GO:0051823", #regulation of synapse structural plasticity
             "GO:1905704", #positive regulation of inhibitory synapse assembly
             "GO:1905703", #negative regulation of inhibitory synapse assembly
             "GO:1904890", #negative regulation of excitatory synapse assembly
             "GO:1904891") #positive regulation of excitatory synapse assembly

#PG_summ[,"Positive_Regulaters_Synapse_Assembly"] <- ifelse(grepl("GO:1905704",
#                                                                 PG_summ$`Gene ontology (biological process)`) |
#                                                             grepl("GO:1904891",
#                                                                   PG_summ$`Gene ontology (biological process)`),
#                                                           T, F)


#PG_summ[,"Negative_Regulaters_Synapse_Assembly"] <- ifelse(grepl("GO:1905703",
#                                                                 PG_summ$`Gene ontology (biological process)`) |
#                                                             grepl("GO:1904890",
#                                                                   PG_summ$`Gene ontology (biological process)`),
#                                                           T, F)






#### Cumulative plots ####


MyCumulative <- function(df, X,
                         xlim.min = -4, xlim.max = 4, xbreak = 1,
                         ylim.min = -1, ylim.max = 1, ybreak = .1,
                         ylab = NULL, xlab = NULL, titlelab = NULL,
                         sizeMain = 4, colorMain = rgb(0,0,0,1)) {
  
  ggplot(working_table_PEP,
         aes(x = delta_funSILAC_Cl_Mean)) +
    coord_cartesian(xlim = c(xlim.min, xlim.max)) +
    scale_x_continuous(breaks=seq(xlim.min, xlim.max, xbreak)) +
    scale_y_continuous(breaks=seq(ylim.min, ylim.max, ybreak)) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(label = titlelab)
  theme_bw() +
    theme(plot.title = element_text(face="bold",
                                    size = 35,
                                    hjust = 0.5,
                                    vjust = 0.4),
          axis.title.x = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.text.x  = element_text(face = "bold", color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          axis.title.y = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold", color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          panel.grid=element_blank()) +
    
    geom_step(stat="ecdf",                         # Main
              size = sizeMain, 
              color = colorMain,
              na.rm = T)
  
  
}
  
