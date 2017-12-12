library(shiny)
library(broom)
library(survMisc)
library(survminer)
library(shinythemes)
library(DT)
library(dplyr)
library(ggrepel)

#Identifying and loading local data frames, there are currently 31 datasets with 7678 patient samples. 
#Each data frame is set up with the first 5 columns as sampleIDs, OS_TIME, OS_IND, RFS_TIME, and RFS_IND, resepctively. 
csv <- list.files(pattern = ".csv")

files <- vector("list", length(csv))
proteins_assayed = data.frame(matrix(NA, 300, length(csv)))  #length of 300 can be any arbitary length, just needs to be long enough for the cancer with most proteins asssayed
colnames(proteins_assayed) <- csv

for (i in seq_along(files)){
  files[[i]] <- read.csv(csv[i], stringsAsFactors = FALSE, check.names = FALSE)
  proteins_assayed[1:length(colnames(files[[i]])),i] = colnames(files[[i]])
}
#proteins_assayed include columns beyond the first 5 that include the sample and survival information
proteins_assayed <- proteins_assayed[-c(1:5),]
proteins_assayed <- unique(proteins_assayed) #shortens the matrix by eliminating some of the NA's
colnames(proteins_assayed) <- gsub(".csv", "", colnames(proteins_assayed))
list_of_proteins = unique(na.omit(unlist(proteins_assayed)))
list_of_datasets <<- seq_along(files) #this populates the drop down menu for cancer types
names(list_of_datasets) <- gsub(pattern = ".csv", replacement = "", x = csv)

#make_plot generates the different variations of the survival curves using the survminer ggsurplot function
make_plot <- function(fit, datOS, input, Group, col, cox, use) {
  #if only one protein is selected
  if (length(input$Protein) == 1){ 
  res <- ggsurvplot(fit, data = datOS,
                    tables.y.text = FALSE,
                    tables.theme = theme_cleantable(),
                    xlab = "Time in days",
                    pval = T,
                    pval.coord = c(0.05, 0.05),
                    legend.title = input$Protein,
                    legend.labs = levels(Group),
                    risk.table = TRUE, 
                    fontsize = 5,
                    palette = col,
                    ggtheme = theme_bw(base_family = "sans", base_size = 20))
  }
  #if more than one protein is selected
  else if (length(input$Protein) >= 1) {
    res <- ggsurvplot(fit, data = datOS,
                      tables.y.text = FALSE,
                      tables.theme = theme_cleantable(),
                      xlab = "Time in days",
                      pval = T,
                      pval.coord = c(0.05, 0.05),
                      legend.labs = levels(Group),
                      risk.table = TRUE, 
                      fontsize = 5,
                      palette = col,
                      ggtheme = theme_bw(base_family = "sans", base_size = 20))
  }
  if (use == 1){
    res$plot <- res$plot + ggplot2::annotate("text", 
                                             x = Inf, y = Inf, hjust = 1.25, vjust = 1.5, # x and y coordinates of the text
                                             label = paste("HR =", signif(1/summary(cox)$coef[,2],3), sep = " "), size = 5)
  }
  
  return(res)
}

cutoff_type <- function(input,data2) {
  if(length(input$Protein) == 0) {
    return(NULL)
  }
  if (input$survival == "Overall") {
    time <- data2[,"OS_TIME"]
    censor <- data2[,"OS_IND"]
  } else {
    time <- data2[,"RFS_TIME"]
    censor <- data2[,"RFS_IND"]
  }
  #xvar is either the single protein selected or the mean of all the proteins selected
  if (length(input$Protein) == 1){
    xvar <- data2[, input$Protein]
  } else {
    xvar <- rowMeans(data2[, input$Protein])
  }
  
  dat <- cbind.data.frame(time, censor, xvar)
  dat <- na.omit(dat)
  
  if (class(xvar) %in% c("integer", "numeric")) {
    if (as.numeric(as.integer(input$cutoff)) == 1) {
      cut <- surv_cutpoint(dat, "time", "censor", variables="xvar", minprop = 0.15)
      Group <- surv_categorize(cut)
      datOS <- cbind.data.frame(dat[,c(1:2)], Group[,-c(1:2)])
      colnames(datOS)[3] <- "Group"
      Group <- as.factor(datOS$Group)
      col = c("#B2182B", "#2166AC")
      
      fit <- survfit(Surv(as.numeric(time), as.numeric(factor(censor))) ~ Group, data = datOS)
      cox <- coxph(Surv(as.numeric(time), as.numeric(factor(censor))) ~ Group, data = datOS)
      
      res <- make_plot(fit, datOS, input, Group, col, cox, 1)
      
    } else if (as.numeric(as.integer(input$cutoff)) == 75) {
      dat$Group <- ntile(dat[, 'xvar'], 4)
      dat$Group[which(dat$Group=="1")] <- "1st (Low)"
      dat$Group[which(dat$Group=="2")] <- "2nd"
      dat$Group[which(dat$Group=="3")] <- "3rd"
      dat$Group[which(dat$Group=="4")] <- "4th (High)"
      Group <- factor(dat$Group, levels=c("4th (High)", "3rd", "2nd", "1st (Low)"))
      Group <- as.factor(dat$Group)
      col = c("#2166AC", "#92C5DE", "#F4A582", "#B2182B")
      
      
      fit <- survfit(Surv(as.numeric(time), as.numeric(factor(censor))) ~ Group, data = dat)
      
      res <- make_plot(fit, dat, input, Group, col, 1, 0)
      
    } else if (as.numeric(as.integer(input$cutoff)) == 66) {
      
      dat$Group <- ntile(dat[, 'xvar'], 3)
      dat$Group[which(dat$Group=="1")] <- "Low"
      dat$Group[which(dat$Group=="2")] <- "Med"
      dat$Group[which(dat$Group=="3")] <- "High"
      Group <- factor(dat$Group, levels=c("High", "Med", "Low"))
      col=c("#B2182B", "#999999", "#2166AC")
      
      fit <- survfit(Surv(as.numeric(time), as.numeric(factor(censor))) ~ Group, data = dat)
      
      res <- make_plot(fit, dat, input, Group, col, 1, 0)
      
    } else {perc <- as.numeric(as.integer(input$cutoff))
    
    dat$Group <- ifelse(dat[, 'xvar'] < quantile(dat[, 'xvar'], perc/100, na.rm= TRUE), "Low", "High")
    Group <- as.factor(dat$Group)
    col = c("#B2182B", "#2166AC")
    
    fit <- survfit(Surv(as.numeric(time), as.numeric(factor(censor))) ~ Group, data = dat)
    cox <- coxph(Surv(as.numeric(time), as.numeric(factor(censor))) ~ Group, data = dat)
    
    res <- make_plot(fit, dat, input, Group, col, cox, 1)
    }
  }
  return(res)
}