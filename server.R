library(shiny)  #Loading required packages for the Shiny App
library(broom)
library(survMisc)
library(survminer)
library(shinythemes)
library(DT)
library(dplyr)
library(ggrepel)

shinyServer(function(input, output, session) {  #Begin shiny server portion of app
  
  CancerDatasets <- read.delim("CancerDatasets.txt", check.names = FALSE, row.names = 1)  #Reads a predetermined file containing dataset descriptions
  output$mytable1 <- DT::renderDataTable({
    DT::datatable(CancerDatasets)
  })
  
  datasetInput <- reactive({  #Variable based on user slected cancer, changing this "resets" the app from the beginning
    files[[as.numeric(input$dataset)]]
  }) 
  
  Input2 <- reactive({
    files[[as.numeric(input$dataset2)]]
  }) 
  
  observe({  #Observe allows the app to realize when a user defined variable has been changed, such as changing the protein of interest
    updateSelectInput(session, "Protein", choices = names(datasetInput()[,-c(1:5)]))
  })
  
  observe({
    updateSelectInput(session, "Protein2", choices = list_of_proteins)
  })
  
  plotInput <- reactive({
    data2 <- datasetInput()
    for (j in 1:ncol(data2)) {
      if (class(data2[, j]) %in% c("character")) 
        data2[,j] <- as.factor(data2[,j])
      else data2[,j] = data2[, j]
    }
    res <- cutoff_type(input, data2)  #Calls function "cutoff_type" which is located in the global environment
})
  
  output$kmplot <- renderPlot({
    print(plotInput())
  })  
  
  plotInput2 <- reactive({
    data3 <- Input2()
    names <- t(colnames(data3)[c(-(1:5))])  #Removes first 5 columns of dataset (survival information), proteins begin in column 6
    data3 <- data.frame(data3)
    
      if (input$survival2 == "Overall") {  #Loop for when user selects "Overall" as survival-type of interest
        time <- data3[,"OS_TIME"]          #"OS_TIME" = overall survival time, "OS_IND" = censor status
        censor <- data3[,"OS_IND"]
          if (as.numeric(as.integer(input$cutoff2)) > 1) {  #Creating cut-point based off user selection
           perc <- as.numeric(as.integer(input$cutoff2))
           Groups <- sapply(data3[,-c(1:5)], function(x) ifelse(x < quantile(x, perc/100, na.rm= TRUE), "Low", "High"))
           dat2 <- cbind.data.frame(time, censor, Groups)
       } else {
          cut <- surv_cutpoint(data3, "OS_TIME", "OS_IND", variables=colnames(data3)[c(-(1:5))], minprop = 0.15)
          Groups <- surv_categorize(cut)
          dat2 <- cbind.data.frame(time, censor, Groups[,-c(1:2)])
       }
      }
     if (input$survival2 != "Overall") {  #Loop for if Recurrence Free Survival is selected
        time <- data3[,"RFS_TIME"]
        censor <- data3[,"RFS_IND"]
        if (as.numeric(as.integer(input$cutoff2)) > 1) {
          perc <- as.numeric(as.integer(input$cutoff2))
          Groups <- sapply(data3[,-c(1:5)], function(x) ifelse(x < quantile(x, perc/100, na.rm= TRUE), "Low", "High"))
          dat2 <- cbind.data.frame(time, censor, Groups)
        } else {
          cut <- surv_cutpoint(data3, "OS_TIME", "OS_IND", variables=colnames(data3)[c(-(1:5))], minprop = 0.15)
          Groups <- surv_categorize(cut)
          dat2 <- cbind.data.frame(time, censor, Groups[,-c(1:2)])
        }
     }
    
 Probes <- names(dat2[,-c(1:2)])
   sumtable <- sapply(Probes, function(i) {
        # STRING INTERPOLATION WITH sprintf, THEN CONVERTED TO FORMULA OBJECT
        iformula <- as.formula(sprintf("Surv(time, censor) ~%s", i))   
        
       # RUN MODEL REFERENCING DYNAMIC FORMULA
        z <- coxph(iformula, data=dat2, na.action=na.omit)
        
       # RETURN COEFF MATRIX RESULTS
        tidy(summary(z)[7][[1]]) 
      })

  GGresults <- data.frame(t(sumtable))
  colnames(GGresults) <- c("Comp", "coef", "HR", "SE", "ZSCORE", "PVAL")
  GGresults2 <- GGresults %>%
    mutate(log10P = -log10(as.numeric(PVAL))) %>%
    mutate(HRH = 1/as.numeric(HR))
  rownames(GGresults2) <- names
  GGresults2$names <- row.names(GGresults2)
  GGresults2 <- subset(GGresults2, log10P > 0.1 & HR <= 20)
  
  data4 <- GGresults2 %>%  #Filter which data points to plot, 1.3 == -log10(0.05)
    filter(log10P >= 1.3 & as.numeric(HRH) >= input$adjusthigh | log10P >= 1.3 & as.numeric(HRH) <= input$adjustlow)
  data4 <- data4 %>%
    mutate(CAT = ifelse(HRH > 1, "Poor", "Good"))  #Mutate adds a column with variable Poor or Good which is dependent on the Hazard Ratio, high-vs-low
  tdat2<- t(dat2) #transpose data that contains the survival information
  Filter <- tdat2[rownames(tdat2) %in% data4$names,] #filter the ggresuts by the names in the dat2 file
  ggProp <- reshape2::melt(Filter) #make a data frame for the porportional bar chart
  ggProp <- ggProp %>%
    filter(!is.na(value))
  ggProp <- merge(ggProp, data4, by.x="Var1", by.y = "names", all.x = TRUE)

  
  if (input$prop == "1") {  #Creation of volcano plot with dynamic labeling based on user selected ranges
    
    res2 <- ggplot(subset(GGresults2, log10P > 0.1), aes(y=log10P, x=as.numeric(HRH))) + 
      geom_jitter(width = 0.05, shape=21, alpha=0.8, aes(size=log10P, fill=as.numeric(HRH))) + 
      geom_hline(yintercept = 1.3, lty=2) +
      geom_vline(xintercept = 1) +
      geom_label_repel(data=subset(GGresults2, log10P >= 1.3 & as.numeric(HRH) >= input$adjusthigh | log10P >= 1.3 & as.numeric(HRH) <= input$adjustlow),
                       aes(label=names),size=3, box.padding = unit(0.3, "lines"), segment.alpha=0.0)  + 
      scale_fill_distiller(palette="RdBu",trans = "log10") +
      guides(size=FALSE, fill=FALSE) + 
      xlab("Cox Regression Hazard Ratio") + 
      ylab("-log10(P-value)") + 
      theme_bw() + 
      theme(text=element_text(family="sans"),
            axis.title=element_text(size=18))
    
  }

  if (input$prop == "2") {  #Creates bar graph that shows the high and low sample proportions from the volcano plot
    
    res2 <- ggplot(subset(ggProp, !is.na(value)), aes(x=reorder(Var1, value=="high", sum), fill=value)) +
      geom_bar(position='fill', alpha=0.85, colour="black") +
      labs(y = "Proportion of Samples") + 
      facet_wrap(~CAT, scales="free") + 
      coord_flip() + 
      scale_fill_manual(values = c("#B2182B", "#2166AC")) + 
      theme_bw() + 
      guides(fill=FALSE) + 
      theme(axis.title.y=element_blank(), 
            strip.text = element_text(face="bold", size=12, family="sans"),
            axis.title=element_text(size=18))
  }
  
return(res2)  
    
  })
  
  output$volplot <- renderPlot({
    print(plotInput2())
  })
  
  plotInput3 <- reactive({
    protein_interest <- as.character(input$Protein2)
    universalProtein <- paste("X", protein_interest, sep = "")  #Some probes have an X to start, this gets around that by searching separately
    if(length(protein_interest) == 0) {
      return(NULL)
    }
    survival_type = input$survival3
    if (survival_type == "Overall"){
      col_time = 2
      col_status = 3   #These define which columns to take survival data from based on what the user wants selects 
    } else {
      col_time = 4
      col_status = 5
    }
    
    cancers = grep(protein_interest, proteins_assayed) #This takes user selected protein and finds which cancer probes used said protein
    num_cancers = length(cancers)
    protein_data = data.frame(matrix(NA, 300, num_cancers))
    protein_hazards = data.frame(matrix(NA, 2, num_cancers))    #row 1 = hazard ratio, row 2 = -log10(pvalue), used later to generate plot
    colnames(protein_data) <- gsub(".csv", "", csv[cancers]) 
    colnames(protein_hazards) <- gsub(".csv", "", csv[cancers]) 
    
    #This loop creates data.frame output of the protein of interest with all the columns being the different cancers that contained user selected protein
    for (i in 1:num_cancers){
      
      cancer_i = colnames(protein_data)[i]  #Determine which cancer file to look into
      cancer_i_file = grep(cancer_i, csv)
      protein_search = paste0("^", protein_interest, "$")  #Require search to match protein exactly, finds the column for selected protein in the cancer file
      protein_searchX = paste0("^", universalProtein, "$") #Variables that start with numbers
      tmpfile = files[[cancer_i_file]]   #Create a temporary duplicate of the cancer file
      tmpfile <- data.frame(tmpfile)
      protein_time = tmpfile[,col_time]      #Pull out the survival time and status before eliminating the first 5 columns (column 6 is where proteins begin)
      protein_status = tmpfile[,col_status]
      protein_i_location = grep(protein_search, colnames(tmpfile))   #Find the protein of interest in the cancer file
      protein_i_locationX = grep(protein_searchX, colnames(tmpfile))
      
      if (length(na.omit(protein_time)) == 0) {
        next
      }
      else {
      
        if (length(protein_i_locationX) < 1) {
          cut <- surv_cutpoint(tmpfile, colnames(tmpfile[col_time]), colnames(tmpfile[col_status]), variables = protein_interest, minprop = 0.15) #Finds optimal cut-point between groups
        }
        
        else {
          cut <- surv_cutpoint(tmpfile, colnames(tmpfile[col_time]), colnames(tmpfile[col_status]), variables=universalProtein, minprop = 0.15)
          
        }
      } 
      
      Group <- surv_categorize(cut)  #Groups data based on cut point found above
      results <- cbind.data.frame(protein_time, protein_status, Group[,3])
      colnames(results)[3] <- "Groups"
      Groups <- as.factor(results$Groups)
      res.cox <- coxph(Surv(as.numeric(protein_time), as.numeric(factor(protein_status))) ~ Groups, data = results, na.action = na.omit)  #Function for HR and p-value
      res.summary <- summary(res.cox)
      protein_hazards[1,i] = 1/(res.summary$coefficients[2])      #HR
      protein_hazards[2,i] = -log10(res.summary$coefficients[5])  #-logt10(pvalue)
      
      rm(tmpfile, results, protein_time, protein_status, res.cox, res.summary)  #Remove variables to prevent forced over-writing
    }
    
    hazards = t(protein_hazards)  #Transpose the matrix for easier plotting
    colnames(hazards) <- c("HR", "log10P")
    hazards <- as.data.frame(hazards)
    hazards$names <- rownames(hazards)
    hazards <- subset(hazards, log10P > 0.1 & HR <= 20)
    Cancer_Type = rownames(hazards)
    
    res3 <- ggplot(hazards, aes(x = HR, y = log10P)) +     #Creation of ggplot graph with pre-determined aesthetics
      geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = HR)) + 
      scale_fill_distiller(palette = "RdBu",trans = "log10") +
      geom_hline(yintercept = 1.3, lty = 2) +
      geom_vline(xintercept = 1) +
      guides(size=FALSE, fill=FALSE) + 
      guides(size=FALSE, fill=FALSE) + 
      xlab("Cox Regression Hazard Ratio") + 
      ylab("-log10(P-value)") + 
      geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh1 | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow1),
                       aes(label=names), size=3, box.padding = unit(0.3, "lines"), segment.alpha=0.0)  + 
      theme_bw() + 
      theme(text=element_text(family="sans"),
            axis.title=element_text(size=18))
    return(res3)
  })
  
  output$volplot2 <- renderPlot({
    print(plotInput3())
  })

  idmap <- list("75" = 'Quartile',
                "50" = 'Median',
                "66" = 'Tertile',
                "1" = "Optimal") 
  
output$downloadPlot <- downloadHandler(
  filename = function() { paste(names(list_of_datasets)[[as.numeric(input$dataset)]], paste(input$Protein, collapse = "_"), idmap[[input$cutoff]], input$survival, 'kmPlot.pdf', sep='_') },
  content = function(file) {
    ggsave(file, print(plotInput(), newpage = FALSE), width = 8,
             height = 6, units = "in")
    })

output$downloadPlot2 <- downloadHandler(
  filename = function() { paste(names(list_of_datasets)[[as.numeric(input$dataset2)]], idmap[[input$cutoff2]], input$survival2, 'allPlot.pdf', sep='_') },
  content = function(file) {
    ggsave(file,plotInput2(), width = 6,
           height = 4.5, units = "in")
  }) 
output$downloadPlot3 <- downloadHandler(
  filename = function() {paste(input$Protein2, input$survival3, 'acrossCancers_Plot.pdf', sep='_') },
  content = function(file) {
    ggsave(file,plotInput3(), width = 8,
           height = 6, units = "in")
  }) 
})