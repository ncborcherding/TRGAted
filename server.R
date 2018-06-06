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
  
  observe({
    if (is.na(clinicalvariables$subtype[as.numeric(input$dataset)]) == T){ #I think this will only go further if it sees the clinicalvariable file says it has subtyping
      return(NULL)
    } else {
    updateSelectInput(session, "subtype", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset)]]$subtype))[order(unique(na.omit(files[[as.numeric(input$dataset)]]$subtype)))])) #unique options
    }
  })
  
  observe({
    if (is.na(clinicalvariables$subtype[as.numeric(input$dataset2)]) == T){ #I think this will only go further if it sees the clinicalvariable file says it has subtyping
      return(NULL)
    } else {
      updateSelectInput(session, "subtype2", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset2)]]$subtype))[order(unique(na.omit(files[[as.numeric(input$dataset2)]]$subtype)))])) #unique options
    }
  })
  
  observe({ 
    surv1.types <- which(subtype_info[as.numeric(input$dataset),] == "y")
    surv1.types <- surv1.types[surv1.types == 11 | surv1.types == 13 | surv1.types == 15 | surv1.types == 17]
    updateSelectInput(session, "survival", choices = colnames(subtype_info)[surv1.types])
  })
  
  observe({  #Observe allows the app to realize when a user defined variable has been changed, such as changing the protein of interest
    updateSelectInput(session, "Protein", choices = names(datasetInput()[,-c(1:17)]))
  })

  observe({
    updateSelectInput(session, "stage", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset)]]$ajcc_pathologic_tumor_stage))[order(unique(na.omit(files[[as.numeric(input$dataset)]]$ajcc_pathologic_tumor_stage)))]))
  })
  observe({
    updateSelectInput(session, "histo_type", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset)]]$histological_type))[order(unique(na.omit(files[[as.numeric(input$dataset)]]$histological_type)))]))
  })
  observe({
    updateSelectInput(session, "histo_grade", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset)]]$histological_grade))[order(unique(na.omit(files[[as.numeric(input$dataset)]]$histological_grade)))]))
  })
  observe({
    updateSelectInput(session, "treat_outcome", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset)]]$treatment_outcome_first_course))[order(unique(na.omit(files[[as.numeric(input$dataset)]]$treatment_outcome_first_course)))]))
  })
  
  observe({
    updateSelectInput(session, "age", choices = append("All", unique(na.omit(files[[as.numeric(input$dataset)]]$age))[order(unique(na.omit(files[[as.numeric(input$dataset)]]$age)))]))
  })
  
  observe({
    updateSelectInput(session, "Protein2", choices = list_of_proteins)
  })
  
  observeEvent(input$log2, {
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    x <- as.numeric(as.integer(input$log2))*0.5
  updateNumericInput(session, "adjustlow", value = x, step=0.1)
  })
  
  observeEvent(input$log3, {
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    x <- as.numeric(as.integer(input$log3))*0.5
    updateNumericInput(session, "adjustlow1", value = x, step=0.1)
  })

  plotInput <- reactive({

    if(any(as.numeric(input$dataset) == list_of_subtypes) == T){
      if(input$subtype!="All"){
        data2<-datasetInput()
        subtype.i <- which(data2[,2] == input$subtype)
        data2 <- data2[subtype.i,]
      } else {data2<-datasetInput()}
    } else {data2<-datasetInput()}
    
    if(any(as.numeric(input$dataset) == list_of_stages) == T){
      if(any(input$stage!="All") == T){
        stage_names <- plyr::count(data2$ajcc_pathologic_tumor_stage)
        stage_use <- which((stage_names$freq > 10) == T)
        stage_names <- droplevels(na.omit(stage_names$x[stage_use]))
        stage.i <- which(data2[,6] %in% input$stage)
        validate(need(length(stage.i) >= 10, paste("This selection of stage(s) produces a cohort in size less than 10. Please use All samples, a different combination of samples, or select", paste(stage_names, collapse = " , "))))
        data2 <- data2[stage.i,]
      } else {data2<-data2}
    } else {data2<-data2}
    
    if(any(as.numeric(input$dataset) == list_of_histo) == T){
      if(any(input$histo_type!="All") == T){
        type_names <- plyr::count(data2$histological_type)
        type_use <- which((type_names$freq > 10) == T)
        type_names <- droplevels(na.omit(type_names$x[type_use]))
        type.i <- which(data2[,7] %in% input$histo_type)
        validate(need(length(type.i) >= 10, paste("This selection of type(s) produces a cohort in size less than 10. Please use All samples, a different combination of samples, or select", paste(type_names, collapse = " , "))))
        data2 <- data2[type.i,]
      } else {data2<-data2}
    } else {data2<-data2}
    
    if(any(as.numeric(input$dataset) == list_of_grade) == T){
      if(any(input$histo_grade!="All") == T){
        grade_names <- plyr::count(data2$histological_grade)
        grade_use <- which((grade_names$freq > 10) == T)
        grade_names <- droplevels(na.omit(grade_names$x[grade_use]))
        grade.i <- which(data2[,8] %in% input$histo_grade)
        validate(need(length(grade.i) >= 10, paste("This selection of grade(s) produces a cohort in size less than 10. Please use All samples, a different combination of samples, or select", paste(grade_names, collapse = " , "))))
        data2 <- data2[grade.i,]
      } else {data2<-data2}
    } else {data2<-data2}
    
    if(any(as.numeric(input$dataset) == list_of_treat) == T){
      if(any(input$treat_outcome!="All") == T){
        treat_names <- plyr::count(data2$treatment_outcome_first_course)
        treat_use <- which((treat_names$freq > 10) == T)
        treat_names <- droplevels(na.omit(treat_names$x[treat_use]))
        treat.i <- which(data2[,9] %in% input$treat_outcome)
        validate(need(length(treat.i) >= 10, paste("This selection of treatment outcome(s) produces a cohort in size less than 10. Please use All samples, a different combination of samples, or select", paste(treat_names, collapse = " , "))))
        data2 <- data2[treat.i,]
      } else {data2<-data2}
    } else {data2<-data2}
    
    if(input$gender!="All"){
      gender_names <- plyr::count(data2$gender)
      gender_use <- which((gender_names$freq > 10) == T)
      gender_names <- droplevels(na.omit(gender_names$x[gender_use]))
      gender.i <- which(data2[,4] %in% input$gender)
      validate(need(length(gender.i) >= 10, paste("This selection of this gender only produces a cohort in size less than 10. Please use All samples or a different combination of samples")))
      data2 <- data2[gender.i,]
    }else {data2<-data2}
    
    if(any(input$age!="All") == T){
      age_names <- plyr::count(data2$age)
      age_use <- which((age_names$freq > 10) == T)
      age_names <- droplevels(na.omit(age_names$x[age_use]))
      age.i <- which(data2[,3] %in% input$age)
      validate(need(length(age.i) >= 10, paste("This selection of this age group only produces a cohort in size less than 10. Please use All samples or a different combination of samples")))
      data2 <- data2[age.i,]
    }else {data2<-data2}

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
  
  dataInput2 <- reactive({
    data3 <- Input2()
    names <- t(colnames(data3)[c(-(1:17))])  #Removes first 5 columns of dataset (survival information), proteins begin in column 6
    data3 <- data.frame(data3)
    
      if(input$subtype2!="All"){
        subtype<-data3[,1:2]
        pts.df<-as.character(data3[,1])
        pts.subtype<-as.character(subtype$patient)
        subtype_class<-subtype[match(pts.df, pts.subtype),2]
        data3<-data.frame(data3, subtype_class)
        subset<-subset(data3, subtype_class==input$subtype2)
        data3<-subset
        data3<-data3[,-ncol(data3)]  
        colnames(data3)[1]<-""
      
      
    }  
    if (input$survival2 == "Overall") {  #Loop for when user selects "Overall" as survival-type of interest
      time <- data3[,"OS.time"]          #"OS.time" = overall survival time, "OS.IND" = censor status
      censor <- data3[,"OS"]
      if (as.numeric(as.integer(input$cutoff2)) > 1) {  #Creating cut-point based off user selection
        perc <- as.numeric(as.integer(input$cutoff2))
        Groups <- sapply(data3[,-c(1:17)], function(x) ifelse(x < quantile(x, perc/100, na.rm= TRUE), "Low", "High"))
        dat2 <- cbind.data.frame(time, censor, Groups)
      } else {
        cut <- surv_cutpoint(data3, "OS.time", "OS", variables=colnames(data3)[c(-(1:17))], minprop = 0.15)
        Groups <- surv_categorize(cut)
        dat2 <- cbind.data.frame(time, censor, Groups[,-c(1:2)])
      }
    }
    if (input$survival2 == "Disease-Specific") {  #Loop for if Recurrence Free Survival is selected
      time <- data3[,"DSS.time"]
      censor <- data3[,"DSS"]
      if (as.numeric(as.integer(input$cutoff2)) > 1) {
        perc <- as.numeric(as.integer(input$cutoff2))
        Groups <- sapply(data3[,-c(1:17)], function(x) ifelse(x < quantile(x, perc/100, na.rm= TRUE), "Low", "High"))
        dat2 <- cbind.data.frame(time, censor, Groups)
      } else {
        cut <- surv_cutpoint(data3, "DSS.time", "DSS", variables=colnames(data3)[c(-(1:17))], minprop = 0.15)
        Groups <- surv_categorize(cut)
        dat2 <- cbind.data.frame(time, censor, Groups[,-c(1:2)])
      }
    }
    
    if (input$survival2 == "Disease-Free") {  #Loop for if Disease-free interval is selected
      time <- data3[,"DFI.time"]
      censor <- data3[,"DFI"]
      if (as.numeric(as.integer(input$cutoff2)) > 1) {
        perc <- as.numeric(as.integer(input$cutoff2))
        Groups <- sapply(data3[,-c(1:17)], function(x) ifelse(x < quantile(x, perc/100, na.rm= TRUE), "Low", "High"))
        dat2 <- cbind.data.frame(time, censor, Groups)
      } else {
        cut <- surv_cutpoint(data3, "DFI.time", "DFI", variables=colnames(data3)[c(-(1:17))], minprop = 0.15)
        Groups <- surv_categorize(cut)
        dat2 <- cbind.data.frame(time, censor, Groups[,-c(1:2)])
      }
    }
    
    if (input$survival2 == "Progression-Free") {  #Loop for if progression-free interval is selected
      time <- data3[,"PFI.time"]
      censor <- data3[,"PFI"]
      if (as.numeric(as.integer(input$cutoff2)) > 1) {
        perc <- as.numeric(as.integer(input$cutoff2))
        Groups <- sapply(data3[,-c(1:17)], function(x) ifelse(x < quantile(x, perc/100, na.rm= TRUE), "Low", "High"))
        dat2 <- cbind.data.frame(time, censor, Groups)
      } else {
        cut <- surv_cutpoint(data3, "PFI.time", "PFI", variables=colnames(data3)[c(-(1:17))], minprop = 0.15)
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
      mutate(HR = 1/as.numeric(HR))
    rownames(GGresults2) <- names
    GGresults2$names <- row.names(GGresults2)
    GGresults2 <- subset(GGresults2, log10P > 0.1 & HR <= 20)
    
    if (input$log2 == "-4") {
      GGresults2$HR <- log(GGresults2$HR)
    }
    
    data4 <- GGresults2 %>%  #Filter which data points to plot, 1.3 == -log10(0.05)
      filter(log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow)
    data4 <- data4 %>%
      mutate(CAT = ifelse(HR > 1, "Poor", "Good"))  #Mutate adds a column with variable Poor or Good which is dependent on the Hazard Ratio, high-vs-low
    tdat2<- t(dat2) #transpose data that contains the survival information
    rownames(tdat2)[3:nrow(tdat2)]<-names
    Filter <- tdat2[rownames(tdat2) %in% data4$names,] #filter the ggresuts by the names in the dat2 file
    ggProp <- reshape2::melt(Filter) #make a data frame for the porportional bar chart
    ggProp <- ggProp %>%
      filter(!is.na(value))
    ggProp <- merge(ggProp, data4, by.x="Var1", by.y = "names", all.x = TRUE)
    
    
    if (input$prop == "1") {  #Creation of volcano plot with dynamic labeling based on user selected ranges
      dataset<-GGresults2
    }
    
    if (input$prop == "2") {  #Creates bar graph that shows the high and low sample proportions from the volcano plot
      dataset<-ggProp
    }
    dataset    
  })
  
  plotInput2 <- reactive({
    data3 <- Input2()
    names <- t(colnames(data3)[c(-(1:17))])
    
    dataset<-dataInput2()
      if (input$log2 == "1") {
        if (input$prop == "1") {  #Creation of volcano plot with dynamic labeling based on user selected ranges
        res2 <- ggplot(subset(dataset, log10P > 0.1), aes(y=log10P, x=as.numeric(HR))) + 
          geom_jitter(width = 0.05, shape=21, alpha=0.8, aes(size=log10P, fill=as.numeric(HR))) + 
          geom_hline(yintercept = 1.3, lty=2) +
          geom_vline(xintercept = 1) +
          geom_label_repel(data=subset(dataset, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow),
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
          
          res2 <- ggplot(subset(dataset, !is.na(value)), aes(x=reorder(Var1, value=="high", sum), fill=value)) +
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
      }
      else {
        if (input$prop == "1") {
          res2 <- ggplot(subset(dataset, log10P > 0.1), aes(y=log10P, x=as.numeric(HR))) + 
            geom_jitter(width = 0.05, shape=21, alpha=0.8, aes(size=log10P, fill=as.numeric(HR))) + 
            geom_hline(yintercept = 1.3, lty=2) +
            geom_vline(xintercept = 0) +
            geom_label_repel(data=subset(dataset, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow),
                             aes(label=names),size=3, box.padding = unit(0.3, "lines"), segment.alpha=0.0)  + 
            scale_fill_distiller(palette="RdBu") +
            guides(size=FALSE, fill=FALSE) + 
            xlab("Natural Log of Cox Regression Hazard Ratio") + 
            ylab("-log10(P-value)") + 
            theme_bw() + 
            theme(text=element_text(family="sans"),
                  axis.title=element_text(size=18))
        }
        if (input$prop == "2") {  #Creates bar graph that shows the high and low sample proportions from the volcano plot
          
          res2 <- ggplot(subset(dataset, !is.na(value)), aes(x=reorder(Var1, value=="high", sum), fill=value)) +
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
      }  
  
    return(res2)  
    
  })  
  
  output$volplot <- renderPlot({
    print(plotInput2())
  })
  
  dataInput3 <- reactive({
    protein_interest <- as.character(input$Protein2)
    universalProtein <- paste("X", protein_interest, sep = "")  #Some probes have an X to start, this gets around that by searching separately
    if(length(protein_interest) == 0) {
      return(NULL)
    }
    survival_type = input$survival3
    if (survival_type == "Overall"){
      col_time = 11
      col_status = 10   #These define which columns to take survival data from based on what the user wants selects 
    } 
    if (survival_type == "Disease-Specific"){
      col_time = 13
      col_status = 12
    }
    if (survival_type == "Disease-Free"){
      col_time = 15
      col_status = 14
    }
    if (survival_type == "Progression-Free"){
      col_time = 17
      col_status = 16
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
    
  })
  
  plotInput3 <- reactive({
    hazards <- dataInput3()
    if (input$log3 == "-4") {
      hazards$HR <- log(hazards$HR)
      if (input$bar == "1") {  #Creation of volcano plot with dynamic labeling based on user selected ranges  
        res3 <- ggplot(hazards, aes(x = HR, y = log10P)) +     #Creation of ggplot graph with pre-determined aesthetics
          geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = HR)) + 
          scale_fill_distiller(palette = "RdBu", limits = c(min(hazards$HR), max(hazards$HR))) +
          geom_hline(yintercept = 1.3, lty = 2) +
          geom_vline(xintercept = 0) +
          guides(size=FALSE, fill=FALSE) + 
          guides(size=FALSE, fill=FALSE) + 
          xlab("Natural Log of Cox Regression Hazard Ratio") + 
          ylab("-log10(P-value)") + 
          geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh1 | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow1),
                           aes(label=names), size=5, box.padding = unit(0.4, "lines"), segment.alpha=0.8, segment.size = 0.5)  + 
          theme_bw() + 
          theme(text=element_text(family="sans"),
                axis.title=element_text(size=18))
      }
    
      if (input$bar == "2") {  #Creates bar graph that shows the high and low sample proportions from the volcano plot
        
        res3 <- ggplot(subset(hazards, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh1 | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow1), aes(x=reorder(names, HR), y=HR, fill=HR)) +
          geom_bar(stat = "identity", color="black", alpha=0.8) + 
          scale_fill_distiller(palette = "RdBu", limits = c(min(hazards$HR), max(hazards$HR))) +
          theme_bw() + 
          coord_flip() + 
          geom_hline(yintercept = 0) + 
          guides(fill=FALSE)+
          xlab("Natural Log of Cox Regression Hazard Ratio") + 
          theme(axis.title.y=element_blank(), 
                strip.text = element_text(face="bold", size=12, family="sans"),
                axis.title=element_text(size=18))
      }
     return(res3)
    }
    if (input$log3 == "1") {
      if (input$bar == "1") {  #Creation of volcano plot with dynamic labeling based on user selected ranges  
        res3 <- ggplot(hazards, aes(x = HR, y = log10P)) +     #Creation of ggplot graph with pre-determined aesthetics
          geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = HR)) + 
          scale_fill_distiller(palette = "RdBu", trans="log", limits = c(min(hazards$HR), max(hazards$HR))) +
          geom_hline(yintercept = 1.3, lty = 2) +
          geom_vline(xintercept = 1) +
          guides(size=FALSE, fill=FALSE) + 
          guides(size=FALSE, fill=FALSE) + 
          xlab("Cox Regression Hazard Ratio") + 
          ylab("-log10(P-value)") + 
          geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh1 | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow1),
                           aes(label=names), size=5, box.padding = unit(0.4, "lines"), segment.alpha=0.8, segment.size = 0.5)  + 
          theme_bw() + 
          theme(text=element_text(family="sans"),
                axis.title=element_text(size=18))
      }
      
      if (input$bar == "2") {  #Creates bar graph that shows the high and low sample proportions from the volcano plot
        
        res3 <- ggplot(subset(hazards, log10P >= 1.3 & as.numeric(HR) >= input$adjusthigh1 | log10P >= 1.3 & as.numeric(HR) <= input$adjustlow1), aes(x=reorder(names, HR), y=HR, fill=HR)) +
          geom_bar(stat = "identity", color="black", alpha=0.8) + 
          scale_fill_distiller(palette = "RdBu", trans="log", limits = c(min(hazards$HR), max(hazards$HR))) +
          theme_bw() + 
          coord_flip() + 
          geom_hline(yintercept = 0) + 
          xlab("Cox Regression Hazard Ratio") + 
          guides(fill=FALSE)+
          theme(axis.title.y=element_blank(), 
                strip.text = element_text(face="bold", size=12, family="sans"),
                axis.title=element_text(size=18))
      }
      return(res3)
    }
  })
  
  output$volplot2 <- renderPlot({
    print(plotInput3())
  })

  idmap <- list("75" = 'Quartile',
                "50" = 'Median',
                "66" = 'Tertile',
                "1" = "Optimal") 
  
output$downloadPlot <- downloadHandler(
  filename = function() { paste(names(list_of_datasets)[[as.numeric(input$dataset)]], "_", paste(input$Protein, collapse = "_"), "_", idmap[[input$cutoff]], "_", input$survival, "_", 'kmPlot', as.character(input$plotformat), sep='') },
  content = function(file) {
    ggsave(file, print(plotInput(), newpage = FALSE), width = 8,
             height = 6, units = "in")
    })

output$downloadPlot2 <- downloadHandler(
  filename = function() { paste(names(list_of_datasets)[[as.numeric(input$dataset2)]],"_", idmap[[input$cutoff2]], "_", input$survival2, "_", 'allPlot', as.character(input$plotformat2), sep='') },
  content = function(file) {
    ggsave(file,plotInput2(), width = 6,
           height = 4.5, units = "in")
  }) 
output$downloadPlot3 <- downloadHandler(
  filename = function() {paste(input$Protein2, "_", input$survival3, "_", 'acrossCancers_Plot', as.character(input$plotformat3), sep='') },
  content = function(file) {
    ggsave(file,plotInput3(), width = 8,
           height = 6, units = "in")
  }) 

output$downloadData <- downloadHandler(
  
  filename = function() {
    paste(names(list_of_datasets)[[as.numeric(input$dataset)]],"_", input$Protein, ".csv", sep='')
  },
  content = function(file) {
    #write.csv(apply(datasetInput(),2,as.character)[,c(1,10:17,(match(input$Protein, proteins_assayed[,as.numeric(input$dataset)])))], file, row.names=FALSE)
    #write.csv(apply(datasetInput(),2,as.character)[,c(1:17,(17 + match(input$Protein, proteins_assayed[,as.numeric(input$dataset)])))], file, row.names=FALSE)
    write.csv(apply(datasetInput(),2,as.character)[,c(1:17,(17 + match(input$Protein, proteins_assayed[,as.numeric(input$dataset)])))], file, row.names=FALSE)
  }) 

output$downloadData2 <- downloadHandler(
  
  filename = function() {
    paste(names(list_of_datasets)[[as.numeric(input$dataset2)]],"_", idmap[[input$cutoff2]], "_", input$survival2, ".csv", sep='')
  },
  content = function(file) {
    write.csv(apply(dataInput2(),2,as.character), file, row.names=FALSE)
  }) 


output$downloadData3 <- downloadHandler(
  
  filename = function() {
    paste(input$Protein2, "_", input$survival3, "_", 'acrossCancers.csv', sep='')
  },
  content = function(file) {
    write.csv(apply(dataInput3(),2,as.character), file, row.names=FALSE)
  }) 
})