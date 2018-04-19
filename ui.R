library(shiny)
library(broom)
library(survMisc)
library(survminer)
library(shinythemes)
library(DT)
library(dplyr)
library(ggrepel)


# ====================================================================================
tagList(
  navbarPage(
    theme = shinythemes::shinytheme("sandstone"),
    "TRGAted",
    tabPanel("Home",
             fluidRow(
               column(12, offset = 1,
                      mainPanel(
                                br(),
                                br(),
                                h3("Methods:"), "Survival Analysis performed in TRGAted is utilizes the survival (v2.41-3) and survminer (v0.4.1) R packages. Optimal cutpoints are based on the surv_cutpoint function in the survminer package,
                                finding the lowest log-rank p-value with the minimum proportional comparison set to 15% versus 85% of samples. Hazard ratio are derived from the Cox Proportional Hazard regression model and are based on the 
                                High-versus-Low comparison.", strong("Important Note:"), "With the smaller datasets in the TCGA, the methodology above can be biased by the law of small numbers, inflating hazard ratio, be careful with interpretations 
                                and conclusions. Currently, to counteract potential extreme biases based on small numbers, hazard ratios above 20 are automatically filtered from the data.",
                                br(),
                                br(),
                                h3("About the Data:"), "Level 4 data from the reverse-phase protein arrays for each cancer type were downloaded from the", a("TCPA Portal", href="http://tcpaportal.org/tcpa/"), "(Date Downloaded: 11/10/17). 
                                Data was scaled into z-scores and overall (OS) and recurrence-free survival (RFS) information was attached for each sample. Recurrence-free survival is defined by the NCI as 'the length of time after primary treatment for 
                                a cancer ends, that a patient survives without any signs or symptoms of that cancer.' Each dataset varies in size in terms of total samples quantified for protein and by survival information. A summary of the each 
                                dataset is found below, including total number of samples, samples with overall survival information, samples with recurrence-free survival information and total number of proteins quantified:",
                                br(),
                                br(),
                                DT::dataTableOutput("mytable1")
                                
                      )
               )
             )
    ),
      tabPanel("Survival Curve",
             sidebarPanel(
            selectInput(
            inputId = "dataset", 'Select Cancer', choices = list_of_datasets,
            selectize = FALSE
            ),
            conditionalPanel(condition = "input.dataset==3", ##DV added
                             selectInput("bc_subtype", "Select Subtype", choices = c("all_breast_carcinoma", "Basal", "LumA", "LumB", "Her2"), selectize = FALSE)), ## dv added 
            selectInput(
              inputId = "survival", 'Select Survival', choices = c("Overall", "Recurrence-Free"),
              selectize = FALSE
            ),
            selectizeInput(
              inputId = "Protein", 'Select Protein(s)', choices = "placeholder1", multiple = TRUE
            ),
            radioButtons("cutoff", "Select Cut-off", choiceNames = list("Quartile","Tertile", "Median(+HR)","Optimal Cutoff(+HR)"), choiceValues = c(75,66,50,1), selected = 75),
            downloadButton('downloadPlot', 'Download Plot'),
            radioButtons("plotformat", label = "", choices=list(".pdf", ".png", ".svg"), inline=TRUE, selected=".pdf")),
            mainPanel(plotOutput("kmplot", height= 600, width = 800), value =2)
    ),
    
    tabPanel("Across Cancer",
      sidebarPanel(
      selectInput(
        inputId = "dataset2", 'Select Cancer', choices = list_of_datasets,
        selectize = FALSE
      ),
      conditionalPanel(condition = "input.dataset2==3", ##DV added
                       selectInput("bc_subtype2", "Select Subtype", choices = c("all_breast_carcinoma", "Basal", "LumA", "LumB", "Her2"), selectize = FALSE)), ## DV added
      selectInput(
        inputId = "survival2", 'Select Survival', choices = c("Overall", "Recurrence-Free"),
        selectize = FALSE
      ),
      radioButtons("cutoff2", "Select Cut-off", list("Quartile" = 75, "Tertile" = 66, "Median" = 50, "Optimal Cutoff" = 1), selected = 1),
      numericInput("adjusthigh", "Adjust Labeling of Poor Prognosis Markers:", min = NA, max = NA, value = 2.0, step=0.1),
      numericInput("adjustlow", "Adjust Labeling of Good Prognosis Markers:", min = NA, max = NA, value = 0.5, step=0.1),
      radioButtons("log2", "Log Transform", list("No" = 1, "Yes" = -4), selected = 1),
      radioButtons("prop", "Show Proportion", list("No" = 1, "Yes" = 2), selected = 1),

      div(style="display: inline-block;vertical-align:top; width: 150px;",downloadButton('downloadPlot2', 'Download Plot')), 
      div(style="display: inline-block;vertical-align:top; width: 150px;", conditionalPanel(condition = "input.prop==1",
        downloadButton('downloadData2', 'Download Data'))),
      radioButtons("plotformat2", label = "", choices=list(".pdf", ".png", ".svg"), inline=TRUE, selected=".pdf")
      ),
      mainPanel(plotOutput("volplot", height= 600, width = 800), value =2,
              br(),
              br(),
              strong("Note:"), "This may take some time depending on the size of the dataset examined. The adaptive labeling settings are based on Cox Hazard Ratio (HR) and will not label proteins that do not reach a significance 
              threshold of P<0.05. If you are experiencing an error message, adjust the cut-off and labeling strategy in order to have variables to graph."
                )
               ),
    
    tabPanel("Across Protein",
             sidebarPanel(
               selectInput(
                 inputId = "survival3", 'Select Survival', choices = c("Overall", "Recurrence-Free"),
                 selectize = FALSE
               ),
               selectizeInput(
                 inputId = "Protein2", 'Select Protein', choices = "placeholder1" 
               ),
               numericInput("adjusthigh1", "Adjust Labeling of Poor Prognosis Markers:",
                           min = NA, max = NA,
                           value = 2.0, step=0.1),
               numericInput("adjustlow1", "Adjust Labeling of Good Prognosis Markers:",
                           min = NA, max = NA,
                           value = 0.5, step=0.1),
               radioButtons("bar", "Show Barchart", list("No" = 1, "Yes" = 2), selected = 1),
               radioButtons("log3", "Log Transform", list("No" = 1, "Yes" = -4), selected = 1),
               div(style="display: inline-block;vertical-align:top; width: 150px;", downloadButton('downloadPlot3', 'Download Plot')),
               div(style="display: inline-block;vertical-align:top; width: 150px;", downloadButton('downloadData3', 'Download Data')),
               
               radioButtons("plotformat3", label = "", choices=list(".pdf", ".png", ".svg"), inline=TRUE, selected=".pdf")),
             mainPanel(plotOutput("volplot2", height= 600, width = 800), value =2,
                       br(),
                       br()
                       
             )
    ),
    tabPanel("How to Guide",
             fluidRow(
               column(12, offset = 1,
                      mainPanel(h1("Survival Curve"),
                                img(src="tab1_1.png", width=1035, height=521),
                                br(),
                                "The survival curve is generated by selecting the cancer type",strong("(1),"), "survival type",strong("(2),"), "and one or more proteins", strong("(3)."),"Selecting a cut-off", strong("(4)"), "will divide the samples into either 4 groups (quartile), 
                                3 groups (tertile) or into 2 groups (median or Optimal Cut-off) based on the protein selected. If multiple proteins are selected, the curve will be generated by taking the mean of 
                                the selected proteins and calculating the survival curve. Reported p-value is based on the log-rank test, while the hazard ratio is calculated for the high-versus-low comparison in 
                                a Cox Regression model (only available for the 2 group cut-off options). Graphs can be downloaded in either .pdf, .png, or .svg formats ", strong("(5)."),
                                hr(),
                                h1("Across Cancer"),
                                img(src="tab2_1.png", width=1066, height=618),
                                br(),
                                "Here the plot will report the hazard ratio and -log10(P-value) across all probes in the given cancer type based on the selected survival type and cut-off. This is performed by iterating through the Cox Regression model. 
                                Outputs are the hazard ratios for the single protein and -Log10P-value based on the Cox Regression calculations. Labeling for the volcano plot can be adjusted for poor prognosis markers", strong("(6)"), "or good prognosis markers", strong("(7)."), "
                                The volcano plot can also be transformed using the natural log of the hazard ratio, which can be used to spread out the x-axis for easy visualization", strong("(8)."),
                                br(),
                                br(),
                                img(src="tab2_2.png", width=1071, height=638),
                                br(),
                                "By clicking the 'Show Proportion' button in", strong("(9),"), "the high (red) and low (blue) sample proportions can be seen for the volcano plot. 
                                Downloading the data, the last three columns are used to generate the plot, the log10P equates to -log10(P-value) and hazard ratio inversed (HR_inverse), meaning the comparison of high versus low grouping.",
                                hr(),
                                h1("Across Protein"),
                                img(src="tab3_1.png", width=1042, height=547),
                                br(),
                                img(src="tab3_2.png", width=1013, height=538),
                                br(),
                                "Here each survival hazard ratios for a single protein can be viewed across all the cancer types. This feature uses the optimal cut-off feature for both overall survival and recurrence-free survival. Data from the volcano plot can also be visualized 
                                as a bar graph", strong("(10)."),
                                br(),
                                hr()
                                )))),
    tabPanel("About Us",
             fluidRow(
               column(12, offset = 1,
                  mainPanel(h3("About Us"), "This was a collaborative coding project between several medical students and post-Docs at the University of Iowa Carver College of Medicine. Principally, the effort was supported by the NIH F30 CA206255 fellowship 
                            and RO1 CA200673. Importantly, these results here are in whole or in part based upon data generated by the TCGA Research Network; you can find out more about the Network and each dataset at the", 
                            a("NCI Genomic Data Commons.", href="https://gdc.cancer.gov/"), "More information on the work conducted in the lab can be found at the", a("Zhang Lab Website.", href="https://zhang.lab.uiowa.edu/"),
                            br(),
                            h3("How to Cite:"), "Currently we are under the submission process. We will update with the citation when we have one.",
                            br(),
                            h3("Questions:"), "Code for this is available on", a("Github", href="https://github.com/ncborcherding"), "with more detailed explanation of methodology or in the citation above. If you still have questions, please email the link below.",
                            br(),
                            h3("Contact:"), "Please email with any questions, comments, or suggestions to", a("Nick Borcherding.", href="mailto:ncborch@gmail.com")
                            )
                  )
               )
            )
  )
)





