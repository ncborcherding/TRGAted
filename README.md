# TRGAted

[TRGAted](https://github.com/ncborcherding/TRGAted)
is an interactive web application that provides 
a graphical user interface combining survival information
and reverse-phase protein array data from 
[the Cancer Genome Atlas](https://cancergenome.nih.gov/).
Cloud-based TRGAted is run using the R Shiny Server
see the [TRGAted Application](https://nborcherding.shinyapps.io/TRGAted/).

## Citation

TRGAted is provided under a free-of-charge, open-source license (A-GPL3).
All we require is that you cite/attribute the following
in any work that benefits from this code or application.

### The App

If using the TRGAted app or derivative work, please cite our F1000research article: TRGAted: A web tool for survival analysis using protein data in the Cancer Genome Atlas. Link [Here.](https://f1000research.com/articles/7-1235/v2)


## Running TRGAted as a Local Session

While it is possible to host the server "back end" somewhere
so that users only need to point their web browser to a link,
it is also possible to launch both the back and front "ends" on your local machine.
The server back end will be an R session on your own machine,
while the front end is your web browser,
pointed to the appropriate local URL.

Make sure that you first have installed [the latest version of R](http://cran.r-project.org/). 

Download repository and open .R file in R Studio. launch TRGAted App. 

TRGAted relies on the following packages, which should be installed if running the code locally.

```r
install.packages("shiny") 
install.packages("ggplots2")
install.packages("broom")
install.packages("survMisc")
install.packages("survminer")
install.packages("shinythemes")
install.packages("DT")
install.packages("dplyr")
install.packages("ggrepel")
install.packages("reshape2")
```
