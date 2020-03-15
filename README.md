# shiny-DEG
A web application to analyze and visualize differentially expressed genes in RNA-seq.

## Information for the shiny-DEG App
Code canbe found on github:https://github.com/344968067/shiny-DEG  
To run this app locally on your machine,download R or Rstudio and run the following command once to set up the environment:  
install.packages(c("shiny","shinythemes","ggplot2","gplots","DESeq2","RColorBrewer","DT","pheatmap","reshape"))  

If you were ready for this packages, You may now run the shiny app with just one command in R:  
library("shiny")  
runApp("shiny-DEG")  
Or,  
shiny::runGitHub("shiny-DEG","344968067")  






