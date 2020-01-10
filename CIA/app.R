#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is a Shiny web application
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
# fda budesonide bioequivalence guidance and more
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(shiny)
library(nlme)
library(knitr)
library(VCA)
library(shinythemes)        # more funky looking apps
options(max.print=1000000)  # allows printing of long listings
fig.width <- 1200
fig.height <- 450
p1 <- function(x) {formatC(x, format="f", digits=1)}
p4 <- function(x) {formatC(x, format="f", digits=4)}
options(width=100)

#-----------------------------------------------
#A Quantitative Example - Carotid Stenosis Data
#see example 4.1 A Quantitative Example - Carotid Stenosis Data
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2856751/

# When the mean squared deviation (MSD) is used as the G function, the coefficient ??N varies between 0 and 1 [16], while psiR may exceed 1. For both coefficients, a value close to unity or above unity indicates an acceptable agreement. Haber and Barnhart [16] and Haber et al. [18] suggested that psi ??? 0.8 indicates acceptable agreement. Alternatively, one can consider agreement as acceptable if the confidence interval for the CIA includes 1.

#copy in data from Statistical Methods in Medical Research 2006; 15: 255-271 paper appendix
dat1 <- read.table(text="
1 100 100 100 100 100 100 100 100 100 2 73 67 77 69 60 72 59 62 62 3 0 0 3 49 0 0 0 12 0 4 17 0 0 0 0 0 0 0 0 5 27 40 0 15 31 2 38 29 19 6 31 29 34 45 51 3 38 42 41 7 18 15 17 21 23 99 57 44 0 8 68 59 100 39 56 33 61 67 77 9 0 0 0 0 0 0 47 42 99 10 18 0 0 37 0 12 25 0 13 11 0 0 0 0 0 0 0 0 0 12 72 65 69 100 99 100 61 60 58 13 27 27 28 42 18 99 0 0 99 14 61 53 51 100 100 71 45 44 53 15 100 100 100 100 100 100 100 100 100 16 0 0 0 0 0 0 16 99 15 17 30 31 22 25 49 9 51 40 44 18 10 27 37 21 38 12 40 28 26 19 29 33 37 50 30 25 37 45 44 20 57 60 67 100 100 100 65 63 65 21 100 100 100 0 100 0 13 17 100 22 28 17 0 22 15 0 49 0 0 23 0 34 0 21 32 0 0 18 0 24 30 34 41 55 53 99 60 46 49 25 46 40 62 56 36 0 63 70 64 26 66 61 71 78 76 55 82 83 100 27 51 43 68 100 28 0 70 0 0 28 23 27 39 45 60 26 24 44 0 29 0 0 0 74 100 100 100 100 100 30 83 100 100 70 99 100 53 100 100 31 0 0 0 26 0 0 12 0 0 32 4 0 0 25 0 0 24 34 23 33 100 100 100 100 100 100 100 100 100 34 60 60 75 45 47 15 62 99 67 35 21 19 0 0 0 0 4 0 0 36 70 79 100 23 36 0 51 42 0 37 6 18 0 7 4 99 99 99 99 38 41 52 49 44 63 0 80 73 0 39 56 45 66 100 53 0 13 23 0 40 5 0 0 30 45 31 70 34 55 41 53 39 54 63 50 31 75 72 55 42 47 56 51 7 99 0 56 58 0 43 75 84 85 100 100 100 100 100 100 44 100 100 100 100 100 100 100 100 100 45 0 100 100 44 0 99 17 100 0 46 58 61 55 46 58 49 71 99 99 47 0 0 0 30 22 0 1 24 99 48 0 0 0 0 0 0 0 0 0 49 5 0 0 0 39 0 0 0 0 50 12 25 9 26 54 0 100 0 100 51 0 0 8 0 0 0 0 0 0 52 28 29 0 100 100 99 100 84 100 53 0 0 0 0 0 0 0 0 0 54 19 31 51 53 100 100 77 83 100 55 33 74 50 69 72 100 79 100 100",
                   header=FALSE, stringsAsFactors=FALSE)


x<-as.vector(unlist(dat1))
d<-matrix(data = x, nrow = 55, ncol = 10, byrow = TRUE, dimnames = NULL)

d<-as.data.frame(d)



####################################################################################################
#Macro Name: CIA_macro                                                Last revision: Apr 4 2009

#Description: The function CIA estimates the CIA's and their standard errors for two observers, X and Y.  
#It can be used for binary data or numerical data with the MSD disagreement function.

#Parameters in the function CIA:

# 1. dataname: name of the data set
# 2. rep_X: number of replications for observer X (K), K > 1
# 3. rep_Y: number of replications for observer Y (L), L > 1 or L = 1
# 4. alpha: one minus the confidence coefficient



# Missing data should be denoted by blank--in R, they will be "NA"
# Way to remove the entire row with missing value
# http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/fSeries/html/na.html
# removeNA{fSeries}

# download fSeries before running this program
####################################################################################################

CIA <-function (dataname,rep_X,rep_Y,alpha){
     sample2<-na.omit(as.matrix(dataname)) #
    if (rep_Y == 1) {
        x <- sample2 [,2:(rep_X+1)]
        y <- sample2[,(rep_X+2)]
        id <- sample2[,1]
        wmean_x <- as.matrix(apply(x,1,mean),ncol=1)
        wstd_x <- as.matrix(apply(x,1,sd),ncol=1)
        est_wmsd_xx=2*wstd_x^2
        
        
        sum_diffsq <- c()
        est_wmsd_xy <- c()
        for (n in 1: length(id)){
            x_obs<-c()
            for (i in 1:rep_X){
                a<-matrix(x[n,i],ncol=1)
                x_temp<-a
                x_obs<-rbind(x_obs,x_temp)
            }
            y_obs<-matrix(rep(y[n],rep_X),ncol=1)
            compare_xy<-cbind(x_obs,y_obs)
            
            sum_temp <- sum(apply(compare_xy,1,function(t) (t[1]-t[2])**2))
            wmsd_temp<- sum_temp/(rep_X*rep_Y)
            sum_diffsq<-rbind(sum_diffsq, sum_temp)
            est_wmsd_xy <- rbind(est_wmsd_xy, wmsd_temp)
        }
        
        mean_G1 <- mean(est_wmsd_xx)
        mean_G3 <- mean(est_wmsd_xy)
        
        
        #calculate the SE and 95% CI
        
        simdata <- cbind(est_wmsd_xx, est_wmsd_xy)
        COV <- cov(simdata)
        var_G1 <- COV[1,1]
        var_G3 <- COV[2,2]
        Cov_G13 <- COV[1,2]
        
        
        N <- length(id) #number of valid observation
        
        A_r <- mean_G1
        B <- mean_G3
        
        Est_Psi_R <- A_r/B
        
        
        V_G_Mean_1 = var_G1/N
        V_G_Mean_3 = var_G3/N
        
        
        Cov_G_Mean_13 = Cov_G13/N
        
        Var_A_r = V_G_Mean_1
        Var_B = V_G_Mean_3
        
        Cov_A_r_B = Cov_G13/N
        
        Var_Est_Psi_R = (A_r^2/B**2) * ( Var_A_r/A_r^2 + Var_B/B^2 - 2*Cov_A_r_B/(A_r*B))
        
        SD_Est_Psi_R = sqrt(Var_Est_Psi_R)
        
        
        alpha2 = 1-0.05/2
        z = qnorm(alpha2)
        
        CI_Lower_Est_Psi_R =  Est_Psi_R - z*SD_Est_Psi_R
        CI_Upper_Est_Psi_R =  Est_Psi_R + z*SD_Est_Psi_R
        
        result_psi_R <- cbind(Est_Psi_R, SD_Est_Psi_R, CI_Lower_Est_Psi_R, CI_Upper_Est_Psi_R )
        print(result_psi_R)
    }
    
    
    if (rep_Y >1) {
        x <- sample2 [,2:(rep_X+1)]
        y <- sample2[,(rep_X+2):(rep_X+rep_Y+1)]
        id <- sample2[,1]
        wmean_x <- as.matrix(apply(x,1,mean),ncol=1)
        wmean_y <- as.matrix(apply(y,1,mean),ncol=1)
        wstd_x <- as.matrix(apply(x,1,sd),ncol=1)
        wstd_y <- as.matrix(apply(y,1,sd),ncol=1)
        est_wmsd_xx=2*wstd_x^2
        est_wmsd_yy=2*wstd_y^2
        
        
        sum_diffsq <- c()
        est_wmsd_xy <- c()
        for (n in 1: length(id)){
            x_obs<-c()
            for (i in 1:rep_X){
                a<-matrix(rep(x[n,i],rep_Y),ncol=1)
                x_temp<-a
                x_obs<-rbind(x_obs,x_temp)
            }
            y_obs<-matrix(rep(y[n,],rep_X),ncol=1)
            compare_xy<-cbind(x_obs,y_obs)
            
            sum_temp <- sum(apply(compare_xy,1,function(t) (t[1]-t[2])**2))
            wmsd_temp<- sum_temp/(rep_X*rep_Y)
            sum_diffsq<-rbind(sum_diffsq, sum_temp)
            est_wmsd_xy <- rbind(est_wmsd_xy, wmsd_temp)
        }
        
        mean_G1 <- mean(est_wmsd_xx)
        mean_G2 <- mean(est_wmsd_yy)
        mean_G3 <- mean(est_wmsd_xy)
        
        
        #calculate the SE and 95% CI
        
        simdata <- cbind(est_wmsd_xx, est_wmsd_yy, est_wmsd_xy)
        COV <- cov(simdata)
        var_G1 <- COV[1,1]
        var_G2 <- COV[2,2]
        var_G3 <- COV[3,3]
        Cov_G12 <- COV[1,2]
        Cov_G13 <- COV[1,3]
        Cov_G23 <- COV[2,3]
        
        N <- length(id) #number of valid observation
        
        A_n <- (mean_G1 + mean_G2)/2
        A_r <- mean_G1
        B <- mean_G3
        
        Est_Psi_N <- A_n/B
        Est_Psi_R <- A_r/B
        
        
        V_G_Mean_1 = var_G1/N
        V_G_Mean_2 = var_G2/N
        V_G_Mean_3 = var_G3/N
        
        Cov_G_Mean_12 = Cov_G12/N
        Cov_G_Mean_13 = Cov_G13/N
        Cov_G_Mean_23 = Cov_G23/N
        
        Var_A_n = (var_G1+ var_G2 + 2*Cov_G12)/(4*N)
        Var_A_n2 = (V_G_Mean_1 + V_G_Mean_2  + 2*Cov_G_Mean_12)/4
        
        Var_A_r = V_G_Mean_1
        Var_B = V_G_Mean_3
        
        Cov_A_n_B = (Cov_G13 + Cov_G23)/(2*N)
        Cov_A_r_B = Cov_G13/N
        
        Var_Est_Psi_N = (A_n^2/B**2) * ( Var_A_n/A_n^2 + Var_B/B^2 - 2*Cov_A_n_B/(A_n*B))
        Var_Est_Psi_R = (A_r^2/B**2) * ( Var_A_r/A_r^2 + Var_B/B^2 - 2*Cov_A_r_B/(A_r*B))
        
        SD_Est_Psi_N = sqrt(Var_Est_Psi_N)
        SD_Est_Psi_R = sqrt(Var_Est_Psi_R)
        
        
        alpha2 = 1-0.05/2
        z = qnorm(alpha2)
        CI_Lower_Est_Psi_N = Est_Psi_N - z*SD_Est_Psi_N
        CI_Upper_Est_Psi_N = Est_Psi_N + z*SD_Est_Psi_N
        CI_Lower_Est_Psi_R =  Est_Psi_R - z*SD_Est_Psi_R
        CI_Upper_Est_Psi_R =  Est_Psi_R + z*SD_Est_Psi_R
        
        if (Est_Psi_N > 1)  {
            CI_Lower_Est_Psi_N=.
            CI_Upper_Est_Psi_N=.
        }
        
        
        result_psi_N <- cbind(Est_Psi_N, SD_Est_Psi_N, CI_Lower_Est_Psi_N, CI_Upper_Est_Psi_N )
        print(result_psi_N)
        
        result_psi_R <- cbind(Est_Psi_R, SD_Est_Psi_R, CI_Lower_Est_Psi_R, CI_Upper_Est_Psi_R )
        print(result_psi_R)
        
    }
}

 
out <- CIA(dataname=d, rep_X = 3, rep_Y = 3, alpha=0.05)

d2<-d[,c(1,2,3,4, 8,9,10)]
out <- CIA(dataname=d2, rep_X = 3, rep_Y = 3, alpha=0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# end rep1 function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ui <- fluidPage(theme = shinytheme("journal"),
                
                shinyUI(pageWithSidebar(
                    
                    headerPanel("CIA Coefficient of Agreement"),
                    
                    sidebarPanel(
                        strong("Introduction"),
                        div(p("xx")),  
                        div(p("xxxxxxxxxxxxxxxx")),
                        
                        div(
                            
                            selectInput("Plot",
                                        strong("Select plot preference "),
                                        choices=c("VCA package plot" , "Base plot")),
                            
                            
                            selectInput("Model",
                                        strong("Select modelling preference "),
                                        choices=c("VCA package" , "nlme package")), 
                            
                            br(),br(),
                            div(strong("Tab 1a Plot the FDA example guidance data")),p("The example FDA dataset is plotted (always plot data), variance components are estimated for non subsetted data, just for interest. Note the replicates are linked, replicates labelled B, M and E have a specific meaning which cannot be ignored."),
                            
                            #  div(strong("Tab 1a Plot the FDA example guidance data")),p("The example FDA dataset is plotted (always plot data), variance components are estimated for non subsetted data, just for interest. Note the replicates are linked, replicates labelled B, M and E have a specific meaning and cannot be ignored."),
                            
                            div(strong("Tab 1b FDA Population Bioequivalence Statistical Analysis")),p("Reproduces the FDA analysis as presented in the guidance, notice differences in decimal places compared to the guidance [1]."),
                            div(strong("Tab 1c List the FDA guidance example data")),p("A simple listing of the FDA example data."),
                            tags$hr(),
                            actionButton("resample", "Simulate a new sample"),
                            br(),br(),
                            div(strong("Tab 2a Simulate data, plot and variance components analysis")),p("Does just that, plot a balanced design based on user inputs, use the slider inputs on this tab"),
                            div(strong("Tab 2b Population Bioequivalence Statistical Analysis - simulated data")),p("Performs the analysis approach based on guidance, but will adjust analysis if negative variance components are encountered. Presents a one way ANOVA on each product, performs PBS analysis and then concludes in a PASS or FAIL"),
                            div(strong("Tab 2c List the simulated data")),p("A simple listing of the simulated data."),
                            tags$hr(),
                            div(strong("Tab 3 Upload your own data for Population Bioequivalence Statistical Analysis")),p("User data can be loaded. PBE analysis takes place, the data plotted and listed. It does not have to be balanced, though it requires a balanced number of replicates within each product. The program will work out the design. Use at your own risk."),
                            tags$hr(),
                            div(strong("Tab 4 Explanation")),p("An explanation of the adjustments to the guidance that are performed if negative variance components are encountered. See the FDA guidance for an explanation of the general approach [1]."),
                            
                            br(),
                            actionButton(inputId='ab1', label="R code here", 
                                         icon = icon("th"), 
                                         onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/FDA-Bioequivalence/master/FDA_Bioequivalence/app.R', '_blank')"),
                            
                            br(),
                            br(),
                            div(p("References:")),  
                            
                            tags$a(href = "https://www.accessdata.fda.gov/drugsatfda_docs/psg/Budesonide_Inhalation_Sus_20929_RC_09-12.pdf", "[1] FDA Draft Guidance on Budesonide"),
                            div(p(" ")),
                            
                            tags$a(href = "https://www.fda.gov/media/70878/download", "[2] Statistical Information from the June 1999 Draft Guidance and Statistical Information for In Vitro Bioequivalence Data Posted on August 18, 1999"),
                            div(p(" "))
                            
                        ))
                    ,  
                    
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                    mainPanel(
                        
                        # https://shiny.rstudio.com/articles/upload.html
                        # https://stackoverflow.com/questions/44222796/shiny-with-multiple-tabs-and-different-sidebar-in-each-tab
                        
                        navbarPage(       
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                            tags$style(HTML(" 
                            .navbar-default .navbar-brand {color: cyan;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: lightgrey;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
                   ")), 
                            ##~~~~~~~~~~~~~~~~~~~~~end of section to add colour
                            
                            
                            #~~~~~~~~~~~~~~~~~~
                            tabPanel("1a Plot the FDA example guidance data", 
                                     
                                 #    div(plotOutput("reg.plot2", width=fig.width, height=fig.height)),  
                                     
                                     p(strong("Arithmetic mean presented above plot when VCA is used otherwise modelled mean
                                      (arithmetic mean and modelled mean will match with a balanced design)")) ,
                                     
                                  #   div( verbatimTextOutput("reg.summaryf"))
                                     
                                     
                            ) ,
                            #~~~~~~~~~~~~          
                            tabPanel("1b FDA Population Bioequivalence Statistical Analysis", 
                                     
                                     div( verbatimTextOutput("bioequivfda1")) #
                                     
                            ) ,
                            #~~~~~~~~~~~~~
                            tabPanel("1c List the FDA guidance example data", 
                                     
                                     div( verbatimTextOutput("bioequivfda2")), #
                                     
                                     p(strong("copy the above data to do your own analysis................."))
                                     
                            ),
                            #~~~~~~~~~~~~          
                            tabPanel("2a Simulate data, plot and variance components analysis", 
                                     
                                     
                                     div(strong("Select true population parameters for simulated data using the sliders below"),
                                         
                                         p("A 3 level nested data set is simulated. A plot of the raw data is generated. 
                    A model is fit to estimate the variance components. or the base plot, the blue dashed lines demark the top level groups. The green dashed lines demark the mid level groups.
                    The boxplots within the green lines demark the lowest level groups. Each boxplot presents the distribution of replicates. 
                    The middle, lowest and replicate numbers are varied randomly based on the slider ranges. The variance components are between blue 'top' groups,
                    between green 'mid' groups (within blue groups), within green 'mid' groups, within 'low' groups (replicates), a.k.a. repeatability. 
                    So we actually estimate 4 components counting the residual error. Product '1' is taken as the reference product.
                    Create a balanced design by reducing all the sliders to one value. 
                    The FDA guidance fits a one way ANOVA to each product's data with the independent variable the lowest level factor. 
                    The mid level 'batch' is not included in the analysis, however 'batch' is taken into account in the evaluation, the degrees of freedom are used in the analysis. 
                    The guidance does not discuss how to proceed in the case in which one or both between variance component(s) is estimated as negative, 
                    that's a bit Pepega, we show what to do in those scenarios, see notes tab for more information. You also have the choice of simulating a new sample. ")),
                                     ## new here -----------------------------------------------------------------------------------
                                     # https://stackoverflow.com/questions/47355043/multiple-column-layout-inside-a-tabpanel
                            #          sidebarLayout(
                            #              # Sidebar panel for inputs ----
                            #              sidebarPanel( 
                            #                  
                            #                  list(
                            #                      fluidRow(
                            #                          
                            #                          column(4,
                            #                                 sliderInput("top",
                            #                                             "Number of levels of top (product) component (demarked by blue or thick lines). Must be 2 for PBE analysis.",
                            #                                             min=2, max=100, step=1, value=2, ticks=FALSE)),
                            #                          column(4,
                            #                                 sliderInput("range1",
                            #                                             "Middle (batch) level: select no. of 'mid' groups within each top level group:", 
                            #                                             min = 2, max = 10, step=1, value = c(3,3), ticks=FALSE)),
                            #                          column(4,
                            #                                 sliderInput("range2", 
                            #                                             "Lower (sector) level: select\n no. of 'low' groups within each mid level group:",
                            #                                             min = 2, max = 10, step=1, value = c(10,10),ticks=FALSE))
                            #                      ),
                            #                      
                            #                      fluidRow(
                            #                          
                            #                          
                            #                          column(4,
                            #                                 sliderInput("replicates", 
                            #                                             "Select number of replicates nested within each boxplot",
                            #                                             min = 2, max = 50, step=1, value = c(3,3), ticks=FALSE)),
                            #                          column(4,
                            #                                 sliderInput("intercept",
                            #                                             "True intercept",
                            #                                             min=0, max=1000, step=.5, value=700, ticks=FALSE)),
                            #                          column(4,
                            #                                 sliderInput("a",
                            #                                             "True top level SD",
                            #                                             min=0, max=100, step=.5, value=75, ticks=FALSE))
                            #                      ),
                            #                      
                            #                      fluidRow(
                            #                          column(4,
                            #                                 sliderInput("b",
                            #                                             "True middle level SD",
                            #                                             min=0, max=100, step=.5, value=1, ticks=FALSE)),
                            #                          column(4,
                            #                                 sliderInput("c",
                            #                                             "True lower level SD",
                            #                                             min=0, max=100, step=.5, value=1, ticks=FALSE)),
                            #                          column(4,
                            #                                 sliderInput("d",
                            #                                             "True error",
                            #                                             min=0, max=100, step=.5, value=75, ticks=FALSE))
                            #                          
                            #                      )
                            #                  )
                            #                  
                            #                  ,   width = 12 ),
                            #              
                            #              # Main panel for displaying outputs ----
                            #              mainPanel(
                            #                  
                            #                  div(plotOutput("reg.plot", width=fig.width, height=fig.height)),
                            #                  
                            #                  p(strong("Arithmetic mean presented above plot when VCA is used otherwise modelled mean
                            #                 (arithmetic mean and modelled mean will match with a balanced design)")) ,
                            #                  
                            #                  div( verbatimTextOutput("reg.summary"))
                            #                  
                            #                  , width = 12 ))  # experiment with this
                            #          
                            ) ,
                            #~~~~~~~~~~~~
                            
                            #~~~~~~~~~~~~
                            tabPanel("2b Population Bioequivalence Statistical Analysis - simulated data", 
                                     
                                     div( verbatimTextOutput("bioequiv0")),
                                     
                                     p(strong(""))
                                     
                            ) ,
                            #~~~~~~~~~~~~~
                            tabPanel("2c List the simulated data", 
                                     
                                     div( verbatimTextOutput("summary2"))
                                     
                            ) ,
                            #~~~~~~~~~~~~~~~~~~
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   
                            tabPanel("3 Upload your own data for Population Bioequivalence Statistical Analysis", fluid = TRUE,
                                     
                                     p(("Upload your own data for determination of population bioequivalence. The top two radio button options are to help load,
                                 the bottom two options are to either show the top six rows of the data or show all the data, and secondly to toggle between PBE analysis or a plot of the data. 
                                 Ensure your data is balanced within product. Use at your own risk.")) ,
                                     
                                     p(("Example data sets (download either file and click 'Browse...' to locate and upload for PBE analysis):")) ,
                                     
                                     tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Three-level-nested-variance-components-analysis2/master/fda", "FDA example data set in guidance [1]"),
                                     div(p(" ")),
                                     
                                     tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Three-level-nested-variance-components-analysis2/master/fda%20B%20rep%20only", "The same FDA data set with only 1 rep per canister"),
                                     div(p(" ")),
                                     
                                     sidebarLayout(
                                         
                                         # Sidebar panel for inputs ----
                                         sidebarPanel(
                                             
                                             # Input: Select a file ----
                                             fileInput("file1", "Choose CSV File",
                                                       multiple = TRUE,
                                                       accept = c("text/csv",
                                                                  "text/comma-separated-values,text/plain",
                                                                  ".csv")),
                                             
                                             # Horizontal line ----
                                             tags$hr(),
                                             
                                             # Input: Checkbox if file has header ----
                                             checkboxInput("header", "Header", TRUE),
                                             
                                             # Input: Select separator ----
                                             radioButtons("sep", "Separator",
                                                          choices = c(Comma = ",",
                                                                      Semicolon = ";",
                                                                      Tab = "\t",
                                                                      Whitespace = ""),
                                                          selected = ""),
                                             
                                             # Input: Select quotes ----
                                             radioButtons("quote", "Quote",
                                                          choices = c(None = "",
                                                                      "Double Quote" = '"',
                                                                      "Single Quote" = "'"),
                                                          selected = ''),
                                             
                                             # Horizontal line ----
                                             tags$hr(),
                                             
                                             # Input: Select number of rows to display ----
                                             radioButtons("disp", "Display",
                                                          choices = c(Head = "head",
                                                                      All = "all"),
                                                          selected = "head"),
                                             
                                             # Horizontal line ----
                                             tags$hr(),
                                             
                                             # Input: Select number of rows to display ----
                                             radioButtons("what", "Output",
                                                          choices = c(Analysis = "Analysis",
                                                                      Plot = "plot"),
                                                          selected = "Analysis")
                                             
                                         ),
                                         
                                         # Main panel for displaying outputs ----
                                         mainPanel(
                                             
                                             # Output: Data file ----
                                             
                                             
                                             div(verbatimTextOutput("contents2")),
                                             plotOutput("plotx"),
                                             tags$hr(),
                                             
                                             tableOutput("contents") 
                                             
                                             
                                         ),
                                     )
                            ) #,
                           #  
                           #  tabPanel("4 Explanation", value=3, 
                           #           
                           #           HTML(" <strong>When the between variance component is estimated to be negative we apply the following adjustments.</strong>"),
                           #           
                           #           
                           #           p("Since three batches is not sufficient to reliably estimate the between batch component, the total
                           #              variances are estimated as the between canister variance of the 'super-batch' consisting of the three
                           #              batches combined [2]. It is not uncommon that the mean square between (MSB) is less than the mean square within (MSW) with a one way ANOVA. This results in a negative estimate 
                           #              for the between variance component. 
                           #   Thus concluding there is no additional variability due to the between variance component. In such cases the FDA PBE equations are adjusted. 
                           #   We have",HTML(" <em>m</em>"),"replicates, 
                           #   ",HTML(" <em>n</em>")," items per batch and",HTML(" <em>l</em>"),"is the no batches per product (test and reference). Refer to the FDA guidance document."),
                           #           br(),
                           #           
                           #           
                           #           HTML(" <strong>Impact to the total variances</strong>"),
                           #           withMathJax(
                           #               helpText('
                           # $${{\\sigma_R = }{\\sqrt{\\frac{MSB_R}{m} + \\frac{(m-1)MSW_R}{m}}}}\\!$$')),   
                           #           
                           #           
                           #           withMathJax(
                           #               helpText('This is equal to $${{}\\sigma_R ={\\sqrt{\\frac{MSB_R-MSW_R}{m} + MSW_R}}}\\!$$ In the event that $$MSB_R < MSW_R$$then $$MSB_R - MSW_R < 0$$and 
                           #               therefore $$\\sigma_R <  {\\sqrt{MSW_R}}$$
                           #    This means the total variance is less than the within variance component which cannot be.
                           #                        If this is encountered the total variance is set equal to the within variance component. 
                           #                        For either or both reference or test product if necessary. This is the first change from the guidance.')),
                           #           
                           #           withMathJax(
                           #               helpText('Therefore if $$MSB_R < MSW_R$$ then')),
                           #           
                           #           
                           #           
                           #           withMathJax(
                           #               helpText(" $$\\sigma_R =  {\\sqrt{MSW_R}}$$")),
                           #           
                           #           withMathJax(
                           #               helpText('and if $$MSB_T < MSW_T$$ then')),
                           #           
                           #           
                           #           withMathJax(
                           #               helpText(" $$\\sigma_T =  {\\sqrt{MSW_T}}$$")),
                           #           
                           #           
                           #           HTML(" <strong>Impact to Delta and HD</strong>"),
                           #           
                           #           withMathJax(
                           #               helpText('$$ \\hat{\\Delta} = \\bar{y}_T - \\bar{y}_R$$')),
                           #           
                           #           br(),
                           #           
                           #           withMathJax(
                           #               helpText("We have")),
                           #           
                           #           withMathJax(
                           #               helpText('$$Var(\\bar{y}_T) = \\frac{\\sigma^2_B}{n.l}  +  \\frac{\\sigma^2_W}{n.l.m}  = \\frac{m\\sigma^2_B + \\sigma^2_W}{n.l.m} = \\frac{MSB_T}{n.l.m} $$')),
                           #           
                           #           withMathJax(
                           #               helpText("with degrees of freedom $$n_T.l_T-1$$")),
                           #           br(),
                           #           withMathJax(
                           #               helpText("We have")),
                           #           withMathJax(
                           #               helpText('$$Var(\\bar{y}_R) = \\frac{\\sigma^2_B}{n.l}  +  \\frac{\\sigma^2_W}{n.l.m}  = \\frac{m\\sigma^2_B + \\sigma^2_W}{n.l.m} = \\frac{MSB_R}{n.l.m} $$')),
                           #           
                           #           withMathJax(
                           #               helpText("with degrees of freedom $$n_R.l_R-1$$")),
                           #           
                           #           withMathJax(
                           #               helpText("The variances of the difference is equal to the sum of the variances, so")),
                           #           
                           #           withMathJax(
                           #               helpText('$$ Var\\hat{\\Delta} =  \\frac{MSB_T}{n.l.m} +  \\frac{MSB_R}{n.l.m} $$')),
                           #           
                           #           withMathJax(
                           #               helpText("with degrees of freedom $$(n_R.l_R-1) + (n_T.l_T-1) = (n_R.l_R + n_T.l_T-2)$$")),
                           #           
                           #           withMathJax(
                           #               helpText("When there is no between variance component for either test or reference then:")),
                           #           
                           #           withMathJax(
                           #               helpText('$$Var(\\bar{y}_T) =  \\frac{MSW_T}{n.l.m} $$')),
                           #           
                           #           withMathJax(
                           #               helpText("with degrees of freedom $$n_T.l_T.(m-1)$$")),
                           #           
                           #           withMathJax(
                           #               helpText("and")),
                           #           
                           #           withMathJax(
                           #               helpText('$$Var(\\bar{y}_R) =  \\frac{MSW_R}{n.l.m} $$')),                                    
                           #           
                           #           withMathJax(
                           #               helpText("with degrees of freedom $$n_R.l_R.(m-1)$$")),
                           #           
                           #           withMathJax(
                           #               helpText("and so depending on if one or both products has a single variance component, there are four scenarios for the calculation of HD including the one in the guidance:")),
                           #           
                           #           withMathJax(
                           #               helpText("I) Both have non negative variance components:")),
                           #           
                           #           withMathJax(
                           #               helpText('$$H_D = \\left(\\lvert\\hat{\\Delta}\\rvert + 
                           #              t_{1-\\alpha,n_T.l_T-1 + n_R.l_R-1)}  
                           #                (\\frac{MSB_T}{n.l.m} + \\frac{MSB_R}{n.l.m})^.5\\right)^2  $$')),
                           #           
                           #           withMathJax(
                           #               helpText("II) Both have negative variance components:")),
                           #           
                           #           withMathJax(
                           #               helpText('$$H_D = \\left(\\lvert\\hat{\\Delta}\\rvert + 
                           #              t_{1-\\alpha,n_T.l_T.(m-1) + n_R.l_R.(m-1)}  
                           #                (\\frac{MSW_T}{n.l.m} + \\frac{MSW_R}{n.l.m})^.5\\right)^2  $$')),
                           #           
                           #           withMathJax(
                           #               helpText("III) Test only has negative variance component:")),
                           #           
                           #           withMathJax(
                           #               helpText('$$H_D = \\left(\\lvert\\hat{\\Delta}\\rvert + 
                           #              t_{1-\\alpha,n_T.l_T.(m-1) + n_R.l_R-1)}  
                           #                (\\frac{MSW_T}{n.l.m} + \\frac{MSB_R}{n.l.m})^.5\\right)^2  $$')),
                           #           
                           #           withMathJax(
                           #               helpText("IV) Reference only has negative variance component:")),
                           #           
                           #           withMathJax(
                           #               helpText('$$H_D = \\left(\\lvert\\hat{\\Delta}\\rvert + 
                           #               t_{1-\\alpha,n_T.l_T-1  + n_R.l_R.(m-1)}  
                           #                (\\frac{MSB_T}{n.l.m} + \\frac{MSW_R}{n.l.m})^.5\\right)^2  $$')),
                           #           
                           #           HTML(" <strong>Impact on FDA parameters E1 and E2</strong>"),
                           #           
                           #           withMathJax(
                           #               helpText("When $$MSB_T >= MSW_T$$ we have")),
                           #           
                           #           withMathJax(
                           #               helpText("$$E1 + E2 = \\frac{MSB_T}{m} + \\frac{(m-1) MSW_T }{m } = \\frac{MSB_T - MSW_T }{m}  +  MSW_T = \\hat{\\sigma^2_B} + \\hat{\\sigma^2_W}$$ 
                           #                       which is the total variance of test products. 
                           #                       If we encounter a negative between variance component for the test product we have shown above the total variance is estimated by the mean squares within. 
                           #                       
                           #                       So let $$E1 = 0, E2 = MSW_T$$ and H2 and U2 are unaltered.")), 
                           #           
                           #           HTML(" <strong>Impact on FDA parameters E3s and E4s, reference scaling</strong>"),
                           #           
                           #           withMathJax(
                           #               helpText("When $$MSB_ r >= MSW_R$$ we have")),
                           #           
                           #           withMathJax(
                           #               helpText("$$E3s + E4s = -(1+\\theta_p)\\frac{MSB_R}{m} -(1+\\theta_p) \\frac{(m-1) MSW_R }{m } = -(1+\\theta_p)\\left(\\frac{MSB_R - MSW_R }{m}  +
                           #              MSW_R\\right) = -(1+\\theta_p)\\left(\\hat{\\sigma^2_B} + \\hat{\\sigma^2_W}\\right)$$ 
                           #              
                           #                       If we encounter a negative between variance component for the reference product we have shown above the total variance is estimated by the mean squares within. 
                           #                       
                           #                       So let $$E3s = 0, E4s = -(1+\\theta_p)MSW_R$$ and H4s and U4s are unaltered.")), 
                           #           
                           #           HTML(" <strong>Impact on FDA parameters E3c and E4c, constant scaling</strong>"),
                           #           
                           #           withMathJax(
                           #               helpText("When $$MSB_ r >= MSW_R$$ we have")),
                           #           
                           #           withMathJax(
                           #               helpText("$$E3c + E4c = -\\frac{MSB_R}{m} + \\frac{(m-1) MSW_R }{m } = -\\left(\\frac{MSB_R - MSW_R }{m}  +  MSW_R\\right) = -\\left(\\hat{\\sigma^2_B} +
                           #              \\hat{\\sigma^2_W}\\right)$$ 
                           #                         
                           #                       If we encounter a negative between variance component for the test product we have shown above the total variance is estimated by the mean squares within. 
                           #                       
                           #                       So let $$E3c = 0, E4c = -MSW_R$$ and H4c and U4c are unaltered.")), 
                           #           
                           #           br(),
                           #           br(),
                           #           br(),
                           #           br())
                            
                            
                            
                            
                        )
                    )
                ) 
                )
)
#~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~


server <- shinyServer(function(input, output) {
    
    # --------------------------------------------------------------------------
    # This is where a new sample is instigated, only random noise is required to be generated
    random.sample <- reactive({
        
        # Dummy line to trigger off button-press
        foo <- input$resample
        
        x1 <- input$range1[1]
        x2 <- input$range1[2]
        x3 <- input$range2[1]
        x4 <- input$range2[2]
        x5 <- input$replicates[1]
        x6 <- input$replicates[2]
        
        top <-  input$top
        
        # seems that I need to use both c(x1,x2) c(x1:x2) so sample function works correctly
        
        if (x1==x2) {
            
            middle <-  sample(c(x1,x2),   top, replace=TRUE)    # ditto groups in each top level 6
            
        } else {
            
            middle <-  sample(c(x1:x2),   top, replace=TRUE)    # ditto groups in each top level 6
        }
        
        
        if (x3==x4) {
            
            lower <-   sample(c(x3,x4),   sum(middle), replace=TRUE )
            
        } else {
            
            lower <-   sample(c(x3:x4),   sum(middle), replace=TRUE )
            
        }
        
        if (x5==x6) {
            
            replicates <-  sample(c(x5,x6),   sum(lower), replace=TRUE )
            
        } else {
            
            replicates <-  sample(c(x5:x6),   sum(lower), replace=TRUE )
            
        }
        
        N <- sum(replicates)
        
        a=input$a
        b=input$b
        c=input$c
        d=input$intercept
        residual <- input$d
        
        # random effects
        top.r <-    rnorm(top,          d,                a)    
        middle.r <- rnorm(sum(middle),  0,                b)    
        lower.r <-  rnorm(sum(lower),   0,                c)    
        
        # ids
        lower.id <- rep(seq_len(sum(lower)), replicates )       
        middle.id <- cut(lower.id, c(0,cumsum(lower)),  labels=FALSE)
        top.id   <- cut(middle.id, c(0,cumsum(middle)), labels=FALSE)
        
        return(list(top.r=top.r, middle.r=middle.r, lower.r=lower.r, intercept=d, residual=residual,
                    middle=middle, top=top, lower=lower, replicates=replicates,
                    N=N, top.id=top.id, middle.id=middle.id, lower.id=lower.id, 
                    a=a, b=b, c=c))
        
    }) 
    
    # --------------------------------------------------------------------------
    # Set up the dataset based on the inputs 
    make.regression <- reactive({
        
        sample <- random.sample()
        
        top <-        sample$top
        middle <-     sample$middle
        lower <-      sample$lower
        replicates <- sample$replicates
        
        # random effects
        top.r <-    sample$top.r 
        middle.r <- sample$middle.r 
        lower.r <-  sample$lower.r 
        
        # ids
        lower.id <- sample$lower.id
        middle.id <- sample$middle.id
        top.id   <- sample$top.id 
        
        Data <- data.frame( top=top.id, mid=middle.id, low=lower.id,
                            y= rnorm( sum(replicates), top.r[top.id] + 
                                          middle.r[middle.id] + lower.r[lower.id], sample$residual ) )
        
        df <- as.data.frame(Data)
        return(list(df=df , middle.id=middle.id, top.id=top.id)) 
        
    })  
    # --------------------------------------------------------------------------
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Fit the specified regression model on simulated data  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fit.regression <- reactive({
        
        data <- make.regression()
        
        df <- data$df
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Conditionally fit the model
        
        if (input$Model == "nlme package") {
            
            fit.res <-  
                tryCatch(intervals(lme(y ~ 1, random = ~1 |  top/mid/low , data=df, method="REML")), 
                         error=function(e) e)
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ###http://stackoverflow.com/questions/8093914
            ###/skip-to-next-value-of-loop-upon-error-in-r-trycatch
            
            if (!inherits(fit.res, "error")) {
                
                modelint <- fit.res
                
                emu      <-p1(modelint[['fixed']][2][[1]])  
                etop     <-p1(modelint[['reStruct']][['top']][2][[1]])
                eday     <-p1(modelint[['reStruct']][['mid']][2][[1]])
                erun     <-p1(modelint[['reStruct']][['low']][2][[1]])
                esigma   <-p1(modelint[['sigma']][2][[1]])
                
            } else  {
                
                fit.res <- NULL
                
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        } else {                   
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            o <- fit.res<- tryCatch(fitVCA(y~top/mid/low, df, "reml"), 
                                    error=function(e) e) 
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            if (!inherits(fit.res, "error")) {
                
                fit.res <- VCAinference(fit.res, ci.method = "sas")
                
                x <- as.matrix(o)
                features <- attributes(x)
                
                emu      <- p1(features$Mean) 
                
                o <- as.matrix(o)
                etop     <-p1(o["top","SD"])
                eday     <-p1(o["top:mid","SD"])
                erun     <-p1(o["top:mid:low","SD"])
                esigma   <-p1(o["error","SD"])
                
            } else  {
                
                fit.res <- NULL
                
            }
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        # Get the model summary
        if (is.null(fit.res)) {
            
            fit.summary <- NULL
            
        } else {
            
            fit.summary <-  (fit.res)
        }
        
        return(list(emu=emu, etop=etop, eday=eday, erun=erun,
                    esigma=esigma, fit.res=fit.res, fit.summary=fit.summary)) 
        
    })     
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plot a scatter of the simulated data  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$reg.plot <- renderPlot({         
        
        # Get the current regression data
        data1 <- make.regression()
        
        d1 <- data1$df
        
        # Conditionally plot
        if (input$Plot == "Base plot") {
            
            #base plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # middle groups for plot
            xxx<- cumsum(as.vector(table(data1$middle.id)))
            xxx <-d1$low[xxx]
            
            # top groups for plot
            yyy<- cumsum(as.vector(table(data1$top.id)))
            yyy <-d1$low[yyy]
            
            plot( y ~ factor(low), data=d1 , 
                  main=paste("Variability Chart. Truth (estimate): intercept ",input$intercept,"(",fit.regression()$emu,"), top level sd=",
                             input$a,"(",fit.regression()$etop,")", ",\n middle level sd=",
                             input$b ,"(",fit.regression()$eday,"), lowest level sd=",
                             input$c, "(",fit.regression()$erun,") & random error sd=", 
                             input$d,"(",fit.regression()$esigma,")"),
                  main="lowest level grouped", xlab="lowest level groups")
            
            abline( v=c(0,xxx)+0.5, lty=2, col='green' )
            abline( v=c(0,yyy)+0.5, lty=2, col='blue' )
            
            
        } else {
            
            #VCA plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            varPlot(y~top/mid/low, d1, 
                    BG=list(var="top", 
                            col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
                                  "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"), 
                            col.table=TRUE), 
                    VLine=list(var=c("top", "mid"), 
                               col=c("black", "mediumseagreen"), lwd=c(2,1), 
                               col.table=c(TRUE,TRUE)), 
                    JoinLevels=list(var="low", col=c("lightblue", "cyan", "yellow"), 
                                    lwd=c(2,2,2)), 
                    MeanLine=list(var="top", col="blue", lwd=2),
                    Title=list(main=paste("Variability Chart. Truth (estimate): intercept ",input$intercept,"(",fit.regression()$emu,"), top level sd=",
                                          input$a,"(",fit.regression()$etop,")", ",\n middle level sd=",
                                          input$b ,"(",fit.regression()$eday,"), lowest level sd=",
                                          input$c, "(",fit.regression()$erun,") & random error sd=", 
                                          input$d,"(",fit.regression()$esigma,")")),
                    
                    # MeanLine=list(var="mid", col="pink", lwd=2),
                    Points=list(pch=list(var="mid", pch=c(21, 22, 24)), 
                                bg =list(var="mid", bg=c("lightblue", "cyan", "yellow")), 
                                cex=1.25))    
            
            
        }
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # analyse fda variance components for plot title
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    vc.fda <-  reactive({
        
        df <- fda.d
        
        #df$y <- log(df$y)  # log the fda data
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Conditionally fit the model
        
        if (input$Model == "nlme package") {
            
            fit.res <-  
                tryCatch(intervals(lme(y ~ 1, random = ~1 |  PRODUCT/BATCH/SECTOR , data=df, method="REML")), 
                         error=function(e) e)
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ###http://stackoverflow.com/questions/8093914
            ###/skip-to-next-value-of-loop-upon-error-in-r-trycatch
            
            if (!inherits(fit.res, "error")) {
                
                modelint <- fit.res
                
                emu      <-p4(modelint[['fixed']][2][[1]])  
                etop     <-p4(modelint[['reStruct']][['PRODUCT']][2][[1]])
                eday     <-p4(modelint[['reStruct']][['BATCH']][2][[1]])
                erun     <-p4(modelint[['reStruct']][['SECTOR']][2][[1]])
                esigma   <-p4(modelint[['sigma']][2][[1]])
                
            } else  {
                
                fit.res <- NULL
                
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        } else {                   
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            o <- fit.res<- tryCatch(fitVCA(y~PRODUCT/BATCH/SECTOR, df, "reml"), 
                                    error=function(e) e) 
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            if (!inherits(fit.res, "error")) {
                
                fit.res <- VCAinference(fit.res, ci.method = "sas")
                
                x <- as.matrix(o)
                features <- attributes(x)
                
                emu      <- p4(features$Mean) 
                
                o <- as.matrix(o)
                etop     <-p4(o["PRODUCT","SD"])
                eday     <-p4(o["PRODUCT:BATCH","SD"])
                erun     <-p4(o["PRODUCT:BATCH:SECTOR","SD"])
                esigma   <-p4(o["error","SD"])
                
            } else  {
                
                fit.res <- NULL
                
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        # Get the model summary
        if (is.null(fit.res)) {
            
            fit.summary <- NULL
            
        } else {
            
            fit.summary <-  (fit.res)
        }
        
        return(list(emu=emu, etop=etop, eday=eday, erun=erun,
                    esigma=esigma, fit.res=fit.res, fit.summary=fit.summary))
        
    })     
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # plot fda data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$reg.plot2 <- renderPlot({         
        
        # Get the fda  data
        
        d1 <- fda.d
        
        
        # Conditionally plot
        if (input$Plot == "Base plot") {
            
            #base plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # middle groups for plot
            xxx<- cumsum(as.vector(table(d1$BATCH)))
            xxx <-d1$SECTOR[xxx]
            
            # top groups for plot
            yyy<- cumsum(as.vector(table(d1$PRODUCT)))
            yyy <-d1$SECTOR[yyy]
            
            plot( y ~ factor(SECTOR), data=d1 , 
                  main=paste("Variability Chart of FDA data. (Estimate in brackets): intercept (",vc.fda()$emu,"), top level sd=(",vc.fda()$etop,")", 
                             ",\n middle level sd=(",vc.fda()$eday,"), lowest level sd=(",vc.fda()$erun,") & random error sd=("
                             ,vc.fda()$esigma,")"),
                  main="lowest level grouped", xlab="lowest level groups", ylab="Guidance response (log y)")
            
            abline( v=c(0,xxx)+0.5, lty=2, col='green' )
            abline( v=c(0,yyy)+0.5, lty=2, col='blue' )
            
            
        } else {
            
            # VCA plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            varPlot(y~PRODUCT/BATCH/SECTOR/REP, d1, 
                    YLabel = list(text =
                                      "Guidance response (log y)", side = 2, line = 3.5, cex = 1.5),
                    BG=list(var="PRODUCT", 
                            col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
                                  "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"), 
                            col.table=TRUE), 
                    VLine=list(var=c("PRODUCT", "BATCH"), 
                               col=c("black", "mediumseagreen"), lwd=c(2,1), 
                               col.table=c(TRUE,TRUE)), 
                    JoinLevels=list(var="SECTOR", col=c("lightblue", "cyan", "yellow"), 
                                    lwd=c(2,2,2)), 
                    MeanLine=list(var="PRODUCT", col="blue", lwd=2),
                    Title=list(main=paste("Variability Chart of FDA data. (Estimate in brackets): intercept (",vc.fda()$emu,"), top level sd=(",vc.fda()$etop,")", 
                                          ",\n middle level sd=(",vc.fda()$eday,"), lowest level sd=(",vc.fda()$erun,") & random error sd=("
                                          ,vc.fda()$esigma,")")),
                    
                    # MeanLine=list(var="mid", col="pink", lwd=2),
                    Points=list(pch=list(var="BATCH", pch=c(21, 22, 24)), 
                                bg =list(var="BATCH", bg=c("lightblue", "cyan", "yellow")), 
                                cex=1.25))    
            
        }
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # execute PBE data analysis on simulation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bioequivf <- reactive({
        
        data1 <- make.regression()
        
        foo <- data1$df
        foo$top <-    ifelse(foo$top %in% 1, "REF", "TEST")  # arbitrarily make 1 the reference product
        foo$top <-    factor(foo$top)
        foo$SECTOR <- foo$low
        foo$BATCH <-  foo$mid
        
        A <- unique(input$range1)*unique(input$range2)
        B <- unique(input$replicates)
        
        res <- bioequiv(foo1=foo, nrXlr=A, mr= B, 
                        ntXlt=A, mt= B,
                        response="y",indep="low", split="top", ref="REF", test="TEST")
        
        return(list(res=res))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # execute fda analysis
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bioequivfdax <- reactive({
        
        
        foo8 <- fda.d
        
        foo8$y <- exp(foo8$y)  # the bioequiv function logs the data, but the fda data is already logged
        
        res.fda <- bioequiv(foo1=foo8 , nrXlr=10*3, mr= 3, 
                            ntXlt=10*3, mt= 3,
                            response="y",indep="SECTOR", split="PRODUCT", ref="REF", test="TEST")
        
        return(res.fda)
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # collect fda data for listing
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    
    
    fda.data <- reactive({
        
        foo99 <- fda2
        
        return(foo99)
        
    })
    
    #---------------------------------------------------------------------------
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # fda data var comp output
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$reg.summaryf <- renderPrint({
        
        summary <- vc.fda()$fit.summary
        
        if (!is.null(summary)) {
            
            return(vc.fda()$fit.summary)
            
        }
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # simulated data var comp output
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$reg.summary <- renderPrint({
        
        summary <- fit.regression()$fit.summary
        
        if (!is.null(summary)) {
            
            return(fit.regression()$fit.summary)
            
        }
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # lsting of simulated data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$summary2 <- renderPrint({
        
        return(make.regression()$df)
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # collect simulation analysis
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$bioequiv0 <- renderPrint({
        
        return(bioequivf()$res)
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # collect the fda data listing
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$bioequivfda2 <- renderPrint({
        
        return(fda2)
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # collect the fda analysis
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$bioequivfda1 <- renderPrint({
        
        return(bioequivfdax()$res.fda)
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # loading in user data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        df<- NULL
        req(input$file1)
        df <- read.csv(input$file1$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
        
        df <- as.data.frame(df)
        
        if(input$disp == "head") {
            
            return(head(df))
        }
        else {
            
            return(df)
        } 
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PBE analysis of user uploaded data 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$contents2 <- renderPrint({
        
        if(input$what == "Analysis"){
            
            df<-NULL
            req(input$file1)
            df <- read.csv(input$file1$datapath,
                           header = input$header,
                           sep = input$sep,
                           quote = input$quote)
            
            df<- as.data.frame(df)
            
            df$y <- exp(df$y) 
            
            REF <-  df[df$PRODUCT %in% "REF",]
            TEST <- df[df$PRODUCT %in% "TEST",]
            
            #work out the design
            mr <- unique(with(REF, table(SECTOR)))
            mt <- unique(with(TEST, table(SECTOR))) 
            
            x <- REF[,c(1,2)]
            x <- unique(x)
            nrXlr <- sum(as.vector(table(x$BATCH)))
            
            
            x <- TEST[,c(1,2)]
            x <- unique(x)
            ntXlt <- sum(as.vector(table(x$BATCH)))
            
            ##new 
            if(mr ==1 & mt==1) {
                
                
                user.analysis<- bioequiv.1.rep(foo1=df , 
                                               nrXlr=nrXlr, mr=mr, 
                                               ntXlt=ntXlt, mt=mt,
                                               response="y",indep="SECTOR", split="PRODUCT", ref="REF", test="TEST")
                
            } else {
                
                user.analysis<- bioequiv(foo1=df , 
                                         nrXlr=nrXlr, mr= mr, 
                                         ntXlt=ntXlt, mt= mr,
                                         response="y",indep="SECTOR", split="PRODUCT", ref="REF", test="TEST")
                
            }
        }})
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # analyse on user data, variance components for plot title
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vc.user <-  reactive({
        
        req(input$file1)
        df <- read.csv(input$file1$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
        
        df<- as.data.frame(df)
        
        REF <-  df[df$PRODUCT %in% "REF",]
        TEST <- df[df$PRODUCT %in% "TEST",]
        
        #work out the design
        mr <- unique(with(REF, table(SECTOR)))
        mt <- unique(with(TEST, table(SECTOR))) 
        
        x <- REF[,c(1,2)]
        x <- unique(x)
        nrXlr <- sum(as.vector(table(x$BATCH)))
        
        x <- TEST[,c(1,2)]
        x <- unique(x)
        ntXlt <- sum(as.vector(table(x$BATCH)))
        
        ##new 
        if(mr ==1 & mt==1) {
            
            #   df$y <- log(df$y)  # log the fda data
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Conditionally fit the model
            
            if (input$Model == "nlme package") {
                
                fit.res <-  
                    tryCatch(intervals(lme(y ~ 1, random = ~1 |  PRODUCT/BATCH , data=df, method="REML")), 
                             error=function(e) e)
                
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ###http://stackoverflow.com/questions/8093914
                ###/skip-to-next-value-of-loop-upon-error-in-r-trycatch
                
                if (!inherits(fit.res, "error")) {
                    
                    modelint <- fit.res
                    
                    emu      <-p4(modelint[['fixed']][2][[1]])  
                    etop     <-p4(modelint[['reStruct']][['PRODUCT']][2][[1]])
                    eday     <-p4(modelint[['reStruct']][['BATCH']][2][[1]])
                    # erun     <-p4(modelint[['reStruct']][['SECTOR']][2][[1]])
                    esigma   <-p4(modelint[['sigma']][2][[1]])
                    
                } else  {
                    
                    fit.res <- NULL
                    
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            } else {                   
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                o <- fit.res<- tryCatch(fitVCA(y~PRODUCT/BATCH, df, "reml"), 
                                        error=function(e) e) 
                
                
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                if (!inherits(fit.res, "error")) {
                    
                    fit.res <- VCAinference(fit.res, ci.method = "sas")
                    
                    x <- as.matrix(o)
                    features <- attributes(x)
                    
                    emu      <- p4(features$Mean) 
                    
                    o <- as.matrix(o)
                    etop     <-p4(o["PRODUCT","SD"])
                    eday     <-p4(o["PRODUCT:BATCH","SD"])
                    #erun     <-p4(o["PRODUCT:BATCH:SECTOR","SD"])
                    esigma   <-p4(o["error","SD"])
                    
                } else  {
                    
                    fit.res <- NULL
                    
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # MORE THAN 1 REP   
            
        } else {
            
            #   df$y <- log(df$y)  # log the fda data
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Conditionally fit the model
            
            if (input$Model == "nlme package") {
                
                fit.res <-  
                    tryCatch(intervals(lme(y ~ 1, random = ~1 |  PRODUCT/BATCH/SECTOR , data=df, method="REML")), 
                             error=function(e) e)
                
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ###http://stackoverflow.com/questions/8093914
                ###/skip-to-next-value-of-loop-upon-error-in-r-trycatch
                
                if (!inherits(fit.res, "error")) {
                    
                    modelint <- fit.res
                    
                    emu      <-p4(modelint[['fixed']][2][[1]])  
                    etop     <-p4(modelint[['reStruct']][['PRODUCT']][2][[1]])
                    eday     <-p4(modelint[['reStruct']][['BATCH']][2][[1]])
                    erun     <-p4(modelint[['reStruct']][['SECTOR']][2][[1]])
                    esigma   <-p4(modelint[['sigma']][2][[1]])
                    
                } else  {
                    
                    fit.res <- NULL
                    
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            } else {                   
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                o <- fit.res<- tryCatch(fitVCA(y~PRODUCT/BATCH/SECTOR, df, "reml"), 
                                        error=function(e) e) 
                
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                if (!inherits(fit.res, "error")) {
                    
                    fit.res <- VCAinference(fit.res, ci.method = "sas")
                    
                    x <- as.matrix(o)
                    features <- attributes(x)
                    
                    emu      <- p4(features$Mean) 
                    
                    o <- as.matrix(o)
                    etop     <-p4(o["PRODUCT","SD"])
                    eday     <-p4(o["PRODUCT:BATCH","SD"])
                    erun     <-p4(o["PRODUCT:BATCH:SECTOR","SD"])
                    esigma   <-p4(o["error","SD"])
                    
                } else  {
                    
                    fit.res <- NULL
                    
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        }  # end
        
        ##new 
        if(mr ==1 & mt==1) { 
            
            # Get the model summary
            if (is.null(fit.res)) {
                
                fit.summary <- NULL
                
            } else {
                
                fit.summary <-  (fit.res)
            }
            
            return(list(emu=emu, etop=etop, eday=eday, #erun=NA,
                        esigma=esigma, fit.res=fit.res, fit.summary=fit.summary))
            
        } else {
            
            if (is.null(fit.res)) {
                
                fit.summary <- NULL
                
            } else {
                
                fit.summary <-  (fit.res)
            }
            
            return(list(emu=emu, etop=etop, eday=eday, erun=erun,
                        esigma=esigma, fit.res=fit.res, fit.summary=fit.summary))
            
        }  
        
    })     
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plotting user uploaded data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plotx <- renderPlot({
        if(input$what == "plot"){
            
            df<-NULL
            req(input$file1)
            df <- read.csv(input$file1$datapath,
                           header = input$header,
                           sep = input$sep,
                           quote = input$quote)
            
            df<- as.data.frame(df)
            
            REF <-  df[df$PRODUCT %in% "REF",]
            TEST <- df[df$PRODUCT %in% "TEST",]
            
            #work out the design
            mr <- unique(with(REF, table(SECTOR)))
            mt <- unique(with(TEST, table(SECTOR))) 
            
            x <- REF[,c(1,2)]
            x <- unique(x)
            nrXlr <- sum(as.vector(table(x$BATCH)))
            
            x <- TEST[,c(1,2)]
            x <- unique(x)
            ntXlt <- sum(as.vector(table(x$BATCH)))
            
            ##new 
            if(mr ==1 & mt==1) {
                # df$y <- log(df$y)  # log the fda data
                
                varPlot(y~PRODUCT/BATCH/SECTOR/REP, df, 
                        BG=list(var="PRODUCT", 
                                col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
                                      "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"), 
                                col.table=TRUE), 
                        VLine=list(var=c("PRODUCT", "BATCH"), 
                                   col=c("black", "mediumseagreen"), lwd=c(2,1), 
                                   col.table=c(TRUE,TRUE)), 
                        JoinLevels=list(var="SECTOR", col=c("lightblue", "cyan", "yellow"), 
                                        lwd=c(2,2,2)), 
                        MeanLine=list(var="PRODUCT", col="blue", lwd=2),
                        Title=list(main=paste("Variability Chart. Estimate: intercept (", vc.user()$emu,"), top level sd= (",vc.user()$etop,")"
                                              , ",\n lower level sd= (",vc.user()$eday,"),  random error sd= (",vc.user()$esigma,")")),
                        
                        # MeanLine=list(var="mid", col="pink", lwd=2),
                        Points=list(pch=list(var="BATCH", pch=c(21, 22, 24)), 
                                    bg =list(var="BATCH", bg=c("lightblue", "cyan", "yellow")), 
                                    cex=1.25)) 
            }  else {   
                
                varPlot(y~PRODUCT/BATCH/SECTOR/REP, df, 
                        BG=list(var="PRODUCT", 
                                col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
                                      "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"), 
                                col.table=TRUE), 
                        VLine=list(var=c("PRODUCT", "BATCH"), 
                                   col=c("black", "mediumseagreen"), lwd=c(2,1), 
                                   col.table=c(TRUE,TRUE)), 
                        JoinLevels=list(var="SECTOR", col=c("lightblue", "cyan", "yellow"), 
                                        lwd=c(2,2,2)), 
                        MeanLine=list(var="PRODUCT", col="blue", lwd=2),
                        Title=list(main=paste("Variability Chart. Estimate: intercept (", vc.user()$emu,"), top level sd= (",vc.user()$etop,")"
                                              , ",\n middle level sd= (",vc.user()$eday,"), lowest level sd= (",vc.user()$erun,") & random error sd= (",vc.user()$esigma,")")),
                        
                        # MeanLine=list(var="mid", col="pink", lwd=2),
                        Points=list(pch=list(var="BATCH", pch=c(21, 22, 24)), 
                                    bg =list(var="BATCH", bg=c("lightblue", "cyan", "yellow")), 
                                    cex=1.25)) 
            }
            
        }  
    })   
    
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run the application 
shinyApp(ui = ui, server = server)
    
 