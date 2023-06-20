#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine performs the purity analysis using the DSC data on melting of the respective substance.
#
#        Δ  ****  __|__  ****  Δ
#       / \      /.....\      / \
#       |Π|_____|_______|_____|Π|
#      /   \  ■ | T T T | ■  /   \
#     | ICP | ■ | | | | | ■ | RAS |
#     |_____|___|_|Π|Π|_|___|_____|
#
#Created by Dr. N.V. Muravyev (2016-2020).
#--------------------------------------------------------------------------------------------------------------------------

library(shiny)
library(shinythemes)

shinyUI(fluidPage(theme = shinytheme("yeti"),
  titlePanel(
    windowTitle = "THINKS",
    title = tags$head(tags$link(rel="icon", 
                                href="favicon.ico", 
                                type="image/x-icon")
    )),
  h2("Single curve and purity analysis"), br(),
  tags$head(
    tags$style(HTML("
                    .MathJax_Display {
                    text-align: left !important;
                    }")),
    tags$style(HTML(".initialParent {
                              background-color: AntiqueWhite ;
                              margin-top: 0.8em;
                              padding: 1em;
                              }
                    .initialChild {
                              text-align: justify;
                              text-indent: -1.5em;
                              margin-left: 2em;
                              margin-top: 0.5em;
                              margin-right: 1em;
                              margin-bottom: 1.2em;
                    }
                      #textPurity {
                              font-size: 12 px;
                              padding: 1em;
                              margin-top: 0 em;
                              margin-bottom: 0em;
                      }
                    hr { 
                          margin-top: 1em !important;
                          margin-bottom: 0.5em !important;
                          border-width: 1px;
                        }"))
  ),
  
  sidebarLayout(
    sidebarPanel(width=3,
      radioButtons("task", label=NULL,c("Simple analysis" = "kq", "Purity analysis"="pa"), inline=T),   
      fileInput('file1', 'Load the .txt file with DSC data', accept=c('text/plain')),
      fluidRow(column(3, helpText("Separator:"),
                      textInput("sep1",label=NULL,value = ";")
                      ),
                column(9,
                      fluidRow(
                            column(6, helpText("Drop X first lines: "),
                                    numericInput("BegSkipL",label=NULL,value = 35, min = 0, max = 100, step = 1)
                                    ),
                            column(6, helpText("Drop X last lines: "),
                                  numericInput("EndSkipL",label=NULL,value = 0, min = 0, step = 1)
                                    )
                              )
                        )
                ),
      fluidRow(column(4, helpText("Column #"),
                      selectInput("tempCol", label=NULL,#width="75px",
                                  choices = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),selected = 1)),
               column(4,
                      helpText("Parameter"),
                      strong("Temperature")),
               column(4,helpText("Units"),
                      selectInput("tempUnits", label=NULL, #width="75px",
                                  choices = c("°C"=1, "K"=2),selected = 1))
      ),
      fluidRow(column(4, selectInput("timeCol", label=NULL,#width="75px",
                                     choices = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),selected = 2)),
               column(4, strong("Time")),
               column(4, selectInput("timeUnits", label=NULL,#width="75px",
                                     choices = c("s"=1, "min"=60, "h"=3600),selected = 60))
      ),
      fluidRow(column(4, selectInput("dscCol", label=NULL,#width="75px",
                                     choices = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),selected = 3)),
               column(4, strong("DSC")),
               column(4, selectInput("dscUnits", label=NULL,#width="95px",
                                     choices = c("W/g" = 1, "mkV/mg" = 2, "mW" = 3, "mkV" = 4),selected = 1))
      ),
      fluidRow(column(6,numericInput("mass1","Sample mass [mg]",value = 3.075,min = 0.05, max = 500, step = 0.1)),
               column(6, numericInput("N_Simp",label="Simpson subintervals:",value = 1000,min = 2, max = 5000, step = 1))
              ),
      #verbatimTextOutput("headFile"),
      actionButton("go", "Start", class="btn btn-primary")
    ),
    mainPanel(
      fluidRow(
            tabsetPanel(id = "MainTabset", selected = "panel3",
              tabPanel("Plot", value="panel1",
                       column(7,
                          plotOutput("distPlot", height = 350),
                          plotOutput("distPlot2", height = 300),
                          #hr(),
                          verbatimTextOutput("text3")),
                          
                       column(4,
                              p(),
                              checkboxInput("ThinData", label = "Thin-out data", value = F),
                              conditionalPanel(
                                condition = "input.ThinData == true", 
                                uiOutput("thinOut")),
                              
                              radioButtons("BLtype", "Baseline type:",c("Linear" = "line","Linear-2" = "linst",
                                                                        "Spline" = "spln", "Spline-2" = "spln2"),inline=T),
                              checkboxInput("BLsmooth", label = "Smoothing data for baseline", value = F),
                              #checkboxInput("TimeConst", label = "Time constant calculation", value = F),
                              uiOutput("BLpointsBeg"),
                              uiOutput("BLpointsEnd"),
                              radioButtons("TRtype", "Thermal resistance:", c("Specify"="us", 
                                                                              "Average slope in 5-95% area range" = "ave5",
                                                                             "Average slope in 10-90%" = "ave10",
                                                                             "Minimal gradient" = "min", 
                                                                             "5% around min.grad." = "near5",
                                                                             "15% around min.grad." = "near15"
                                                                             ), selected="ave5", inline=T),
                              conditionalPanel(
                                condition = "input.TRtype == \"us\"", 
                                    div(style="display:inline-block", helpText("Thermal resistance [K/mW]:")),
                                    div(style="display:inline-block; margin-left: 0.2em", 
                                        numericInput("usTR",label=NULL, value = 0.068,min = 0.001, max = 0.999, width='100%')
                                        )
                                ),
                              #actionButton("goButton", "Add to Table"),
                              hr(),
                              img(src="thinks.png", height = 75, width = 100),#, align='right'),
                              h5(strong("Dr. Nikita V. Muravyev, Version 01.07.20")),
                              p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
                                   a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
                                   a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
                              ))
                       ),
              tabPanel("Purity", value="panel2",
                       column(7,
                              plotOutput("PurityPlot", height = 300),
                              verbatimTextOutput("textPurity"),
                              plotOutput("purityPlot2", height = 240),
                              conditionalPanel(
                                condition = "input.ShowRegr == true", 
                                verbatimTextOutput("textPurity2"))
                              ),
                       column(4,
                              fluidRow(column(6,
                                              radioButtons("PurCorr", "Apply Correction:",c("yes" = "y", "no" = "n"),inline=T),
                                              numericInput("MolarMass","Molar mass [g/mol]:",value = 122,min = 1, max = 1500, step = 0.1)
                                              ),
                                       column(6,
                                              radioButtons("PurrThin","Thin-out data:",c("yes" = "y", "no" = "n"), selected="y", inline=T),
                                              conditionalPanel(
                                                condition = "input.PurrThin == \"y\"", 
                                                numericInput("PT_pts","N points:",value = 20,min = 5, max = 500, step = 1)
                                              )
                                              )
                                       ),
                              sliderInput("Flim", "Final Fraction :",min = 1, max = 98, value = c(10,50), step=1),
                              checkboxInput("ShowRegr", label = "Show regression summary", value = F),
                              
                              hr(),
                              h4("Purity analysis"),
                              withMathJax(),
                              div("$$ \\Large
                                   {
                                    
                                    T_{fus} = T_{0} - \\frac{RT_{0}^2 }{\\Delta H_{fus}} x_{2} \\frac{1}{F}
                                    
                                   }
                                   $$"),
                              #ln\\left(\\frac{d\\alpha}{dt}\\right)_{\\alpha,i} = ln[f(\\alpha) A_\\alpha] - \\frac{E_\\alpha}{RT_{\\alpha,i}}
                              
                              h4("References:"),
                              p(h6("1. R. Blaine and C. Schoff, eds., Purity Determinations by Thermal Methods. 
                                   (ASTM International, 1984), https://doi.org/10.1520/STP838-EB"
                                   )
                              ),
                              hr(),
                              img(src="thinks.png", height = 75, width = 100),#, align='right'),
                              h5(strong("Dr. Nikita V. Muravyev, Version 01.07.20")),
                              p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
                                   a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
                                   a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
                              )
                              )
                       ),
              tabPanel("Help", value="panel3",
                       column(7,
                              div(class = "initialParent", 
                                  div( class="initialChild",
                                       p(strong("/ Instructions /")), 
                                       p("1. Load the DSC data of melting of the target substance. The input file should contain the temperature, time, and heat flow data (
                                        here is an ",
                                       tags$a(href = "example_Benz.acid_20bar.txt", target="_blank", "example"),")."),
                                       p("2. Define the number of columns with signals, other file details, and press \"Start\" button. 
                                       Select the interested region (skip lines at the beginning and the end if needed),
                                          tune the baseline and other details."),
                                       p("3. If the experiment is for the standard sample (e.g., indium), get the thermal resistance to use for the target substance.
                                         If it is already the target sample, specify the thermal resistance, select the \"Purity analysis\" task
                                         and go to the \"Purity\" tab."),
                                        p("4. Define the molar mass of the substance and the area range that is used for calculation.
                                          Note, that the method has a certain limitations concerning the nature of impurities (see relevant literature)
                                          and has to be applied for >95 (preferably >98.5) mol.% pure samples. Moreover, check the correction parameters 
                                          in the results, it should be less than 15%, otherwise the Vant-Hoff equation fails to describe the process and 
                                          the results are not reliable."),
                                       p("5.  Do thermal analysis with pure samples regularly)")
                                      
                              ))),
                       column(4,
                              br(),
                              img(src="thinks.png", height = 75, width = 100),#, align='right'),
                              h5(strong("Dr. Nikita V. Muravyev, Version 01.07.20")),
                              p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
                                   a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
                                   a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
                              ))
                         )
              )
          )
  )
)
))
