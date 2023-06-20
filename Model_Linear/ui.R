#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: https://doi.org/10.3390/molecules24122298
#This routine performs the model-fitting linear regression analysis.
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
    windowTitle = "THINKS/Model_Linear",
    title = tags$head(tags$link(rel="icon", 
                                href="favicon.ico", 
                                type="image/x-icon")
    )),
  h2("Model fitting (linear regression)", style="margin-bottom: 0em"), br(),
  tags$head(
    tags$style(HTML("
                    .MathJax_Display {
                    text-align: left !important;
                    }
                    .initialParent {
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
                    .multicol{
                              height:auto;
                              -webkit-column-count: 4;
                              -moz-column-count: 4;
                              column-count: 4;
                    }    
                    #text {
                              width: 100%;
                              height: 220px;
                              overflow: auto;
                              font-size: 10px;
                              margin-top: 2em; 
                              margin-bottom: 2em;
                              white-space: pre-wrap;
                    }
                    hr { 
                          margin-top: 0em !important;
                          margin-bottom: 1em !important;
                          border-width: 1px;
                        } "))
  ),

  sidebarLayout(
    sidebarPanel(width=3,
        tabsetPanel(
            tabPanel("Basic",
                 br(),       
                 #h4("Input data"),
                 helpText("Load the .txt prepared with ", tags$a(href = "http://www.thinks.chemphys.ru/shiny/Import/", target="_blank", "Import")," app"),
                 fileInput('file1', label=NULL, accept=c('text/plain')),
                 conditionalPanel(
                    condition = "input.go != 0",
                    helpText("Data file structure: row i - T[°C], row i+1 - t[s], row i+2 - α[0..1], row i+3 - dα/dt[1/s]. 
                        Separator: space")),
                 fluidRow(column(6, numericInput("aMin", label="Use conversion from:", 0.05, min = 0.01, max = 0.3, step=0.01)),
                          column(6, numericInput("aMax", label="Use conversion up to:", 0.95, min = 0.6, max = 0.99, step=0.01))
                          ),
                 
                 actionButton("go", "Calculate", class="btn btn-primary"), br(),br(),
                 checkboxInput("summaryLR", "Show full regression summary", value=F),
                 helpText("Select models to plot [2]:"),
                 tags$div(align = "left", 
                          class = "multicol",
                          checkboxGroupInput("reMod", label=NULL,
                                             c("F0" = "F0","F1" = "F1","F2" = "F2","F3" = "F3","R2" = "R2","R3" = "R3",
                                               "A2" = "A2","A3" = "A3","A4" = "A4",
                                               "P2" = "P2","P3" = "P3","P4" = "P4","P23" = "P23",
                                               "L2" = "L2","B1"="B1",
                                               "D1" = "D1","D2" = "D2","D3" = "D3","D4" = "D4"
                                               ), selected=c("F1","A2","R2","D2"))
                          )
            ),
            tabPanel("Advanced",
                br(),
                checkboxInput("DL", "Define constraints for kinetic parameters", value=F),
                conditionalPanel(
                  condition = "input.DL == true",
                  fluidRow(column(7, 
                      sliderInput("lA","Preexponent lncA, 1/s:",min = 1,max = 100,value = c(1,100), step=1),
                      sliderInput("Ea","Activation energy Ea, kJ/mol:",min = 10,max = 500,value = c(50,300), step=10)),
                      column(5, sliderInput("n1","Reaction order n:",min = 0,max = 3,value = c(0,2), step=0.1),
                      sliderInput("m1","Reaction order m:",min = -2,max = 3,value = c(-1,2), step=0.1)))
                ),
                hr(),
                helpText("Plot setting"),
                checkboxInput("Compact", "Compact mode (adjusted plot heights)", value=T),
                conditionalPanel(
                  condition = "input.Compact == false", 
                  fluidRow(
                    column(7,helpText("Basic height (px):")),
                    column(5,numericInput("BaseHeight", label=NULL, 350, min = 100, max = 1000, step=50, width='100%'))
                  )),
                radioButtons("TitleFont", label="Title font style",c("normal"=1,"bold"=2,"italic"=3, "bold italic"=4), 
                             selected =1, inline =T),
                helpText("Font size for labels (cex):"),
                fluidRow(column(4, helpText("Title"),
                                numericInput("TitleSize", label=NULL, 1.6, min = 0.8, max = 2.5, step=0.1, width='100%')),
                         column(4, helpText("Axis title"),
                                numericInput("AxisSize", label=NULL, 1.8, min = 0.8, max = 2.5, step=0.1, width='100%')),
                         column(4, 
                                helpText("Legend"),
                                numericInput("LegendSize", label=NULL, 1, min = 0.5, max = 2.5, step=0.1, width='100%'))
                )
            ))
    ),
    
    mainPanel(
      column(width=7,
       tabsetPanel(id = "MainTabset", selected = "panel2",
          tabPanel("Results", value = "panel1",   
                   conditionalPanel(
                     condition = "input.go != 0",
                     plotOutput("distPlot", height = "auto"),
                     plotOutput("distPlot2", height = "auto"))
          ),
          tabPanel("Help", value = "panel2",  
                   div(class = "initialParent", 
                       div( class="initialChild",
                            p(strong("/ Instructions /")), 
                  p("1.  Load the input data that were created with ", 
                              tags$a(href = "http://www.thinks.chemphys.ru/shiny/Import/", target="_blank", "Import"), 
                              "app or manually. It should be a single .txt file with rows of equal length, four rows for each experiment:
                    (i) temperature [°C], (i+1) time [s], (i+2) conversion [0..1], (i+3) conversion rate [1/s], 
                    , separator - space. An example can be found ",
                              tags$a(href = "example_tnaa_tga_2-20.txt", target="_blank", "here"),
                    
                    ". Note, that more than four experiments is needed to perform the kinetic analysis (see the ",
                    tags$a(href = "https://doi.org/10.1016/j.tca.2011.03.034", target="_blank", "ICTAC Kinetic committee recommendations"),
                    ")."),
                              p("2. Press the \"Calculate\" button and see the results of the linear regression for the flexible reaction model
                 with two exponents (also known as the reduced Sestak-Berggren, shown on the right) [1]. During the calculation
                 the conversion degree range can be specified. The allowed range for the kinetic parameters can be set using the 
                 \"Advanced\" tab, once the \"Define 
                 constraints\" is checked."),
                              p("3.  See the results. The full information on the regression is available once the relevant checkbox is selected. 
                    Various ideal reaction types (e.g., from [2]) can be selected to visually compare it on the master plot with the optimized 
                    flexible reaction model."
                     )
                       ))
          ))
      ),
      column(width=5,
             
             h5("Model fitting (linear regression) [1]"),
             #h6("Idea: linearization of the following equation [1]:"),
             withMathJax(),
             div("$$ \\large  
                 {
                 ln\\left(\\frac{d\\alpha}{dt}\\right) - ln[(1-\\alpha)^n\\alpha^m]= ln(c A) - \\frac{E_a}{RT}
                 }
                 $$"),
             
             h5("References:"),
             p(h6("1. L.A. Pérez-Maqueda, J.M. Criado, P.E. Sánchez-Jiménez, Combined Kinetic Analysis 
                  of Solid-State Reactions: A Powerful Tool for the Simultaneous Determination of Kinetic 
                  Parameters and the Kinetic Model without Previous Assumptions on the Reaction Mechanism, 
                  J. Phys. Chem. A. 110 (2006) 12456–12462. doi:10.1021/jp064792g.",br(),
                  "2. S.V. Vyazovkin, A.K. Burnham, J.M. Criado, L.A. Pérez-Maqueda, C. Popescu, N. Sbirrazzuoli,
                  ICTAC Kinetics Committee recommendations for performing kinetic computations on thermal analysis data, Thermochimica Acta 520 (2011) 1-19,
                  doi:10.1016/j.tca.2011.03.034.")),
             
             conditionalPanel(
               condition = "input.summaryLR == true",
               verbatimTextOutput("text")),
             conditionalPanel(
               condition = "input.summaryLR == false", br(),hr()),
             img(src="thinks.png", height = 75, width = 100),#, align='right'),
             h5(strong("Dr. Nikita V. Muravyev, Version 26.06.20")),
             p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
                  a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
                  a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
             )
             
             )
    )
  )
))
