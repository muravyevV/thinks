#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine is used to perform the model fitting (specifically, the nonlinear regression)
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
library(shinyjs)
library(shinythemes)

shinyUI(fluidPage(theme = shinytheme("yeti"),useShinyjs(),
  titlePanel(
    windowTitle = "THINKS/Model_Nonlinear",
    title = tags$head(tags$link(rel="icon", 
                                href="favicon.ico", 
                                type="image/x-icon")
    )),
  h2("Model fitting (nonlinear regression)"),
  br(),
  
  tags$head(
    tags$style(HTML("
                    .MathJax_Display {
                    text-align: left !important;
                    }")),
    tags$style(HTML("        #div.checkbox {margin-top: 0px;}
                               .initialParent {
                              background-color: AntiqueWhite ;
                              margin-top: 0.8em;
                              padding: 0.8em;
                              }
                    .initialChild {
                              text-align: justify;
                              font-size: 12px;
                              text-indent: -1.5em;
                              margin-left: 2em;
                              margin-top: 0em;
                              margin-right: 1em;
                              margin-bottom: 0em;
                    }
                    #reggrSummary {
                              width: 100%;
                              height: 250px;
                              overflow: auto;
                              font-size: 12px;
                              margin-top: 1.5em; 
                              margin-bottom: 1em;
                              white-space: pre-wrap;
                    }
                    .radioDiff {
                              margin-left: 0em;
                              padding: 0em;
                              font-size: 9px;
                              }
                    hr { 
                          margin-top: 1em !important;
                          margin-bottom: 1em !important;
                          border-width: 1px;
                        }"))
  ),
  
  sidebarLayout(
    sidebarPanel(width=3,
        tabsetPanel(
           tabPanel("Initial",
                 #h4("Input data"),
                 br(),
                 helpText("Load the input data prepared with ", tags$a(href = "http://www.thinks.chemphys.ru/shiny/Import/", target="_blank", "Import")," app"),
                 fileInput('file1', label=NULL, accept=c('text/plain')),
                 conditionalPanel(
                   condition = "input.datasetsCnt == 2",
                   fileInput('file2', label=NULL, accept=c('text/plain'))
                 ),
                 helpText("Specify the initial guesses"),
                 fluidRow(column(7, conditionalPanel(
                                      condition = "input.lnA1f == false",
                                      sliderInput("lnA1","Preexponent lnA1, 1/s:",min = 5,max = 120,value = 43, step=0.1)
                                 )),
                          column(5, conditionalPanel(
                            condition = "input.taskType != 3&&input.taskType != 2&&input.nf1 == false&&input.taskType != 8",
                            sliderInput("n1","Reaction order n1:",min = 0,max = 3,value = 2, step=0.01)
                          )
                                 )
                 ),
                 fluidRow(column(7, conditionalPanel(
                                      condition = "input.Ea1f == false",
                                      sliderInput("Ea1","Activation energy Ea1, kJ/mol:",min = 50,max = 400,value = 181, step=0.1))),
                          column(5, conditionalPanel(
                                      condition = "input.taskType != 3&&input.mf1 == false&&input.taskType != 2&&input.taskType != 8",
                                      sliderInput("m1","Exponent m1:",min = -1,max = 2,value = 0, step=0.01)
                          )
                                 )
                 ), 
                 fluidRow(column(7, conditionalPanel(
                                      condition = "input.taskType != 1&&input.taskType != 2&&input.lnA2f == false",
                                      sliderInput("lnA2","Preexponent lnA2, 1/s:",min = 5,max = 120,value = 17.4, step=0.1))),
                          column(5, conditionalPanel(
                                      condition = "input.taskType == 3&&input.datasetsCnt == 2&&input.lnA22f == true",
                                      sliderInput("lnA22","lnA22, 1/s:",min = 5,max = 120,value = 17.4, step=0.1)
                                      ), 
                                 conditionalPanel(
                                      condition = "input.taskType >= 4&&input.taskType != 8&&input.nf2 == false",
                                      sliderInput("n2","Reaction order n2:",min = 0,max = 3,value = 0.1, step=0.01)
                          ) )
                 ),
                 fluidRow(column(7, conditionalPanel(
                                      condition = "input.taskType > 2&&input.Ea2f == false",
                                      sliderInput("Ea2","Activation energy Ea2, kJ/mol:",min = 50,max = 400,value = 96, step=0.1))),
                          column(5, conditionalPanel(
                                      condition = "(input.taskType == 3||input.taskType == 8)&&input.muf == false",
                                      sliderInput("mu","Volume factor mu:",min = 0,max = 1,value = 0.1, step=0.01)
                                    ), 
                                    conditionalPanel(
                                      condition = "input.taskType >= 4&&input.taskType != 8&&input.mf2 == false",
                                      sliderInput("m2","Exponent m2:",min = -1,max = 3,value = 0.22, step=0.01)
                          ))
                 ),
                 fluidRow(column(7, conditionalPanel(
                   condition = "input.taskType > 7&&input.lnA3f == false",
                   sliderInput("lnA3","Preexponent lnA3, 1/s:",min = 5,max = 120,value = 17.4, step=0.1))),
                   column(5, conditionalPanel(
                     condition = "input.taskType >= 8&&input.nf3 == false",
                     sliderInput("n3","Reaction order n3:",min = 0,max = 3,value = 0.1, step=0.01)
                   ) )
                 ),
                 fluidRow(column(7, conditionalPanel(
                   condition = "input.taskType > 7&&input.Ea3f == false",
                   sliderInput("Ea3","Activation energy Ea3, kJ/mol:",min = 50,max = 400,value = 96, step=0.1))),
                   column(5, conditionalPanel(
                     condition = "input.taskType >= 8&&input.mf3 == false",
                     sliderInput("m3","Exponent m3:",min = -1,max = 3,value = 0.22, step=0.01)
                   ))
                 ),
                 
                 conditionalPanel(
                   condition = "input.taskType > 4&&input.cd1f == false",
                   sliderInput("cd1","Fraction of first reaction cd1:",min = -3,max = 4,value = 0.5, step=0.01)
                 ),
                 
                 conditionalPanel(
                   condition = "input.taskType > 4&&input.datasetsCnt == 2&&input.taskType !=7",
                   sliderInput("cd2","Fraction of first reaction for 2nd dataset cd2:",min = 0,max = 1,value = 0.5, step=0.01)
                 ),
                
                 fluidRow(column(5, actionButton("go", "Run analysis", class="btn btn-primary")), #icon("brain"),
                          column(6, uiOutput("setButton")
                            )#,
                        #  actionButton("stop", "Stop", class="btn btn-danger")
                 )
           ),
           tabPanel("Plot settings",
                    br(),
                    checkboxInput("normC", "Normalize curves (x weights)", value=F),
                    checkboxInput("iG", "Show initial guess", value=T),
                    radioButtons("drawX", "Draw as X:", c("T[K]"=2,"T[°C]"=1,"time[min]"=3,"time[s]"=4), 
                                 selected =2, inline =T),
                    conditionalPanel(
                      condition = "input.datasetsCnt == 2",
                      radioButtons("drawX2", "(Second plot) Draw as X:", c("T[K]"=2,"T[°C]"=1,"time[min]"=3,"time[s]"=4), 
                                   selected =2, inline =T)
                    ),       
                    fluidRow(column(8, helpText("Show every X point (rarefy):")),
                             column(4, numericInput("PointThin", label=NULL, 1, min = 1, max = 20, step=1, width='100%')
                             )),
                    checkboxInput("Compact", "Compact mode (adjusted plot heights)", value=T),
                    conditionalPanel(
                      condition = "input.Compact == false", 
                      fluidRow(
                        column(7,helpText("Basic height (px):")),
                        column(5,numericInput("BaseHeight", label=NULL, 350, min = 100, max = 1000, step=50, width='100%'))
                      )),
                   # radioButtons("TitleFont", label="Title font style",c("normal"=1,"bold"=2,"italic"=3, "bold italic"=4), 
                  #               selected =1, inline =T),
                    
                  fluidRow(column(4, numericInput("AxisSize", label="Axis title size", 1.8, min = 0.8, max = 2.5, step=0.1, width='100%'), 
                                  numericInput("PointSize", label="Point size", 1, min = 0.1, max = 6, step=0.1, width='100%')),
                          column(4, numericInput("LegendSize", label="Labels size", 1.1, min = 0.5, max = 2.5, step=0.1, width='100%'),
                                 numericInput("PointType", label="Point type", 1, min = 1, max = 18, step=1, width='100%')),
                          column(4, numericInput("nlrLwd", label="Line width", 3, min = 0.1, max = 6, step=0.1, width='100%'),
                                 numericInput("marginTitle", label="Margin(axis)", 0, min = 0, max = 6, step=0.1, width='100%'))
                  ),
                  fluidRow(column(4, numericInput("marginLeft", label="Margin(left)", 0, min = 0, max = 6, step=0.1, width='100%')),
                           column(4, numericInput("marginBottom", label="Margin(bottom)", 0, min = 0, max = 6, step=0.1, width='100%')),
                           column(4, numericInput("maxAl", label="Maximal alpha", 1, min = 0, max = 3, step=0.05, width='100%'))
                           ),
                  
                  fluidRow(column(5, downloadButton('download1', 'Download Plot')),
                           column(6, conditionalPanel(
                             condition = "input.datasetsCnt == 2",
                             downloadButton('download2', 'Download Plot2')
                           )))
                  
                  
                  
                  
           )
        )
    ),
    
    mainPanel(
      column(width=7,
       tabsetPanel(id = "MainTabset", selected = "panel2",
          tabPanel("Results", value = "panel1",  
             plotOutput("Plot", height="auto"),
             conditionalPanel(
             condition = "input.datasetsCnt == 2",
               plotOutput("Plot2", height="auto")
             ),
             conditionalPanel(
               condition = "input.go != 0",
             verbatimTextOutput("reggrSummary"))
          ),
          tabPanel("Help", value = "panel2",  
              div(class = "initialParent", 
                 div( class="initialChild",
                    p(strong("/ Instructions /")), 
                    p("1. Load the input file created with ", 
                       tags$a(href = "http://www.thinks.chemphys.ru/shiny/Import/", target="_blank", "Import"), 
                       "app or manually. It should be a single .txt file with rows of equal length, four rows for each experiment:
                        (i) temperature [°C], (i+1) time [s], (i+2) conversion [0..1], (i+3) conversion rate [1/s], 
                        , separator - space. An example can be found ",
                       tags$a(href = "example_tnaa_tga_2-20.txt", target="_blank", "here."),
                      ". Note, that more than four experiments is needed to perform the kinetic analysis (see the ",
                      tags$a(href = "https://doi.org/10.1016/j.tca.2011.03.034", target="_blank", "ICTAC Kinetic committee recommendations"),
                      ")."),
                    p("2. Select the kinetic model, the type of signal which will be optimized, and define the initial guess of the kinetic parameters.
                      Six kinetic models are available so far: the single-step in flexible form (ePT), the single-step of a theoretical form,
                      autocatalytic model (DMM), and two reactions in flexible form (ePT) that are parallel, consecutive or independent. Three types
                      of signals that can be optimized are the integral, the differential, and the weighted differential (data is scaled by the 
                      amplitude to minimize the difference in contribution from the fast and slow heating rates). The general recommendation is to optimize
                      the same type of data that was originally measured.
                      Define the initial guesses for the kinetic parameters, the more close they are, the lesser the calculation time will be and the 
                      higher is the probability of converging
                      (use some insights from the other types of the kinetic analysis, e.g.,",
                      tags$a(href = "http://www.thinks.chemphys.ru/shiny/Isoconversional/", target="_blank", "Isoconversional"), 
                      "). The ODE solver type is selected further, start with a low-order solvers due to its robustness."),
                    p("3. Press the \"Run analysis\" button and see the results of the nonlinear regression for the 
                        selected formal model. The general guideline is to move from the simple kinetic models to the complex.
                      Same for ODE solvers - use solvers above \"Euler\" for the finishing kinetic analysis, but apply the high-order 
                       only if they benefit to the accuracy of the kinetic parameters, converging of the optimization, or the stability of the
                       results (e.g., for stiff v5 models at a certain set of kinetic parameters the results can be unstable with the low-order solvers).
                       Please, keep in mind, that increasing the computational complexity reduces the server resources available the other users.
                      When the close result is available, use \"Set the current results\" button to make the current results of nonlinear regression to be the initial guesses
                      for the next calculation. Use \"Fix parameters\" tab to fix some parameters throghout the analysis. To compare the various models use the Bayes IC in the \"Advanced\" tab,
                      when the models built for the same dataset are compared, the lower BIC value is preferred."),
                    p("4. Tune the plot settings (in the respective tab) to get the pleasant figure, right click to save or copy. 
                      Have a good day and do kinetics regularly!")
              ))
          )
       )
      ),
      column(width=5,
        tabsetPanel(
          tabPanel("Settings",
             br(),
             radioButtons("dataTp", "Optimize signal:", c("integral"=2,"differential"=1,"weighted differential"=3), 
                          selected =2, inline =T),
             #radioButtons("methodDiff", "Solver for ODEs:", c("Euler"=1,"Improved Euler"=2,"Runge-Kutta 4th"=4,"Runge-Kutta 5th"=5), 
            #              selected =1, inline =T),
             helpText("ODE solvers [order in brackets]:"),
            div(class="radioDiff",radioButtons("methodDiff", label=NULL, c("Euler [1]"=1,"Runge-Kutta [2]"=2,"Runge-Kutta [4]"=4,
                                                      #"5"=5,"6"=6,"Runge-Kutta Cash-Karp [4(5)]"=7, 
                                                      "Dormand-Prince [4(5) local 7]"=8, 
                                                      "Petzold-Hindmarsh, lsoda"=9), 
                          selected =1, inline =T)),
             selectInput("taskType", "Kinetic model:", c("v1: Flexible single step (ePT)"=1,
                                                  "v2: Simple reaction types"=2,
                                                  "v3: Autocatalysis (DMM)"=3,
                                                  "v4: Parallel ePT || ePT"=4,
                                                  "v5: Consecutive ePT -> ePT"=5,
                                                  "v6: Independent ePT + ePT"=6,
                                                  "v7: ePT sw(a) ePT"=7,
                                                  "v8: Consecutive DMM -> ePT"=8,
                                                  "v9: ePT -> ePT || ePT"=9,
                                                  "v10: ePT sw(a) ePT -> ePT"=10
                                                  ), selected = 1),
             conditionalPanel(
               condition = "input.taskType == 2",
               selectInput("reMod", "Simple reaction model:", c("Zero-order reaction (F0)" = "F0","First-order reaction (F1)" = "F1",
                                     "Second-order reaction (F2)" = "F2","Third-order reaction (F3)" = "F3",
                                     "KJMAE nucleation-growth (A2)" = "A2","KJMAE nucleation-growth (A3)" = "A3","KJMAE nucleation-growth (A4)" = "A4",
                                     "Contracting cylinder (R2)" = "R2","Contracting sphere (R3)" = "R3",
                                     "Power law (P2)" = "P2","Power law (P3)" = "P3","Power law (P4)" = "P4","Power law (P23)" = "P23",
                                     "One-dimensional diffusion (D1)" = "D1","Two-dimensional diffusion (D2)" = "D2",
                                     "3D Jander diffusion (D3)" = "D3","3D Ginstling-Brounshtein diffusion (D4)" = "D4",
                                     "Polymer random scission (PRS2)" = "L2","Classical Prout-Tompkins (B1)"="B1"), selected = "F1")
             ),
            conditionalPanel(
              condition = "(input.taskType !=2&&input.taskType !=3)&&!(input.taskType==1&&input.mf1 == true&&input.m1fixed==10)",
              h4("Flexible single step (ePT):"),
              withMathJax(),
              div("$$ \\Large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k\\left( 1 - \\alpha \\right)^n \\left[ 1 - q \\left( 1 - \\alpha \\right) \\right]^m
                 }
                 $$"),
              p(h6("1. Burnham, A. K.; Zhou, X.; Broadbelt, L. J. Critical Review of the Global Chemical Kinetics of Cellulose Thermal Decomposition. Energy Fuels 2015, 29 (5), 2906–2918. doi:10.1021/acs.energyfuels.5b00350."))
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"F0\"",
              h4("Zero-order reaction:"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k
                 }
                 $$") ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"F1\"",
              h4("First-order reaction:"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)
                 }
                 $$") ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"F2\"",
              h4("Second-order reaction:"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)^2
                 }
                 $$") ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"F3\"",
              h4("Third-order reaction:"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)^3
                 }
                 $$") ),
            conditionalPanel(
              condition = "(input.taskType == 2 && (input.reMod == \"A2\"||input.reMod == \"A3\"||input.reMod == \"A4\")) ||(input.taskType != 7&&input.taskType != 10&&input.mf1 == true&&input.m1fixed==10)",
              h4("KJMAE nucleation-growth (An, n = 2, 3, 4):"),
              withMathJax(),
              div("$$ \\large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k n( 1- \\alpha) [-ln(1-\\alpha)]^{(n-1)/n}
                 }
                 $$"),
              p(h6("1.	Kolmogorov A. A statistical theory for the recrystallization of metals. Izv Akad Nauk SSSR Ser Mat 1937:355–9.", br(),
                   "2. Johnson WA, Mehl RF. Reaction kinetics in processes of nucleation and growth. Trans AIME 1939;135:416–42.",br(),
                   "3. Avrami M. Kinetics of Phase Change. I General Theory. J Chem Phys 1939;7:1103. doi:10.1063/1.1750380.",br(),
                   "4. Erofeev BV. Generalized Equation of Chemical Kinetics and its Application to Reactions involving solid phase components. Dokl Akad Nauk USSR 1946;52:515–8."
              ))
            ),
            
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"R2\"",
              h4("Contracting cylinder (R2):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 2k ( 1- \\alpha)^{1/2}
                 }
                 $$")
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"R3\"",
              h4("Contracting sphere (R3):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 3k ( 1- \\alpha)^{2/3}
                 }
                 $$")
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && (input.reMod == \"P2\"|| input.reMod == \"P3\"||input.reMod == \"P4\"|| input.reMod == \"P23\")",
              h4("Power law (Pn, n = 2/3, 2, 3, 4):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k n \\alpha^{(n-1)/n}
                 }
                 $$")
            ), #------------------------------------------add reference here
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"D1\"",
              h4("One-dimensional diffusion (D1):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{1}{2\\alpha}
                 }
                 $$")
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"D2\"",
              h4("Two-dimensional diffusion (D2):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{1}{-ln(1-\\alpha)}
                 }
                 $$")
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"D3\"",
              h4("3D Jander diffusion (D3):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{3(1-\\alpha)^{2/3}}{2[1-(1-\\alpha)^{1/3}]}
                 }
                 $$"),
              p(h6("1.Jander W. Reaktionen im festen Zustande bei höheren Temperaturen. Reaktionsgeschwindigkeiten endotherm verlaufender Umsetzungen. Z Für Anorg Allg Chem 1927;163:1–30. doi:10.1002/zaac.19271630102."))
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"D4\"",
              h4("3D Ginstling-Brounshtein diffusion (D4):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{3}{2[(1-\\alpha)^{-1/3}-1]}
                 }
                 $$"),
              p(h6("1.Ginstling AM, Brounshtein BI. The diffusion kinetics of reactions in spherical particles. J Appl Chem USSR 1950;23:1249–59."))
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"L2\"",
              h4("Random scission of polymer chain (PRS2):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 2k (\\alpha^{1/2}-\\alpha)
                 }
                 $$"),
              p(h6("1. Flynn JH, Wall LA. General treatment of the thermogravimetry of polymers. J Res Natl Bur Stand Sect Phys Chem 1966;70A:487. doi:10.6028/jres.070A.043.",br(),
                   "2. Simha R, Wall LA. Kinetics of Chain Depolymerization. J Phys Chem 1952;56:707–15. doi:10.1021/j150498a012."))
            ),
            conditionalPanel(
              condition = "input.taskType == 2 && input.reMod == \"B1\"",
              h4("Classical Prout-Tompkins (B1):"),
              withMathJax(),
              div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\alpha (1-\\alpha)
                 }
                 $$"),
              p(h6("1.Prout EG, Tompkins FC. The thermal decomposition of potassium permanganate. Trans Faraday Soc 1944;40:488. doi:10.1039/tf9444000488."))
            ),
            conditionalPanel(
              condition = "input.taskType == 3 || input.taskType == 8",
              h4("Autocatalysis (DMM):"),
              withMathJax(),
              div("$$ \\large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k_1 \\left( 1 - \\alpha \\right) + k_2 \\left( 1 - \\mu \\right) \\frac{\\alpha \\left(1 - \\alpha \\right)}{1 - \\mu \\alpha }
                 }
                 $$"),
              p(h6("1. Dubovitskii, F. I.; Manelis, G. B.; Merzhanov, A. G. Formal Kinetic Model of Thermal Decomposition of Explosives in Liquid State. Trans Acad. Nauk SSSR Dokl. 1958, 121 (4), 668–670."))
            ),
             hr()
             
          ),
          tabPanel("Advanced",
              br(),
              div(style="display:inline-block", helpText("Amount of datasets:")),
              div(style="display:inline-block; margin-left: 0.2em", 
                  radioButtons("datasetsCnt", label=NULL, c("one"=1,"two"=2), selected =1, inline =T)),
              helpText("ODE solvers used for plot of initial guess:"),
              div(class="radioDiff",radioButtons("methodIni", label=NULL, c("Euler [1]"=1,"Runge-Kutta [2]"=2, 
                                                                            #"Runge-Kutta [4]"=4,
                                                                            #"5"=5,"6"=6,"Runge-Kutta Cash-Karp [4(5)]"=7, 
                                                                            "Dormand-Prince [4(5) local 7]"=8, 
                                                                             "Petzold-Hindmarsh, lsoda"=9),
                                                 selected =1, inline =T)),
              checkboxInput("Trunc_al", "Use all data (0..1 alpha's)", value=T),
              conditionalPanel(
                condition = "input.Trunc_al == false",
                fluidRow(column(5,numericInput("Trunc_al_min", label="Use data from alpha:", 0,min = 0, max = 0.4, width= '100%')),
                         column(5,numericInput("Trunc_al_max", label="Use data till alpha:", 1,min = 0.6, max = 1, width= '100%'))
                         )
                ),
              fluidRow(column(5,numericInput("al_min", "Starting alpha for NLR:", 0.00005,min = 0.0000005, max = 0.4)),
                       column(5,numericInput("al_max", "Limiting alpha for NLR:", 0.99995,min = 0.6, max = 0.9999995))
                       ),
              checkboxInput("suppression", "Suppress NaN's and outliers after solver", value=T),
              radioButtons("optimType", "Nonlinear regression algorithm:", c("Levenberg-Marquardt[nlsLM]"=1,"Bates-Chambers[port]"=2#,
                                                                   #"nlfb(Y)"=3, "nlxb(Y)"=4
                                                                   ), selected =1, inline =T),
              checkboxGroupInput("statM", "Statistical parameters of the optimized model:", c("RSS"=1,
                                                           #"Fisher"=2,
                                                           "Akaike IC"=3, "Bayes IC"=4), 
                                 selected = c(1,3,4), inline = T),
              
              #hr(),
              verbatimTextOutput("reggrSummary2"),
              helpText(textOutput("reggrSummaryTime"))
          ),
          tabPanel("Fixed parameters",
                   br(),
                   fluidRow(column(6,conditionalPanel(
                     condition = "input.taskType >= 9",
                     checkboxInput("lnA1f", "Fix preexponent lnA1")
                   )),
                   column(6, conditionalPanel(
                     condition = "input.lnA1f == true&&input.taskType >= 9",
                     div(style="display:inline-block", strong("lnA1:")),
                     div(style="display:inline-block; margin-left: 0.2em",numericInput("lnA1fixed", label=NULL, 15,min = 5, max = 120, width= '100%'))
                   ))),
                   fluidRow(column(6,checkboxInput("Ea1f", "Fix activation energy Ea1")),
                            column(6, conditionalPanel(
                              condition = "input.Ea1f == true",
                              div(style="display:inline-block", strong("Ea1:")),
                              div(style="display:inline-block; margin-left: 0.65em", numericInput("Ea1fixed", label=NULL, 150,min = 50, max = 350, width= '100%'))
                            ))),
                   fluidRow(column(6,conditionalPanel(
                                       condition = "input.taskType != 3&&input.taskType != 2&&input.taskType != 8",
                                       checkboxInput("nf1", "Fix reaction order n1")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.nf1 == true&&input.taskType != 3&&input.taskType != 2&&input.taskType != 8",
                                        div(style="display:inline-block", strong("n1:")),
                                        div(style="display:inline-block; margin-left: 1.15em", numericInput("n1fixed", label=NULL, 2,min = 0, max = 3, width= '100%'))
                            ))),
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType != 3&&input.taskType != 2&&input.taskType != 8",
                                        checkboxInput("mf1", "Fix exponent m1 (set as 10 to change ePT to KJMAE(n))")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.mf1 == true&&input.taskType != 3&&input.taskType != 2&&input.taskType != 8",
                                        div(style="display:inline-block", strong("m1:")),
                                        div(style="display:inline-block; margin-left: 0.9em", numericInput("m1fixed", label=NULL, 0.05,min = -1, max = 3, width= '100%'))
                            ))),
                   fluidRow(column(6,conditionalPanel(
                     condition = "input.taskType >= 9",
                     checkboxInput("lnA2f", "Fix preexponent lnA2")
                   )),
                   column(6, conditionalPanel(
                     condition = "input.lnA2f == true&&input.taskType >= 9",
                     div(style="display:inline-block", strong("lnA2:")),
                     div(style="display:inline-block; margin-left: 0.2em",numericInput("lnA2fixed", label=NULL, 15,min = 5, max = 120, width= '100%'))
                   ))),
                   fluidRow(column(6,conditionalPanel(
                     condition = "input.taskType != 1&&input.taskType != 2",
                     checkboxInput("Ea2f", "Fix activation energy Ea2 (if set as 0, Ea2=Ea1)")
                   )),
                   column(6, conditionalPanel(
                     condition = "input.Ea2f == true&&input.taskType != 1&&input.taskType != 2",
                     div(style="display:inline-block", strong("Ea2:")),
                     div(style="display:inline-block; margin-left: 0.65em",numericInput("Ea2fixed", label=NULL, 150,min = 50, max = 350, width= '100%'))
                   ))),
                   
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType >= 4&&input.taskType != 8",
                                        checkboxInput("nf2", "Fix reaction order n2")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.nf2 == true&&input.taskType >= 4&&input.taskType != 8",
                                        div(style="display:inline-block", strong("n2:")),
                                        div(style="display:inline-block; margin-left: 1.15em", numericInput("n2fixed", label=NULL, 1,min = 0, max = 3, width= '100%'))
                            ))),
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType >= 4&&input.taskType != 8",
                                        checkboxInput("mf2", "Fix exponent m2")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.mf2 == true&&input.taskType >= 4&&input.taskType != 8",
                                        div(style="display:inline-block", strong("m2:")),
                                        div(style="display:inline-block; margin-left: 0.9em", numericInput("m2fixed", label=NULL, 0,min = -1, max = 3, width= '100%'))
                            ))),
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType == 3||input.taskType == 8",
                                        checkboxInput("muf", "Fix volume factor mu")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.muf == true&&(input.taskType == 3||input.taskType == 8)",
                                        div(style="display:inline-block", strong("mu:")),
                                        div(style="display:inline-block; margin-left: 0.4em",numericInput("mufixed", label=NULL, 0,min = 0, max = 1, width= '100%'))
                            ))),
                   fluidRow(column(6,conditionalPanel(
                     condition = "input.taskType >= 8",
                     checkboxInput("lnA3f", "Fix preexponent lnA3")
                   )),
                   column(6, conditionalPanel(
                     condition = "input.lnA3f == true&&input.taskType >= 8",
                     div(style="display:inline-block", strong("lnA3:")),
                     div(style="display:inline-block; margin-left: 0.2em",numericInput("lnA3fixed", label=NULL, 15,min = 5, max = 120, width= '100%'))
                   ))),
                   fluidRow(column(6,conditionalPanel(
                     condition = "input.taskType >= 8",
                     checkboxInput("Ea3f", "Fix activation energy Ea3 (if set as 0, Ea3=Ea1)")
                   )),
                   column(6, conditionalPanel(
                     condition = "input.Ea3f == true&&input.taskType >= 8",
                     div(style="display:inline-block", strong("Ea3:")),
                     div(style="display:inline-block; margin-left: 0.65em",numericInput("Ea3fixed", label=NULL, 150,min = 50, max = 350, width= '100%'))
                   ))),
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType >= 8",
                                        checkboxInput("nf3", "Fix reaction order n3")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.nf3 == true&&input.taskType >= 8",
                                        div(style="display:inline-block", strong("n3:")),
                                        div(style="display:inline-block; margin-left: 1.15em", numericInput("n3fixed", label=NULL, 1,min = 0, max = 3, width= '100%'))
                            ))),                   
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType >= 8",
                                        checkboxInput("mf3", "Fix exponent m3")
                                  )),
                            column(6, conditionalPanel(
                                        condition = "input.mf3 == true&&input.taskType >= 8",
                                        div(style="display:inline-block", strong("m3:")),
                                        div(style="display:inline-block; margin-left: 0.9em", numericInput("m3fixed", label=NULL, 0,min = -1, max = 3, width= '100%'))
                            ))),

                   fluidRow(column(6, checkboxInput("TrigTf", "Set triggering temperature")),
                            column(6, conditionalPanel(
                                          condition = "input.TrigTf == true",
                                          div(style="display:inline-block", strong("T*:")),
                                          div(style="display:inline-block; margin-left: 1.25em",numericInput("TrigTfixed", label=NULL, 400,min = 50, max = 1050, width= '100%'))
                   ))),
                   conditionalPanel(
                     condition = "input.taskType == 3&&input.datasetsCnt == 2",
                     checkboxInput("lnA22f", "Allow different lnA2 for datasets")
                   ),
                   tags$table(style = "width: 100%; white-space:nowrap; valign: bottom",
                              tags$tr(tags$td(valign="middle", 
                                              conditionalPanel(
                                                condition = "input.taskType != 3&&input.taskType != 2",
                                              helpText("Initiation parameter "))),
                                      tags$td(valign="bottom",
                                        conditionalPanel(
                                        condition = "input.taskType >= 4&&input.taskType != 8&&input.taskType != 9&&input.taskType != 10",
                                          div(style="display:inline-block; margin-left:0.3em", helpText(" (")),
                                          div(style="display:inline-block; margin-top:0.2em", checkboxInput("q2f", "Specify q2 separately)", value=F))
                                        ),
                                        conditionalPanel(
                                          condition = "input.taskType == 8||input.taskType == 10",
                                          div(style="display:inline-block; margin-left:0.3em", helpText(" (")),
                                          div(style="display:inline-block; margin-top:0.2em", checkboxInput("q3f", "Specify q3 separately)", value=F))
                                        )
                                      ))),
                   
                   fluidRow(column(6,conditionalPanel(
                     condition = "input.taskType != 3&&input.taskType != 2&&input.taskType != 8",
                     div(style="display:inline-block", strong("q:")),
                     div(style="display:inline-block; margin-left: 0.2em", numericInput("q", label=NULL, 0.999,min = 0.9, max = 1, width= '100%'))
                   )),
                   column(6, 
                     conditionalPanel(
                     condition = "input.q2f == true",
                     div(style="display:inline-block", strong("q2:")),
                     div(style="display:inline-block; margin-left: 0.75em", numericInput("q2", label=NULL, 0.999,min = 0.9, max = 1, width= '100%'))
                   ),
                   conditionalPanel(
                     condition = "input.q3f == true",
                     div(style="display:inline-block", strong("q3:")),
                     div(style="display:inline-block; margin-left: 0.75em", numericInput("q3", label=NULL, 0.999,min = 0.9, max = 1, width= '100%'))
                   )
                   )
                   ),
                   
                   fluidRow(column(6,conditionalPanel(
                                        condition = "input.taskType > 4",
                                        checkboxInput("cd1f", "Fix fraction of first reaction cd1 or define its limits:")
                                  )),
                            column(6, conditionalPanel(
                              condition = "input.taskType > 4&&input.cd1f == false",
                              sliderInput("cd1Lim",label=NULL,min = -3,max = 4,value = c(0.01,0.99), step=0.01)
                            ),
                                      conditionalPanel(
                                        condition = "input.cd1f == true&&input.taskType > 4",
                                        div(style="display:inline-block", strong("cd1:")),
                                        div(style="display:inline-block; margin-left: 0.2em", numericInput("cd1fixed", label=NULL, 0.1,min = 0, max = 1, width= '100%'))
                            ))),
                   conditionalPanel(
                     condition = "input.taskType ==9",
                        fluidRow(column(6,div(style="display:inline-block", strong("cr21:")),
                                          div(style="display:inline-block; margin-left: 0.2em", numericInput("cr21", label=NULL, 1.0,min = 0.01, max = 20, width= '100%'))
                                ),
                                column(6, div(style="display:inline-block", strong("cr22:")),
                                          div(style="display:inline-block; margin-left: 0.75em", numericInput("cr22", label=NULL, 1.0,min = 0.01, max = 20, width= '100%'))
                          ))),
                   conditionalPanel(
                     condition = "input.taskType ==10",
                     div(style="display:inline-block", strong("sw(a):")),
                     div(style="display:inline-block; margin-left: 0.75em", numericInput("swA", label=NULL, 0.5,min = 0.01, max = 0.99, width= '100%'))
                   ),
                   hr(),
                   conditionalPanel(
                     condition = "(input.taskType !=2&&input.taskType !=3)&&!(input.taskType==1&&input.mf1 == true&&input.m1fixed==10)",
                     h4("Flexible single step (ePT):"),
                     withMathJax(),
                     div("$$ \\Large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k\\left( 1 - \\alpha \\right)^n \\left[ 1 - q \\left( 1 - \\alpha \\right) \\right]^m
                 }
                 $$"),
                     p(h6("1. Burnham, A. K.; Zhou, X.; Broadbelt, L. J. Critical Review of the Global Chemical Kinetics of Cellulose Thermal Decomposition. Energy Fuels 2015, 29 (5), 2906–2918. doi:10.1021/acs.energyfuels.5b00350."))
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"F0\"",
                     h4("Zero-order reaction:"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k
                 }
                 $$") ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"F1\"",
                     h4("First-order reaction:"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)
                 }
                 $$") ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"F2\"",
                     h4("Second-order reaction:"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)^2
                 }
                 $$") ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"F3\"",
                     h4("Third-order reaction:"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)^3
                 }
                 $$") ),
                   conditionalPanel(
                     condition = "(input.taskType == 2 && (input.reMod == \"A2\"||input.reMod == \"A3\"||input.reMod == \"A4\"))||(input.taskType != 7&&input.taskType != 10&&input.mf1 == true&&input.m1fixed==10)",
                     h4("KJMAE nucleation-growth (An, n = 2, 3, 4):"),
                     withMathJax(),
                     div("$$ \\large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k n( 1- \\alpha) [-ln(1-\\alpha)]^{(n-1)/n}
                 }
                 $$"),
                     p(h6("1.	Kolmogorov A. A statistical theory for the recrystallization of metals. Izv Akad Nauk SSSR Ser Mat 1937:355–9.", br(),
                          "2. Johnson WA, Mehl RF. Reaction kinetics in processes of nucleation and growth. Trans AIME 1939;135:416–42.",br(),
                          "3. Avrami M. Kinetics of Phase Change. I General Theory. J Chem Phys 1939;7:1103. doi:10.1063/1.1750380.",br(),
                          "4. Erofeev BV. Generalized Equation of Chemical Kinetics and its Application to Reactions involving solid phase components. Dokl Akad Nauk USSR 1946;52:515–8."
                     ))
                   ),
                   
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"R2\"",
                     h4("Contracting cylinder (R2):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 2k ( 1- \\alpha)^{1/2}
                 }
                 $$")
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"R3\"",
                     h4("Contracting sphere (R3):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 3k ( 1- \\alpha)^{2/3}
                 }
                 $$")
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && (input.reMod == \"P2\"|| input.reMod == \"P3\"||input.reMod == \"P4\"|| input.reMod == \"P23\")",
                     h4("Power law (Pn, n = 2/3, 2, 3, 4):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k n \\alpha^{(n-1)/n}
                 }
                 $$")
                   ), #------------------------------------------add reference here
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"D1\"",
                     h4("One-dimensional diffusion (D1):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{1}{2\\alpha}
                 }
                 $$")
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"D2\"",
                     h4("Two-dimensional diffusion (D2):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{1}{-ln(1-\\alpha)}
                 }
                 $$")
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"D3\"",
                     h4("3D Jander diffusion (D3):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{3(1-\\alpha)^{2/3}}{2[1-(1-\\alpha)^{1/3}]}
                 }
                 $$"),
                     p(h6("1.Jander W. Reaktionen im festen Zustande bei höheren Temperaturen. Reaktionsgeschwindigkeiten endotherm verlaufender Umsetzungen. Z Für Anorg Allg Chem 1927;163:1–30. doi:10.1002/zaac.19271630102."))
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"D4\"",
                     h4("3D Ginstling-Brounshtein diffusion (D4):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{3}{2[(1-\\alpha)^{-1/3}-1]}
                 }
                 $$"),
                     p(h6("1.Ginstling AM, Brounshtein BI. The diffusion kinetics of reactions in spherical particles. J Appl Chem USSR 1950;23:1249–59."))
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"L2\"",
                     h4("Random scission of polymer chain (PRS2):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 2k (\\alpha^{1/2}-\\alpha)
                 }
                 $$"),
                     p(h6("1. Flynn JH, Wall LA. General treatment of the thermogravimetry of polymers. J Res Natl Bur Stand Sect Phys Chem 1966;70A:487. doi:10.6028/jres.070A.043.",br(),
                          "2. Simha R, Wall LA. Kinetics of Chain Depolymerization. J Phys Chem 1952;56:707–15. doi:10.1021/j150498a012."))
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 2 && input.reMod == \"B1\"",
                     h4("Classical Prout-Tompkins (B1):"),
                     withMathJax(),
                     div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\alpha (1-\\alpha)
                 }
                 $$"),
                     p(h6("1.Prout EG, Tompkins FC. The thermal decomposition of potassium permanganate. Trans Faraday Soc 1944;40:488. doi:10.1039/tf9444000488."))
                   ),
                   conditionalPanel(
                     condition = "input.taskType == 3||input.taskType == 8",
                     h4("Autocatalysis (DMM):"),
                     withMathJax(),
                     div("$$ \\large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k_1 \\left( 1 - \\alpha \\right) + k_2 \\left( 1 - \\mu \\right) \\frac{\\alpha \\left(1 - \\alpha \\right)}{1 - \\mu \\alpha }
                 }
                 $$"),
                     p(h6("1. Dubovitskii, F. I.; Manelis, G. B.; Merzhanov, A. G. Formal Kinetic Model of Thermal Decomposition of Explosives in Liquid State. Trans Acad. Nauk SSSR Dokl. 1958, 121 (4), 668–670."))
                   ),
                   hr()
                   )
        ),
        
        img(src="thinks.png", height = 75, width = 100),#, align='right'),
        h5(strong("Dr. Nikita V. Muravyev, Version 11.05.21")),
        p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
             a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
             a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
        )
      )
    )
  )
))
