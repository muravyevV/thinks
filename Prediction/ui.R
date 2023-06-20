#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: https://doi.org/10.3390/molecules24122298
#This routine performs the predition of thermal behavior using the model-fitting or the isoconversional kinetic data.
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
                    windowTitle = "THINKS/Prediction",
                    title = tags$head(tags$link(rel="icon", 
                                                href="favicon.ico", 
                                                type="image/x-icon")
                    )),
                  h2("Kinetic prediction", style="margin-bottom: 0em"),
                  br(),
  
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
                              margin-bottom: 0.8em;
                    }
                        .radioDiff {
                              margin-left: 0em;
                              padding: 0em;
                              font-size: 9px;
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
                 h5("Temperature program"),
                 radioButtons("Tprg", label=NULL, c("linear heating"=1,"isothermal"=2, "user-specified"=3), 
                                selected =1, inline =T),
                 conditionalPanel(
                   condition = "input.Tprg == 1",
                   fluidRow(column(4, helpText("T(initial), °C:"),
                                      numericInput("THRini", label=NULL, 35, min = -500, max = 1500, step=50, width='100%')),
                            column(4, helpText("T(final), °C:"),
                                      numericInput("THRfin", label=NULL, 300, min = 35, max = 3500, step=50, width='100%')),
                            column(4, helpText("HR, K/min:"),
                                   numericInput("THR", label=NULL, 1, min = 1E-9, max = 1E+9, step=1, width='100%'))
                                   )
                   ),
                 conditionalPanel(
                   condition = "input.Tprg == 2",
                   fluidRow(
                     column(7,helpText("Isothermal temperature, °C:")),
                     column(5,numericInput("Tiso", label=NULL, 100, min = 35, max = 3500, step=10, width='100%'))),
                   fluidRow(
                     column(7,helpText("Duration of experiment, min:")),
                     column(5,numericInput("Timefin", label=NULL, 300, min = 5, max = 1E+9, step=5, width='100%')))
                   ),
                 conditionalPanel(
                   condition = "input.Tprg == 3",
                   fileInput('file2', label="Load .txt file*", accept=c('text/plain')),
                   p("*Data file structure: row 1 - T[°C], row 2 - t[s], row 3 (optional) - alpha[-], row 4 (optional) - da/dt[1/s], Separator: space", style="font-size: 70%"),
                   checkboxInput("showAlTprg", "Show conversion data from file", value=F)
                 ),
                
                 hr(),
                 h5("Basic plot settings"),
                 fluidRow(
                   column(3,
                          helpText("X-axis: "),
                          helpText("Y-axis: ")),
                   column(9,
                          radioButtons("showTP", label=NULL, c("Time"=2, "Temperature"=3), selected =2, inline =T),
                          radioButtons("showY", label=NULL, c("Conversion"=1,"Conversion rate"=2),selected =1, inline =T))
                   ),
                 fluidRow(
                   column(4,actionButton("go", "Calculate", class="btn btn-primary")),
                   column(5, conditionalPanel(
                     condition = "input.go != 0",
                          downloadButton('downloadData2', 'Export prediction', class="btn btn-success")))
                 ),
                 br()
            ),
            tabPanel("Plot",
                     br(),
                     div(style="display:inline-block", helpText("Temperature in: ")),
                     div(style="display:inline-block; margin-left: 0.5em", radioButtons("PlotT", label=NULL, c("C"=1,"K"=2),selected =1, inline =T)),
                     checkboxInput("showTt", "Show T(t) curve", value=F),
                     checkboxInput("Compact", "Compact mode (adjusted plot heights)", value=T),
                     conditionalPanel(
                       condition = "input.Compact == false", 
                       fluidRow(
                         column(7,helpText("Basic height (px):")),
                         column(5,numericInput("BaseHeight", label=NULL, 350, min = 100, max = 1000, step=50, width='100%'))
                       )),
                     radioButtons("TitleFont", label="Title font style",c("normal"=1,"bold"=2,"italic"=3, "bold italic"=4), 
                                  selected =2, inline =T),
                     helpText("Font size (cex):"),
                     fluidRow(column(4, helpText("Title font:"),
                                     numericInput("TitleSize", label=NULL, 1.2, min = 0.8, max = 2.5, step=0.1, width='100%')),
                              column(4, helpText("Axis font:"),
                                     numericInput("AxisSize", label=NULL, 1.7, min = 0.8, max = 2.5, step=0.1, width='100%')),
                              column(4, 
                                     helpText("Legend font:"),
                                     numericInput("LegendSize", label=NULL, 1.3, min = 0.5, max = 2.5, step=0.1, width='100%'))
                     ),
                     
                     conditionalPanel(
                       condition = "input.KinMod >= 4&&input.KinT == 1",
                       helpText("Change the plot with concentrations:"),
                       fluidRow(
                         column(6,selectInput("ConcLeft", "On left Y-axis", c("A"=0,"B"=1,"C"=2), 
                                              selected = c(0,1,2), multiple = T)),
                         column(6,selectInput("ConcRight", "On right Y-axis", c("A"=0,"B"=1,"C"=2), 
                                              selected = NULL, multiple = T, width='100%'))
                       ))
            ),
            tabPanel("Advanced",  
                  br(),
                     #h5("Computational settings"),
                     fluidRow(
                       column(7,helpText("Number of points (limited in online version to 2000):")),
                       column(5,numericInput("Npts", label=NULL, 1000, min = 200, max = 2000, step=50, width='100%'))
                     ),
                     fluidRow(
                       column(7,helpText("Starting value of conversion degree:")),
                       column(5,numericInput("a0",label=NULL,1E-4, min = 1E-9,max = 0.2, width='100%'))
                     ),
                     div(style="display:inline-block", helpText("Time step: ")),
                     div(style="display:inline-block; margin-left: 0.2em", radioButtons("TimeSc", label=NULL, c("linear"=1,"logarithmic"=2),selected =1, inline =T)),
                     conditionalPanel(
                       condition = "input.KinT == 1", 
                       helpText("ODE solvers [order in brackets]:"),
                       div(class="radioDiff",radioButtons("methodDiff", label=NULL, c("Euler [1]"=1,"Runge-Kutta [2]"=2,"Runge-Kutta [4]"=4,
                                                                                      #"Runge-Kutta Cash-Karp [4(5)]"=7, 
                                                                                      "Dormand-Prince [4(5) local 7]"=5, 
                                                                                      "Petzold-Hindmarsh, lsoda"=6), 
                                    selected =4, inline =T)),
                       checkboxInput("suppression", "Suppress NaN's and outliers after solver", value=T),
                       checkboxInput("CompSolv", "Compare the result with another ODE solver:", value=F)
                     ),
                     conditionalPanel(
                       condition = "input.KinT == 1&&input.CompSolv == true",     
                       div(class="radioDiff",radioButtons("ODEs", label=NULL, c("Euler [1]"=1,"Runge-Kutta [2]"=2,"Runge-Kutta [4]"=4,
                                                                                      #"Runge-Kutta Cash-Karp [4(5)]"=7, 
                                                                                      "Dormand-Prince [4(5) local 7]"=5, 
                                                                                      "Petzold-Hindmarsh, lsoda"=6), 
                                                          selected =1, inline =T)))
                  
            )
        ),
        hr(),
        img(src="thinks.png", height = 75, width = 100),#, align='right'),
        h5(strong("Dr. Nikita V. Muravyev, Version 22.09.20")),
        p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
             a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
             a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
        )
    ),
    
    mainPanel(
      column(width=7,
        tabsetPanel(id = "MainTabset", selected = "panel2",
           tabPanel("Results", value = "panel1",   
             conditionalPanel(
               condition = "input.go != 0",
               plotOutput("Plot3", height = "auto")),
             conditionalPanel(
               condition = "input.go != 0&&input.KinT==2",
               plotOutput("Plot1", height= "auto")),
             conditionalPanel(
               condition = "input.go != 0&&input.KinMod >= 4&&input.KinT == 1", 
               plotOutput("Plot4",height = "auto")),
             conditionalPanel(
               condition = "input.go != 0&&input.KinT == 1&&input.CompSolv == true",
               plotOutput("Plot5",height = "auto"))
           ),
           tabPanel("Help", value="panel2",
               div(class = "initialParent", 
                   div( class="initialChild",
                        p(strong("/ Instructions /")),
                        p("1.    Specify the details of target thermal program in the left block (linear heating, isothermal, 
                        or that loaded from the file). Note, that the prediction is purely kinetic, 
                         heat accumulation, heat loss phenomena are not considered specifically."), 
                        p("2.    Select the type of kinetic parameters that will be used, i.e., the model fitting results or 
                        isoconversional data. The model fitting data are defined by user (e.g., obtained via ",
                          tags$a(href = "http://www.thinks.chemphys.ru/shiny/Model_Linear/", target="_blank", "Model fitting (linear)"), ", ",
                          tags$a(href = "http://www.thinks.chemphys.ru/shiny/Model_Nonlinear/", target="_blank", "Model fitting (nonlinear)"), 
                          " routines), the isoconverional results are loaded from the file that is exported from ",
                          tags$a(href = "http://www.thinks.chemphys.ru/shiny/Isoconversional/", target="_blank", "Isoconversional"), 
                          "routine)."),
                        p("3.    Once the results are plotted, use the \"Plot\" tab to optimize its graphical representation. 
                        If necessary tune the advanced parameters of computation (number of points, the type of solver used for ODEs) in tab \"Advanced\". 
                        Note, that the number of points and the time step type (logarithmic one can be beneficial for long 
                        isothermal experiments) affect only the results for first three solvers,
                        the high-accuracy solvers (from \"rk45dp7\") use their own selection of steps. Note that some stiff models (e.g., v5)
                        require using of the strong solvers (e.g., \"lsoda\"). The accuracy of various solvers can be compared once the relevant checkbox
                        is set and two solvers are selected.
                        Please, keep in mind, that increasing the computational complexity reduces the server resources 
                        available the other users."),
                        p("4.    Export the predicted data in text file to work with."),
                        p("5.    If ready to start, press \"Calculate\" and have a good day!")
                   ))
           )
           
           ) ),
      
      column(width=5,           
             h5("Kinetic parameters", style="margin-top: 0em"),  
              radioButtons("KinT", NULL, c("Model fitting"=1,"Isoconversional"=2), selected =2, inline =T),
              conditionalPanel(
                condition = "input.KinT == 1",
                selectInput("KinMod", "Kinetic model:", c("v1: Flexible single step (ePT)"=1,
                                                          "v2: Simple reaction types"=2,
                                                          "v3: Autocatalysis (DMM)"=3,
                                                          "v4: Parallel ePT || ePT"=4,
                                                          "v5: Consecutive ePT -> ePT"=5,
                                                          "v6: Independent ePT + ePT"=6,
                                                          "v7: ePT sw(a) ePT"=7,
                                                          "v8: Consecutive DMM -> ePT"=8,
                                                          "v9: ePT -> ePT || ePT"=9
                ), selected = 1)),
              conditionalPanel(
                condition = "input.KinMod == 2&&input.KinT == 1",
                selectInput("reMod", "Simple reaction model:", c("Zero-order reaction (F0)" = "F0","First-order reaction (F1)" = "F1",
                                                                 "Second-order reaction (F2)" = "F2","Third-order reaction (F3)" = "F3",
                                                                 "KJMAE nucleation-growth (A2)" = "A2","KJMAE nucleation-growth (A3)" = "A3","KJMAE nucleation-growth (A4)" = "A4",
                                                                 "Contracting cylinder (R2)" = "R2","Contracting sphere (R3)" = "R3",
                                                                 "Power law (P2)" = "P2","Power law (P3)" = "P3","Power law (P4)" = "P4","Power law (P23)" = "P23",
                                                                 "One-dimensional diffusion (D1)" = "D1","Two-dimensional diffusion (D2)" = "D2",
                                                                 "3D Jander diffusion (D3)" = "D3","3D Ginstling-Brounshtein diffusion (D4)" = "D4",
                                                                 "Polymer random scission (PRS2)" = "L2","classical Prout-Tompkins (B1)"="B1"), selected = "F1")
              ),
              
              conditionalPanel(
                condition = "input.KinT == 1",   
                tabsetPanel(id = "tabs", 
                          tabPanel("First stage", 
                                     p(style="margin-left: 0.2em;margin-top: 0.5em",
                                       div(style="display:inline-block", helpText("Activation energy Ea1, kJ/mol:")),
                                       div(style="display:inline-block; margin-left: 0.2em;", numericInput("Ea", label=NULL, 170,min = 50, max = 350, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                       div(style="display:inline-block", helpText("Preexponent lnA1, 1/s:")),
                                       div(style="display:inline-block; margin-left: 2.8em;", numericInput("lA1",label=NULL, 40, min = 1,max = 100, width='100%'))
                                     ),
                                   conditionalPanel(
                                     condition = "input.KinMod != 2&&input.KinMod != 3",
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Reaction order n1:")),
                                      div(style="display:inline-block; margin-left: 4.4em", numericInput("n1", label=NULL, 1,min = 0, max = 3, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Exponent m1:")),
                                      div(style="display:inline-block; margin-left: 6.2em", numericInput("m1", label=NULL, 0,min = -3, max = 3, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Initiation parameter q1:")),
                                      div(style="display:inline-block; margin-left: 2.6em", numericInput("q", label=NULL, 1,min = 0.9, max = 1, width= '100%'))
                                   )),
                                   conditionalPanel(
                                     condition = "input.KinMod >= 5&&input.KinT == 1",
                                     p(style="margin-left: 0.2em;margin-top: -1em;;margin-bottom: -1em",
                                       div(style="display:inline-block", helpText("Fraction of first reaction cd1:")),
                                       div(style="display:inline-block; margin-left: 0.9em", numericInput("cd1", label=NULL, 0.5,min = -3, max = 3, width='100%'))
                                     )),
                                   conditionalPanel(
                                     condition = "input.KinMod == 4&&input.KinT == 1",
                                     p(style="margin-left: 0.2em;margin-top: -1em;;margin-bottom: -1em",
                                       div(style="display:inline-block", helpText("Fraction of first reaction cr1:")),
                                       div(style="display:inline-block; margin-left: 0.9em", numericInput("cr1", label=NULL, 1.0,min = -3, max = 3, width='100%'))
                                     ))
                                   ),
                          tabPanel("Second stage", 
                                   conditionalPanel(
                                     condition = "input.KinMod >= 3", 
                                     p(style="margin-left: 0.2em;margin-top: 0.5em",
                                      div(style="display:inline-block", helpText("Activation energy Ea2, kJ/mol:")),
                                      div(style="display:inline-block; margin-left: 0.2em", numericInput("Ea2", label=NULL, 180,min = 50, max = 350, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Preexponent lnA2, 1/s:")),
                                      div(style="display:inline-block; margin-left: 2.8em", numericInput("lA2",label=NULL, 40, min = 1,max = 100, width='100%'))
                                   )),
                                   conditionalPanel(
                                     condition = "input.KinMod == 3",
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("mu-parameter:")),
                                      div(style="display:inline-block; margin-left: 5.5em", numericInput("mu", label=NULL, 0.5,min = 0, max = 1, width='100%'))
                                   )),
                                   conditionalPanel(
                                     condition = "input.KinMod >= 4",
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Reaction order n2:")),
                                      div(style="display:inline-block; margin-left: 4.4em", numericInput("n2", label=NULL, 0.3,min = 0, max = 3, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Exponent m2:")),
                                      div(style="display:inline-block; margin-left: 6.2em", numericInput("m2", label=NULL, 1,min = -3, max = 3, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                      div(style="display:inline-block", helpText("Initiation parameter q2:")),
                                      div(style="display:inline-block; margin-left: 2.6em", numericInput("q2", label=NULL, 1,min = 0.9, max = 1, width= '100%'))
                                   )),
                                   conditionalPanel(
                                     condition = "input.KinMod == 4&&input.KinT == 1",
                                     p(style="margin-left: 0.2em;margin-top: -1em;;margin-bottom: -1em",
                                       div(style="display:inline-block", helpText("Fraction of second reaction cr2:")),
                                       div(style="display:inline-block; margin-left: 0.9em", numericInput("cr2", label=NULL, 1.0,min = -3, max = 3, width='100%'))
                                     )),
                                   conditionalPanel(
                                     condition = "input.KinMod == 9&&input.KinT == 1",
                                     p(style="margin-left: 0.2em;margin-top: -1em;;margin-bottom: -1em",
                                       div(style="display:inline-block", helpText("Fraction of second reaction cr21:")),
                                       div(style="display:inline-block; margin-left: 0.9em", numericInput("cr21", label=NULL, 1.17,min = -3, max = 3, width='100%'))
                                     ))
                                   ),
                          tabPanel("Third stage", 
                                   conditionalPanel(
                                     condition = "input.KinMod >= 8", 
                                     p(style="margin-left: 0.2em;margin-top: 0.5em",
                                       div(style="display:inline-block", helpText("Activation energy Ea3, kJ/mol:")),
                                       div(style="display:inline-block; margin-left: 0.2em", numericInput("Ea3", label=NULL, 192.8,min = 50, max = 350, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                       div(style="display:inline-block", helpText("Preexponent lnA3, 1/s:")),
                                       div(style="display:inline-block; margin-left: 2.8em", numericInput("lA3",label=NULL, 40, min = 1,max = 100, width='100%'))
                                     )),
                                   conditionalPanel(
                                     condition = "input.KinMod >= 8",
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                       div(style="display:inline-block", helpText("Reaction order n3:")),
                                       div(style="display:inline-block; margin-left: 4.4em", numericInput("n3", label=NULL, 0.25,min = 0, max = 3, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                       div(style="display:inline-block", helpText("Exponent m3:")),
                                       div(style="display:inline-block; margin-left: 6.2em", numericInput("m3", label=NULL, 0,min = -3, max = 3, width='100%')),
                                     ),
                                     p(style="margin-left: 0.2em;margin-top: -1em",
                                       div(style="display:inline-block", helpText("Initiation parameter q3:")),
                                       div(style="display:inline-block; margin-left: 2.6em", numericInput("q3", label=NULL, 1,min = 0.9, max = 1, width= '100%'))
                                     )),
                                   conditionalPanel(
                                     condition = "input.KinMod == 9&&input.KinT == 1",
                                     p(style="margin-left: 0.2em;margin-top: -1em;;margin-bottom: -1em",
                                       div(style="display:inline-block", helpText("Fraction of first reaction cr22:")),
                                       div(style="display:inline-block; margin-left: 0.9em", numericInput("cr22", label=NULL, 0.98,min = -3, max = 3, width='100%'))
                                     ))
                          )
              )),
             
             conditionalPanel(
               condition = "input.KinT == 2",
               helpText("Load results of the ", tags$a(href = "http://www.thinks.chemphys.ru/shiny/Isoconversional/", target="_blank", "isoconversional")," analysis"),
               fileInput('file1', label=NULL, accept=c('text/plain')),
               fluidRow(
                  column(5,numericInput("aMin", "Use conversion from:", 0.01,min = 0.005, max = 0.3)),
                  column(5,numericInput("aMax", "Use conversion up to:", 0.99,min = 0.7, max = 1))
               ),
               radioButtons("KinFT", "Spline type for kinetic function:", c("fmm"=1,"natural"=2, 
                                                                 "linear"=3), selected =3, inline =T)),
            
          conditionalPanel(condition = "input.KinT == 1", hr()
          ),
          
             conditionalPanel(
               condition = "input.KinMod != 3 && input.KinMod != 2&&input.KinT == 1&&!(input.KinMod==1&&input.m1==10)",
               h5("Flexible single step (ePT):"),
               withMathJax(),
               div("$$ \\Large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k\\left( 1 - \\alpha \\right)^n \\left[ 1 - q \\left( 1 - \\alpha \\right) \\right]^m
                 }
                 $$"),
               p(h6("1. Burnham, A. K.; Zhou, X.; Broadbelt, L. J. Critical Review of the Global Chemical Kinetics of Cellulose Thermal Decomposition. Energy Fuels 2015, 29 (5), 2906–2918. doi:10.1021/acs.energyfuels.5b00350."))
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"F0\"&&input.KinT == 1",
               h5("Zero-order reaction:"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k
                 }
                 $$") ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"F1\"&&input.KinT == 1",
               h5("First-order reaction:"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)
                 }
                 $$") ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"F2\"&&input.KinT == 1",
               h5("Second-order reaction:"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)^2
                 }
                 $$") ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"F3\"&&input.KinT == 1",
               h5("Third-order reaction:"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k ( 1- \\alpha)^3
                 }
                 $$") ),
             conditionalPanel(
               condition = "(input.KinMod == 2 && (input.reMod == \"A2\"||input.reMod == \"A3\"||input.reMod == \"A4\")&&input.KinT == 1) ||(input.m1==10&&input.KinT == 1&&(input.KinMod==1||input.KinMod==4||input.KinMod==5||input.KinMod==6||input.KinMod==9))",
               h5("KJMAE nucleation-growth (An, n = 2, 3, 4):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
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
               condition = "input.KinMod == 2 && input.reMod == \"R2\"&&input.KinT == 1",
               h5("Contracting cylinder (R2):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 2k ( 1- \\alpha)^{1/2}
                 }
                 $$")
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"R3\"&&input.KinT == 1",
               h5("Contracting sphere (R3):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = 3k ( 1- \\alpha)^{2/3}
                 }
                 $$")
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && (input.reMod == \"P2\"|| input.reMod == \"P3\"||input.reMod == \"P4\"|| input.reMod == \"P23\")&&input.KinT == 1",
               h5("Power law (Pn, n = 2/3, 2, 3, 4):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k n \\alpha^{(n-1)/n}
                 }
                 $$")
             ), #------------------------------------------add reference here
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"D1\"&&input.KinT == 1",
               h5("One-dimensional diffusion (D1):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{1}{2\\alpha}
                 }
                 $$")
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"D2\"&&input.KinT == 1",
               h5("Two-dimensional diffusion (D2):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{1}{-ln(1-\\alpha)}
                 }
                 $$")
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"D3\"&&input.KinT == 1",
               h5("3D Jander diffusion (D3):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{3(1-\\alpha)^{2/3}}{2[1-(1-\\alpha)^{1/3}]}
                 }
                 $$"),
               p(h6("1.Jander W. Reaktionen im festen Zustande bei höheren Temperaturen. Reaktionsgeschwindigkeiten endotherm verlaufender Umsetzungen. Z Für Anorg Allg Chem 1927;163:1–30. doi:10.1002/zaac.19271630102."))
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"D4\"&&input.KinT == 1",
               h5("3D Ginstling-Brounshtein diffusion (D4):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\frac{3}{2[(1-\\alpha)^{-1/3}-1]}
                 }
                 $$"),
               p(h6("1.Ginstling AM, Brounshtein BI. The diffusion kinetics of reactions in spherical particles. J Appl Chem USSR 1950;23:1249–59."))
             ),
             conditionalPanel(
               condition = "input.KinMod == 2 && input.reMod == \"L2\" &&input.KinT == 1",
               h5("Random scission of polymer chain (PRS2):"),
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
               condition = "input.KinMod == 2 && input.reMod == \"B1\" &&input.KinT == 1",
               h5("classical Prout-Tompkins (B1):"),
               withMathJax(),
               div("$$ \\LARGE \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k \\alpha (1-\\alpha)
                 }
                 $$"),
               p(h6("1.Prout EG, Tompkins FC. The thermal decomposition of potassium permanganate. Trans Faraday Soc 1944;40:488. doi:10.1039/tf9444000488."))
             ),
             conditionalPanel(
               condition = "input.KinMod == 3&&input.KinT == 1",
               h5("Autocatalysis (DMM):"),
               withMathJax(),
               div("$$ \\Large \\bbox[white] 
                 {
                 \\frac{d\\alpha}{dt} = k_1 \\left( 1 - \\alpha \\right) + k_2 \\left( 1 - \\mu \\right) \\frac{\\alpha \\left(1 - \\alpha \\right)}{1 - \\mu \\alpha }
                 }
                 $$"),
               p(h6("1. Dubovitskii, F. I.; Manelis, G. B.; Merzhanov, A. G. Formal Kinetic Model of Thermal Decomposition of Explosives in Liquid State. Trans Acad. Nauk SSSR Dokl. 1958, 121 (4), 668–670."))
             )
             
             )
             )
             )
             ))