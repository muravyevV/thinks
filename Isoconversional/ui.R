#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine performs the isoconversional analysis.
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
    windowTitle = "THINKS/Isoconversional",
    title = tags$head(tags$link(rel="icon", 
                                href="favicon.ico", 
                                type="image/x-icon")
    )),
  h2("Isoconversional analysis", style="margin-bottom: 0em"), br(),
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
                    hr { 
                          margin-top: 0em !important;
                          margin-bottom: 1em !important;
                          border-width: 1px;
                        } "))
    ),

  sidebarLayout(
    sidebarPanel(width=3,
            tabsetPanel(
               tabPanel("General",
                 br(),       
                 #h4("Input data"),
                 helpText("Load the .txt prepared with ", tags$a(href = "http://www.thinks.chemphys.ru/shiny/Import/", target="_blank", "Import")," app"),
                 fileInput('file1', label=NULL, accept=c('text/plain')),
                  helpText("Data file structure: row i - T[°C], row i+1 - t[s], row i+2 - α[0..1], row i+3 - dα/dt[1/s]. 
                        Separator: space"),
                 div(style="display:inline-block", helpText("Select method: ")),
                 div(style="display:inline-block; margin-left: 0.2em", radioButtons("Method", label=NULL, 
                                        c("Friedman"=1, "Vyazovkin"=2
                                          ), selected =1, inline =T)), 
                 fluidRow(column(4, helpText("α(min):"),
                                 numericInput("aMin", label=NULL, 0.05, min = 0.001, max = 0.5, step=0.01)
                                 ),
                          column(4, helpText("α(max):"),
                                 numericInput("aMax", label=NULL, 0.95, min = 0.5, max = 0.999, step=0.01)
                                ),
                          column(4, helpText("Δα:"),
                                 numericInput("aStep", label=NULL, 0.05, min = 0.005, max = 0.1, step=0.01)
                                )
                          ),         
                 hr(),
                 helpText("Show additionally:"),
                 div(style="display:inline-block; margin-left: 0.2em",checkboxInput("showR2", "R^2(Ω) plot", value=T)),
                 div(style="display:inline-block",checkboxInput("showKiss", "Kissinger result", value=F)), 
                 checkboxInput("showBoth", "Plot the results of both isoconversional methods", value=F), 
                 div(style="display:inline-block",checkboxInput("showFri", "Friedman plot", value=F)), br(),
                 conditionalPanel(
                   condition = "input.showFri == true",
                   #radioButtons("FriD", label=NULL, c("2D"=1,"3D"=2), selected =1, inline =T),
                   div(style="display:inline-block","Show every: "),
                   div(style="display:inline-block",numericInput("daX", label=NULL, 3, min = 1, max = 10, step=1)),
                   div(style="display:inline-block"," conversion")
                   ),
              #   conditionalPanel(
              #     condition = "input.go != 0",
              #     radioButtons("KinFT", "Kinetic function type:", c("spline-fmm"=1,"spline-natural"=2, 
              #                                                       "linear"=3), selected =2, inline =T)
              #   )
                fluidRow(
                  column(4,actionButton("go", "Calculate", class="btn btn-primary")),
                  column(5, conditionalPanel(
                    condition = "input.go != 0",
                    downloadButton('downloadData1', 'Export results', class="btn btn-success")))
                )
              ),
            tabPanel("Plot settings",
                     br(),
                     checkboxInput("Compact", "Compact mode (adjusted plot heights)", value=T),
                     conditionalPanel(
                       condition = "input.Compact == false", 
                       fluidRow(
                         column(7,helpText("Basic height (px):")),
                         column(5,numericInput("BaseHeight", label=NULL, 350, min = 100, max = 1000, step=50, width='100%'))
                       )),
                    # radioButtons("TitleFont", label="Title font style",c("normal"=1,"bold"=2,"italic"=3, "bold italic"=4), 
                    #              selected =1, inline =T),
                    # numericInput("TitleSize", label="Title font size", 1.4, min = 0.8, max = 2.5, step=0.1, width='100%')
                     fluidRow(column(4, numericInput("AxisSize", label="Axis title size", 1.8, min = 0.8, max = 2.5, step=0.1, width='100%')),
                              column(4, numericInput("LegendSize", label="Labels size", 1.1, min = 0.5, max = 2.5, step=0.1, width='100%')),
                              column(4, numericInput("LegendInset", label="Inset(legend)", -0.3, min = -10, max = 10, width='100%'))
                     ),
                     fluidRow(column(4, numericInput("marginTitle", label="Margin(axis)", 0, min = 0, max = 6, step=0.1, width='100%')),
                              column(4, numericInput("marginLeft", label="Margin(left)", 0, min = 0, max = 6, step=0.1, width='100%')),
                              column(4, numericInput("marginBottom", label="Margin(bottom)", 0, min = 0, max = 6, step=0.1, width='100%')))
                     )
            )
    ),
    
    mainPanel(
      
      column(width=7,
             tabsetPanel(id = "MainTabset", selected = "panel2",
                         tabPanel("Results", value = "panel1",
                                  conditionalPanel(
                                    condition = "input.go != 0",
                                    plotOutput("Plot1", height="auto")),
                                  conditionalPanel(
                                    condition = "input.showKiss == true",
                                    verbatimTextOutput("textKiss")),
                                  conditionalPanel(
                                    condition = "input.go != 0",
                                    plotOutput("Plot2", height="auto"))
                         ),
                         tabPanel("Help", value="panel2",
                                  div(class = "initialParent", 
                                      div( class="initialChild",
                        p(strong("/ Instructions /")), 
                        p("1. Load thermal analysis data. The input file can be created with ", 
                        tags$a(href = "http://www.thinks.chemphys.ru/shiny/Import/", target="_blank", "import"), 
                        " app or manually. It should be the single file with the rows of equal length, four rows for each experiment, i.e.,
                        (i) temperature [°C], (i+1) time [s], (i+2) conversion degree [0..1], (i+3) conversion rate [1/s], 
                        , separator - space. Example can be found ",
                        tags$a(href = "example_nano-Ti_tga.txt", target="_blank", "here"),
                        
                        ". Note, that more than four experiments is needed to perform the kinetic analysis (see the ",
                                tags$a(href = "https://doi.org/10.1016/j.tca.2011.03.034", target="_blank", "ICTAC Kinetic committee recommendations"),
                                ")."),
                        p("2. Press \"Calculate\" button. The range of conversion degree values that used in calculations can be 
                        adjusted defining α-min, α-max and the step dα (by default the 0.05..0.95 range is set to eliminate the minor effect of the baseline selection). 
                        Check the accuracy of the isoconversional analysis - 
                        correlation coefficient (R^2) in case of the Friedman method or the optimization indicator (Ω) for the advanced Vyazovkin method.
                        Note, that while the first value has to be close to unity, the latter should approach the n(n-1) value (where n is the number or experiments).
                         Additionally, the results of the Kissinger method can be shown on the plot and in text from
                          (consider it only as a first estimation, not the as final result!). The results of both isoconversional methods, i.e., the differential Friedman
                          and integral advanced Vyazovkin, can be plotted together for visual comparison."),
                        p("3. Finally, use the \"Export results\" button
                          to generate the text file with the results of the selected isoconversional method. This file can be used later for thermal ", 
                          tags$a(href = "http://www.thinks.chemphys.ru/shiny/Prediction/", target="_blank", "prediction."))
               ))))
               ),
      
      column(width=5,
             conditionalPanel(
               condition = "input.showFri ==true",
               plotOutput("Plot4", height="auto")),
             conditionalPanel(
               condition = "input.showKiss ==true",
               plotOutput("Plot5", height="auto")),
             withMathJax(),
             conditionalPanel(
               condition = "input.showBoth ==true||input.Method==1",
              h4("Isoconversional differential Friedman method [1]"),
              div("$$ \\large \\bbox[white] 
                 {
                 ln\\left(\\frac{d\\alpha}{dt}\\right)_{\\alpha,i} = ln[f(\\alpha) A_\\alpha] - \\frac{E_\\alpha}{RT_{\\alpha,i}}
                 }
                 $$")
              ),
             conditionalPanel(
               condition = "input.showBoth ==true||input.Method==2",
               h4("Advanced integral Vyazovkin method [2]"),
               div("$$ \\large \\bbox[white] 
                 {
                 \\Phi(E_\\alpha) = \\sum_{i=1}^n \\sum_{j \\neq i}^n \\frac{J[E_\\alpha,T_i(t_\\alpha)]}{J[E_\\alpha,T_j(t_\\alpha)]}\\;,
                 }
                 $$
                 "),
               div("$$ \\large \\bbox[white]
                   {
                   \\quad J[E_\\alpha,T(t_\\alpha)] = \\int_{t_{\\alpha-\\Delta\\alpha}}^{t_\\alpha} exp\\left[\\frac{-E_\\alpha}{RT(t)}\\right] dt 
                   }
                   $$")
             ),
             conditionalPanel(
               condition = "input.showKiss ==true",
               h4("Kissinger method [3,4]"),
               div("$$ \\large \\bbox[white] 
                 {
                 ln\\left(\\frac{\\beta_i}{T_{p,i}}\\right) = ln\\left(\\frac{AR}{E_a}\\right) - \\frac{E_a}{RT_{p,i}}
                 }
                 $$")
             ),
             h4("References:"),
             conditionalPanel(
               condition = "input.showBoth ==true||input.Method==1",
               h6("1. Friedman, H. L. Kinetics of Thermal Degradation of Char-Forming Plastics from Thermogravimetry. 
                  Application to a Phenolic Plastic. J. Polym. Sci. Part C Polym. Symp. 1964, 6 (1), 183–195. DOI:10.1002/polc.5070060121")),
             conditionalPanel(
               condition = "input.showBoth ==true||input.Method==2",
               h6("2. Vyazovkin, S. Modification of the Integral Isoconversional Method to Account for Variation in the 
                  Activation Energy. J. Comput. Chem. 2001, 22 (2), 178–183. DOI:10.1002/1096-987X(20010130)22:2<178::AID-JCC5>3.0.CO;2-#")),
             conditionalPanel(
               condition = "input.showKiss ==true",
               h6("3. Kissinger, H. E. Reaction Kinetics in Differential Thermal Analysis. Anal. Chem. 1957, 29 (11), 1702–1706.",br(),
                  "4. ASTM E698-05, Standard Test Method for Arrhenius Kinetic Constants for Thermally Unstable Materials, ASTM International, West Conshohocken, PA, 2005. 10.1520/E0698-11.")),
               
             hr(),
             img(src="thinks.png", height = 75, width = 100),#, align='right'),
             h5(strong("Dr. Nikita V. Muravyev, Version 31.10.21")),
             p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
                  a("Cite",href = "https://doi.org/10.1039/C6CP06498A")," / ",
                  a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
             )
             )
             )
             )
             ))
