#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine is used to import the experimental data into the kinetic project.
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
    windowTitle = "THINKS/Import",
    title = tags$head(tags$link(rel="icon", 
                                href="favicon.ico", 
                                type="image/x-icon")
    )),
  h2("Creating the kinetic project", style="margin-bottom: 0em"),
  br(),
  tags$head(
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
                    hr { 
                          margin-top: 0em !important;
                          margin-bottom: 1em !important;
                          border-width: 1px;
                        } "))),
  sidebarLayout(
    sidebarPanel(width=3,
        tabsetPanel(
            tabPanel("File info",
                 br(),
                 div(style="display:inline-block", helpText("Signal type: ")),
                 div(style="display:inline-block; margin-left: 1.4em", radioButtons("signal", label=NULL, 
                            c("differential"="dif","integral"="int"), selected = "int", inline=T)),
                 br(),
                 div(style="display:inline-block", helpText("Measurement: ")),
                 div(style="display:inline-block; margin-left: 0.2em", radioButtons("tprg", label=NULL, 
                                                                c("nonisothermal"="niso","isothermal"="iso"), selected = "niso", inline=T)),
                 fluidRow(column(3, helpText("Separator:"),
                                 textInput("sep1",label=NULL,value = ";")
                                ),
                          column(9,
                            fluidRow(
                              column(6, helpText("Drop X first lines: "),
                                 numericInput("SkipL",label=NULL,value = 35, min = 0, max = 100, step = 1)
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
                fluidRow(column(4, uiOutput("thirdCol")),
                         column(4, strong(uiOutput("thirdName"))),
                         column(4, uiOutput("thirdUnits"))
                ),
                                
                 conditionalPanel(
                    condition = "input.tgaUnits == 2 || input.dscUnits == 2",  
                    fluidRow(column(8, helpText("Sample mass, mg:")),
                             column(4, numericInput("mass1",label=NULL,value = 5,min = 0.05, max = 500, step = 0.1, width='100%')))
                    ),
                 fileInput('file1', label="Load .txt file (limited to 5 MB)", accept=c('text/plain')),
                 actionButton("go1", "Add file to project", class="btn btn-primary"),
                 conditionalPanel(
                   condition = "input.go1 != 0",
                   br(),
                   uiOutput("prjData"),
                   downloadButton('downloadData', 'Export kinetic data', class="btn btn-success"))               
            ),
            tabPanel("Advanced",
                 br(),    
                 conditionalPanel(
                       condition = "input.go1 == 0",
                       fluidRow(column(8, helpText("Number of points to be applied for all experiments in this project:*")),
                                column(4, numericInput("Npts",label=NULL,value = 1000,min = 200, max = 2000, step = 50, width='100%')))
                       ),
                 fluidRow(column(8, helpText("Number of Simpson subintervals used for integration:*")),
                          column(4, numericInput("N_Simp",label=NULL,value = 1000,min = 2, max = 2000, step = 1, width='100%'))),
                 
                 div(style="margin-top:0em, margin-bottom:0em", helpText("*limited in online version to 2000 pts")),
                 div(style="display:inline-block", checkboxInput("DSCLess0", "Suppress DSC < 0", value=F)),
                 div(style="display:inline-block", checkboxInput("aLess0", "Allow α < 0", value=F)),
                 div(style="display:inline-block; margin-left: 1em", checkboxInput("aHigher1", "Allow α > 1", value=F)),
                 fluidRow(column(8, helpText("Final conversion degree value:")),
                          column(4, numericInput("afin",label=NULL,value = 1,min = 0.05, max = 1, step = 0.1, width='100%'))),

                  hr(),
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
                                  numericInput("TitleSize", label=NULL, 1.2, min = 0.8, max = 2.5, step=0.1, width='100%')),
                           column(4, helpText("Axis"),
                                  numericInput("AxisSize", label=NULL, 1.7, min = 0.8, max = 2.5, step=0.1, width='100%')),
                           column(4, 
                                  helpText("Legend"),
                                  numericInput("LegendSize", label=NULL, 1.1, min = 0.5, max = 2.5, step=0.1, width='100%'))
                  ),
                  helpText("For plot with the loaded in project data"),
                  div(style="display:inline-block", helpText("Y-axis: ")),
                  div(style="display:inline-block", radioButtons("showA", label=NULL, c("Conversion"=1,"Conversion rate" =0),selected =1,inline=T)),
                  br(),
                  div(style="display:inline-block", helpText("X-axis: ")),
                  div(style="display:inline-block", radioButtons("showT", label=NULL, c("Temperature"=1,"Time" =2),selected =1,inline=T))
            )
        )
                 
    ),
    
    mainPanel(
      #shinythemes::themeSelector(),
      column(width=7,
        tabsetPanel(id = "MainTabset", selected = "panel3",
          tabPanel("This file", value = "panel1",
             conditionalPanel(
               condition = "output.cond == false",
                plotOutput("Plot1", height = "auto"),
                plotOutput("Plot2", height = "auto"))
          ),
          tabPanel("Kinetic project", value = "panel2",
             plotOutput("Plot3", height = "auto"),
             plotOutput("Plot4", height = "auto")
             ),
          tabPanel("Help", value="panel3",
                     div(class = "initialParent", 
                         div( class="initialChild",
                              p(strong("/ Instructions /")),
                              p("1. Prepare the thermal analysis results. The input .txt file should contain 
                        the following quantities - temperature, time and conversion degree (i.e., integral signal) 
                        or conversion rate data (i.e., differential signal), or the experimental signals that give conversion data.
                        Define the file details (separators, column numbers etc.) in the left block and upload the file."),
                              p("2. Correct the baseline for DSC data, define the mass loss range for TGA data using the right block.
                        Note, that the number of points to which the data will be thin-out (shown in 
                        \"Advanced\" tab) is set globally, i.e., after applying to the first file it will be used for others in the project.
                        Usually, after the exclusion of the unimportant regions far from the process (using the \"Drop X lines\" parameters)
                        the interested region is adequately represented by 1000-2000 points. Control visually the data after all manipulations
                        over the original data in the first plot."),
                              p("3. Go to \"Kinetic project\" tab  to see the conversion degree data prepared for loading to kinetic project. 
                        If necessary, smooth the data (again controlling the adequacy over the original data on the first plot). 
                          Press \"Add file to project\" button after all treatment and upload the data for the next experiment."),
                              p("4.  Once all files are in the project (more than four experiments is recommended by ",
                                tags$a(href = "https://doi.org/10.1016/j.tca.2011.03.034", target="_blank", "ICTAC Kinetic committee!"),
                                ") download the resulting file (contains the runs that are selected) pressing the \"Export kinetic data\". Enjoy your coffee and move to the ",
                                tags$a(href = "https://www.thinks.chemphys.ru", target="_blank", "kinetic analysis."))
                         ))
                   )
          )),
      column(width=4,
             conditionalPanel(
               condition = "input.signal == 'dif'", 
               fluidRow(column(6, helpText("Select baseline type:")),
                        column(6, style='padding:0px;',
                        selectInput("BLtype", label =NULL, choices = c("No correction"="no",
                                                                                  "Linear" = "line","Linear-2" = "linst",
                                                                                  "Horizontal-right"="horRight", "Tangential"="tang",
                                                                                  "Spline" = "spln2", #"Spline" = "spln", 
                                                                                  "Spline-2-center" = "spln2c", "Glass transition" = "glasstr"),
                                    selected ="no", width='100%'))
                        )),
             uiOutput("BLpointsBeg"),
             conditionalPanel(
               condition = "input.BLtype == 'spln2c'", 
               uiOutput("BLpointsCntr"),
               fluidRow(column(8, helpText("Weight of the central magnet:")),
                        column(4, numericInput("FNcentrPts", label=NULL, value = 0.2, min = 0, max =20, step = 0.001, width='100%'))
               )),
             uiOutput("BLpointsEnd"),
             uiOutput("Yscale"),
             conditionalPanel( condition = "output.cond == false",
                radioButtons("smType", "Smoothing of raw data:",c("No"="no","Cubic spline" = "cub",
                                                               "Polynomial" = "pol"),selected ="no",inline=T),hr()),  
         #    conditionalPanel( condition = "output.cond2 == false",
        #                       div(style="display:inline-block", "no enough data for smoothing"), br()
        #     ),
             conditionalPanel(
               condition = "input.smType == 'cub'", 
               sliderInput("CubSpar", "Smoothing parameter :",min = 0, max = 1, value = 0.2, step=0.01)),
             conditionalPanel(
               condition = "input.smType == 'pol'", 
               sliderInput("PolSpan", "Smoothing parameter (log):",min = -3, max = 0, value = -2, step=0.01)),
             #hr(),
             img(src="thinks.png", height = 75, width = 100),#, align='right'),
             h5(strong("Dr. Nikita V. Muravyev, Version 08.07.22")),
             p(h6(a("Download codes or feedback", href = "https://www.researchgate.net/project/THINKS")," / ",
                  a("Cite",href = "https://doi.org/10.3390/molecules24122298")," / ",
                  a("Mailto",href = "mailto:n.v.muravyev@chemphys.ru"))
             )
             )
             
      )
             )
             ))
