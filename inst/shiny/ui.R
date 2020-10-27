

## Main UI page for Microbiome Explorer
shinyUI(
    dashboardPage(
        title = "Microbiome Explorer",
        ## create the three bars icon square (STACK OVERFLOW)
        dashboardHeader(
            titleWidth = 72,
            title = a(
                href = "#", class = "sidebar-toggle",
                `data-toggle` = "offcanvas", role = "button",
                span(
                    class = "sr-only",
                    "Toggle navigation"
                ),
                icon("bars")
            )
        ),
        ## main sidebar menu items
        dashboardSidebar(
            width = 220, collapsed = TRUE,
            sidebarMenu(
                id = "sidebar_menu",
                menuItem("Data Input",
                    icon = icon("folder"),
                    tabName = "datamenu"
                ),
                menuItem("  Analysis",
                    icon = icon("dna"),
                    tabName = "analysismenu"
                ),
                menuItem("    Report",
                    icon = icon("file"),
                    tabName = "reportmenu"
                ),
                menuItem("User Guide",
                         icon = icon("far fa-hand-point-right"),
                         tabName = "howtomenu"
                ),
                menuItem("     About",
                         icon = icon("fas fa-info-circle"),
                         tabName = "aboutmenu"
                )
            )
        ),
        dashboardBody(
            tags$head(
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "me.css"
                )
            ),
            ## capture interaction on phenorow remove buttons (data input page)
            tags$head(tags$script(HTML(
            "$(document).on('click', '.phenoremoverow .action-button', 
                function () {
                    Shiny.setInputValue('loadnfilter-last_btn', 
                        this.id,{priority: 'event'});
                });"))),
            tags$script(HTML(
            "$('.ellipsis').tooltip({show: {effect:'none', delay:0}});")),


            ## Place the header to the right of the menu button.
            ## From STACK OVERFLOW
            tags$script(HTML('
                    $(document).ready(function() {
                      $("header").find("nav").append(\'<div  class="meIcon"><img src="dishart.png" height = "50"></img><span class="meTitle">Microbiome Explorer</span></div>\');
                    })
                  ')),
            shinyjs::useShinyjs(),
            shinyjs::extendShinyjs(script = "JS/parsePhenoFilters.js",
                                   functions = c("parsePhenoFilters")),
            shinyjs::extendShinyjs(script = "JS/enableTabs.js",
                                   functions = c("disableTab",
                                                 "enableTab",
                                                 "enableAll",
                                                 "collapse")),
            shinyjs::extendShinyjs(script = "JS/traceName.js",
                                   functions = c("getTraceName")),
            shinyjs::extendShinyjs(script = "JS/handleModebar.js",
                                   functions = c("resetAxes")),
            shinyjs::extendShinyjs(script = "JS/changeTable.js",
                                   functions = c("moveButton",
                                                 "removeInputElems")),

            tabItems(
                ##  DATA INPUT PAGES 
                tabItem(
                    tabName = "datamenu",
                    tabBox(
                        id = "data_input_tabbox", 
                        title = "Data Input", width = 12, side = "left",
                        tabPanel(
                            "LOAD & FILTER",
                            dataInputUI("loadnfilter")
                        ),
                        tabPanel(
                            "PHENOTYPE",
                            phenotypeTableUI("phenotab")
                        ),
                        tabPanel(
                            "FEATURES",
                            featureTableUI("featuretab")
                        )
                    )
                ),
                ## ANALYSIS PAGES 
                tabItem(
                    tabName = "analysismenu",
                    column(width = 12,
                           aggregationTabUI("aggregateTab"),
                           tabBox(
                               id = "analysismenu_tabbox", 
                               title = "Analysis", width = 12, side = "left",
                               tabPanel(
                                   "INTRA SAMPLE",
                                   intraAnalysisUI("intraAnalysis")
                               ),
                               tabPanel(
                                   "INTRA FEATURE",
                                   featureAnalysisUI("featureAnalysis")
                               ),
                               tabPanel(
                                   "INTER SAMPLE",
                                   interAnalysisUI("interAnalysis")
                               ),
                               tabPanel(
                                   "CORRELATION",
                                   corrAnalysisUI("corrAnalysis")
                               ),
                               tabPanel(
                                   "DIFFERENTIAL",
                                   diffAnalysisUI("diffAnalysis")
                               ),
                               tabPanel(
                                   "LONGITUDINAL",
                                   longAnalysisUI("longAnalysis")
                               )
                           )
                    )
                ),
                ## REPORT PAGE  
                tabItem(
                    tabName = "reportmenu",
                    tabBox(
                        id = "report_tabbox", title = "Report generation", 
                        width = 12, side = "left",
                        tabPanel(
                            "GENERATE REPORT",
                            reportListUI("reportlist")
                        )
                    )
                ),
                tabItem(
                    tabName = "howtomenu",
                    tabBox(
                        id = "howto_tabbox", title = "User Guide", 
                        width = 12, side = "left",
                        tabPanel(
                            "HOW TO ...",
                            shiny::includeMarkdown("Tutorial.md")
                        )
                    )
                ),
                tabItem(
                    tabName = "aboutmenu",
                    tabBox(
                        id = "about_tabbox", title = "About", 
                        width = 12, side = "left",
                        tabPanel(
                            "PEOPLE",
                            h4("Janina Reeder (Genentech - gRED OMNI Bioinformatics)"),
                            p("Janina developed the Microbiome Explorer under the mentorship of Joseph Paulson from a prototype version by Mo Huang."),
                            h4("Mo Huang (PhD student - Wharton University)"),
                            p("Mo implemented an initial prototype of the application as an intern project under the supervision of Joseph."),
                            h4("Joseph N. Paulson (Genentech - PDBB)"),
                            p("Joseph wrote the metagenomeSeq R package upon which many of the statistical analyses in this application are based.")
                        )
                    )
                )
            )
        )
    )
)
