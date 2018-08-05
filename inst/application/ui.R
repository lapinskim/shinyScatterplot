shinyUI(fluidPage(
    tags$div(
        tags$h2("Polysome association dynamics"),
        tags$body(
            tags$p(
                "Figure 2. published in",
                tags$a(
                    href = 'http://dev.biologists.org/content/early/2017/12/08/dev.159566',
                    'Cytoplasmic polyadenylation-mediated translational control \
                    of maternal mRNAs directs maternal to zygotic transition'
                ),
                "Winata CL, Łapiński M et al., Development, 11 December 2017"
                ),
            tags$br(),
            tags$p(
                "Measurements of transcription (\u2206TS) and translation \
                (\u2206TL) rates were obtained as a fraction of total and \
                polysome-bound expression values at particular stages compared \
                with the egg stage as baseline. Translational regulation is \
                defined as a non-linear relationship between transcription and \
                translation rates."
            ),
            tags$br(),
            tags$div(
                align = "center",
                tags$em("transcription rate (\u2206TS) = total(stage x) / total(egg)"),
                tags$br(),
                tags$em("translation rate(\u2206TL) = polysome(stage x) / polysome(egg)"),
                tags$br(),
                tags$br()
            )

            )
        ),

    sidebarLayout(
        sidebarPanel(
            selectInput("data", "Developmental stage:",
                        names(data),
                        selected = "1Cell"),
            colourpicker::colourInput(
                "CPA",
                "Maternal CPA colour:",
                value = "#FDBF6F",
                palette = "limited",
                allowedCols = RColorBrewer::brewer.pal(12,
                                                       "Paired")
            ),
            colourpicker::colourInput(
                "NoCPA",
                "Maternal Non-CPA colour:",
                value = "#1F78B4",
                palette = "limited",
                allowedCols = RColorBrewer::brewer.pal(12,
                                                       "Paired")
            ),
            sliderInput(
                "xlim",
                "X axis limits:",
                min = -10,
                max = 10,
                value = c(-5, 5)
            ),
            sliderInput(
                "ylim",
                "Y axis limits:",
                min = -10,
                max = 10,
                value = c(-5, 5)
            ),
            selectInput(
                "type",
                "Density plot type:",
                c(Histogram = "histogram",
                  Boxplot = "boxplot"),
                selected = "boxplot"
            ),
            conditionalPanel(
                condition = "input.type == 'histogram'",
                sliderInput(
                    "bins",
                    "Number of bins:",
                    min = 10,
                    max = 200,
                    step = 10,
                    value = 100
                )
            )
        ),

        mainPanel(
            plotOutput(outputId = "scatterplot"),
            textOutput(outputId = "desc")
        )
    )
    ))
