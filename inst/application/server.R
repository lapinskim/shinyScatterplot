shinyServer(function(input, output) {
    output$scatterplot <- renderPlot({
        shinyScatterplot::plotLogMat(
            data[[input$data]],
            title = input$data,
            xlim = input$xlim,
            ylim = input$ylim,
            type = input$type,
            bins = input$bins,
            colour = c(input$CPA, input$NoCPA)
        )
    },
    width = "auto",
    height = "auto")

    output$desc <- renderText({
        "Relationship between ∆TS and ∆TL for maternal non-CPA and CPA genes, \
        expressed as a log2 fold change value. Boxplots that represent ∆TS and \
        ∆TL are depicted for the X- and Y-axes, respectively."
    })

    })
