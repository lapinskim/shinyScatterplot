library(shiny)


# load the external data
load(
    system.file(
        "extdata",
        "bound_ncounts.rda" ,
        package = "shinyScatterplot",
        mustWork = TRUE
    )
)
load(
    system.file(
        "extdata",
        "unbound_ncounts.rda" ,
        package = "shinyScatterplot",
        mustWork = TRUE
    )
)
load(
    system.file(
        "extdata",
        "total_ncounts.rda" ,
        package = "shinyScatterplot",
        mustWork = TRUE
    )
)
load(
    system.file(
        "extdata",
        "totalFpkm.rda" ,
        package = "shinyScatterplot",
        mustWork = TRUE
    )
)
load(
    system.file(
        "extdata",
        "clustering.rda" ,
        package = "shinyScatterplot",
        mustWork = TRUE
    )
)

# sample only 5000 genes for faster processing of the example
ngenes <- max(vapply(
    list(bound_ncounts,
         total_ncounts,
         totalFpkm,
         unbound_ncounts),
    nrow,
    numeric(1)
))
set.seed(1234)
gsample <- sample(ngenes, 5000)

stages <- c("1Cell", "16Cell", "128Cell", "3.5hpf", "5.3hpf")

data <- lapply(
    stages,
    getRTvsRB_FPKM,
    total_ncounts[gsample, ],
    bound_ncounts[gsample, ],
    unbound_ncounts[gsample, ],
    totalFpkm[gsample, ],
    clustering
)

names(data) <- stages
