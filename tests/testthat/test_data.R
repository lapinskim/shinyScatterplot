context("Dataset processing")
library(shinyScatterplot)

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

# sample only 100 genes for testing
ngenes <- max(vapply(
    list(bound_ncounts,
         total_ncounts,
         totalFpkm,
         unbound_ncounts),
    nrow,
    numeric(1)
))
set.seed(1234)
gsample <- sample(ngenes, 100)

# Create a test dataset
stage_test <- "1Cell"
total_test <- total_ncounts[gsample,]
bound_test <- bound_ncounts[gsample,]
unbound_test <- unbound_ncounts[gsample,]
FPKM_test <- totalFpkm[gsample,]

test_that("getRTvsRB_FPKM function produces correct values", {
    expect_equal_to_reference(
        getRTvsRB_FPKM(
            stage = stage_test,
            totalCounts = total_test,
            boundCounts = bound_test,
            unboundCounts = unbound_test,
            totalFPKM = FPKM_test,
            cluster = clustering,
            cutoff = 1
        ),
        file = "../cache/data_test.cache"
    )
})


test_that("getRTvsRB_FPKM function merges the clustring data properly", {
    expect_output(
        getRTvsRB_FPKM(
            stage = stage_test,
            totalCounts = total_test,
            boundCounts = bound_test,
            unboundCounts = unbound_test,
            totalFPKM = FPKM_test,
            cluster = clustering,
            cutoff = 1
        ),
        regexp = "unclustered genes: 16"
    )
})
