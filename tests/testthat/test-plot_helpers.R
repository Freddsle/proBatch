test_that(".pb_split_arg_by_assay keeps atomic vectors intact", {
    assays <- paste0("assay", 1:4)
    vec <- c("tech1", "tech2", "tech3")
    res <- proBatch:::`.pb_split_arg_by_assay`(vec, assays)

    expect_equal(length(res), length(assays))
    expect_true(all(vapply(res, identical, logical(1), vec)))
})

test_that(".pb_split_arg_by_assay honours named vectors", {
    assays <- c("assay1", "assay2", "assay3")
    vec <- c(assay1 = "x", assay3 = "z")
    res <- proBatch:::`.pb_split_arg_by_assay`(vec, assays)

    expect_equal(res[[1]], "x")
    expect_equal(res[[2]], "z")
    expect_equal(res[[3]], "z")
})
