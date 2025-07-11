test_that("geom_split_violin returns custom layer", {
    layer <- geom_split_violin()
    expect_s3_class(layer, "LayerInstance")
    expect_true(inherits(layer$geom, "GeomSplitViolin"))
})


test_that("split violin plot builds expected layers", {
    set.seed(123)
    df <- data.frame(
        x = rep(c("A", "B"), each = 100),
        m = rep(c("C", "D"), 100),
        y = rnorm(200)
    )
    p <- plot_split_violin_with_boxplot(
        df,
        y_col = "y",
        col_for_color = "m",
        col_for_box = "x",
        colors_for_plot = c("orange", "purple"),
        hlineintercept = 5,
        plot_title = "Title"
    )

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::get_expr(p$mapping$x), rlang::sym("x"))
    expect_equal(rlang::get_expr(p$mapping$fill), rlang::sym("m"))
    expect_equal(p$labels$title, "Title")

    expect_true(any(sapply(p$layers, function(x) inherits(x$geom, "GeomSplitViolin"))))
    expect_true(any(sapply(p$layers, function(x) inherits(x$geom, "GeomBoxplot"))))
    expect_true(any(sapply(p$layers, function(x) inherits(x$geom, "GeomHline"))))

    hline_ind <- which(sapply(p$layers, function(x) inherits(x$geom, "GeomHline")))
    expect_equal(p$layers[[hline_ind]]$data$yintercept, 5)
    expect_equal(p$scales$scales[[1]]$palette(2), c("orange", "purple"))
})


test_that("no hline when argument NULL", {
    df <- data.frame(x = "A", m = "C", y = 1)
    p <- plot_split_violin_with_boxplot(df, hlineintercept = NULL)
    expect_false(any(sapply(p$layers, function(x) inherits(x$geom, "GeomHline"))))
})
