test_that("dates_to_posix", {
    data(example_sample_annotation, package = "proBatch")

    sample_annotation <- example_sample_annotation[1:5, ]
    sample_annotation$DateTime <- NULL

    new_annotation <- dates_to_posix(
        sample_annotation,
        time_column = c("RunDate", "RunTime"),
        new_time_column = "DateTime",
        dateTimeFormat = c("%b_%d", "%H:%M:%S"),
        tz = "GMT"
    )

    expect_equal(new_annotation$RunDate[1], "Oct_05")
    expect_equal(new_annotation$RunTime[1], "18:35:00")
    expect_equal(format(new_annotation$DateTime[1], "%d"), "05")
    expect_equal(format(new_annotation$DateTime[1], "%H"), "18")

    expect_s3_class(new_annotation$DateTime, "POSIXct")
})


test_that("date_to_sample_order", {
    data(example_sample_annotation, package = "proBatch")

    sample_annotation_test <- example_sample_annotation[1:5, ]
    sample_annotation_test$DateTime <- NULL
    sample_annotation_test$order <- NULL

    new_annotation_worder <- date_to_sample_order(
        sample_annotation_test,
        time_column = c("RunDate", "RunTime"),
        new_time_column = "new_DateTime",
        dateTimeFormat = c("%b_%d", "%H:%M:%S"),
        new_order_col = "new_order",
        instrument_col = NULL
    )

    expect_equal(new_annotation_worder$RunDate[1], "Oct_05")
    expect_equal(new_annotation_worder$RunTime[1], "18:35:00")

    expect_s3_class(new_annotation_worder$new_DateTime, "POSIXct")
    expect_length(new_annotation_worder$new_order, nrow(sample_annotation_test))

    expect_equal(new_annotation_worder$new_order[1], 1)
    expect_equal(new_annotation_worder$new_order[2], 2)
})


test_that("dates_to_posix: single-column fallback works", {
    data(example_sample_annotation, package = "proBatch")
    df <- example_sample_annotation[1:3, ]

    # keep only RunDate, parse only that
    df$RunTime <- NULL

    out <- dates_to_posix(
        sample_annotation = df,
        time_column = "RunDate",
        # leave new_time_column NULL -> should default to time_column
        new_time_column = NULL,
        dateTimeFormat = "%b_%d",
        tz = "GMT"
    )

    # original RunDate still character
    expect_type(df$RunDate, "character")

    # new column created (same name) and class POSIXct
    expect_s3_class(out$RunDate, "POSIXct")

    # check that day is parsed correctly
    expect_equal(
        format(out$RunDate[1], "%d"),
        format(as.Date(df$RunDate[1], "%b_%d"), "%d")
    )
})

test_that("dates_to_posix: error when lengths mismatch", {
    data(example_sample_annotation, package = "proBatch")
    df <- example_sample_annotation[1:3, ]

    expect_error(
        dates_to_posix(
            sample_annotation = df,
            time_column = c("RunDate", "RunTime"),
            new_time_column = "DT",
            dateTimeFormat = "%b_%d", # only length 1
            tz = "GMT"
        ),
        "`dateTimeFormat` must match length of `time_column`"
    )
})

test_that("dates_to_posix restores LC_TIME after success and error", {
    original_locale <- Sys.getlocale("LC_TIME")
    on.exit(
        suppressWarnings(Sys.setlocale("LC_TIME", original_locale)),
        add = TRUE
    )

    locale_candidates <- unique(c(
        original_locale,
        "C.UTF-8",
        "C.utf8",
        "en_US.UTF-8",
        "en_US.utf8",
        "English_United States.1252"
    ))
    alternate_locale <- NULL
    for (candidate in locale_candidates) {
        changed <- suppressWarnings(Sys.setlocale("LC_TIME", candidate))
        if (!is.na(changed) && Sys.getlocale("LC_TIME") != "C") {
            alternate_locale <- Sys.getlocale("LC_TIME")
            break
        }
    }
    suppressWarnings(Sys.setlocale("LC_TIME", original_locale))
    skip_if(is.null(alternate_locale), "No non-C LC_TIME locale is available")

    suppressWarnings(Sys.setlocale("LC_TIME", alternate_locale))
    locale_before_call <- Sys.getlocale("LC_TIME")

    out <- dates_to_posix(
        sample_annotation = data.frame(date = "2026-01-02"),
        time_column = "date",
        new_time_column = "parsed",
        dateTimeFormat = "%Y-%m-%d",
        locale = "C"
    )
    expect_s3_class(out$parsed, "POSIXct")
    expect_identical(Sys.getlocale("LC_TIME"), locale_before_call)

    expect_error(
        dates_to_posix(
            sample_annotation = data.frame(date = "2026-01-02"),
            time_column = c("date", "missing"),
            new_time_column = "parsed",
            dateTimeFormat = "%Y-%m-%d",
            locale = "C"
        ),
        "`dateTimeFormat` must match length of `time_column`",
        fixed = TRUE
    )
    expect_identical(Sys.getlocale("LC_TIME"), locale_before_call)
})

test_that("date_to_sample_order: grouping by instrument resets ranks", {
    data(example_sample_annotation, package = "proBatch")
    df <- example_sample_annotation[1:6, ]

    # create two instruments repeating
    df$instrument <- rep(c("A", "B"), each = 3)

    out <- date_to_sample_order(
        sample_annotation = df,
        time_column = c("RunDate", "RunTime"),
        new_time_column = "DT",
        dateTimeFormat = c("%b_%d", "%H:%M:%S"),
        new_order_col = "ord",
        instrument_col = "instrument"
    )

    # Within each instrument group, the minimum ord must be 1
    ord_by_inst <- split(out$ord, out$instrument)
    expect_equal(ord_by_inst$A, 1:3)
    expect_equal(ord_by_inst$B, 1:3)
})

test_that("date_to_sample_order resolves tied times by arranged order", {
    tied <- data.frame(
        sample = c("a2", "a1", "b2", "b1"),
        instrument = c("A", "A", "B", "B"),
        run_time = rep("2026-01-02 10:00:00", 4)
    )

    ungrouped <- suppressWarnings(date_to_sample_order(
        sample_annotation = tied,
        time_column = "run_time",
        new_time_column = "parsed",
        dateTimeFormat = "%Y-%m-%d %H:%M:%S",
        new_order_col = "run_order",
        instrument_col = NULL
    ))
    expect_equal(ungrouped$sample, tied$sample)
    expect_equal(ungrouped$run_order, seq_len(nrow(tied)))

    grouped <- suppressWarnings(date_to_sample_order(
        sample_annotation = tied,
        time_column = "run_time",
        new_time_column = "parsed",
        dateTimeFormat = "%Y-%m-%d %H:%M:%S",
        new_order_col = "run_order",
        instrument_col = "instrument"
    ))
    expect_equal(
        unname(split(grouped$run_order, grouped$instrument)),
        list(1:2, 1:2)
    )
})
