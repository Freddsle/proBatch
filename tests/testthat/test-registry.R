test_that("registry preserves legacy Core step registrations", {
    expect_true(pb_has_step("log"))
    expect_true(pb_has_step("log2"))
    expect_true(pb_has_step("medianNorm"))
    expect_identical(pb_list_steps("^median"), "medianNorm")
    expect_true(pb_has_step("loesslimmarbe"))
    expect_identical(
        pb_list_steps("^loesslimmarbe$"),
        "loessLimmaRBE"
    )
    expect_false("loesslimmarbe" %in% pb_list_steps())

    details <- pb_list_steps("^median", details = TRUE)
    expect_s4_class(details, "DataFrame")
    expect_identical(as.character(details$step), "medianNorm")
    expect_identical(as.character(details$name), "medianNorm")
    expect_identical(as.character(details$pkg), "proBatch")
    expect_identical(as.character(details$package), "proBatch")
})

test_that("canonical registration and aliases are collision-safe", {
    provider <- "stats"
    on.exit(pb_unregister_steps(provider), add = TRUE)

    step_fun <- function(m) m
    expect_true(isTRUE(pb_register_step(
        "registry_alpha",
        step_fun,
        package = provider,
        kind = "normalization",
        requires = "utils",
        aliases = c("registry_a", "registry_shared"),
        label = "Registry alpha"
    )))
    before <- pb_list_steps(details = TRUE)

    expect_true(isTRUE(pb_register_step(
        "registry_alpha",
        step_fun,
        package = provider,
        kind = "normalization",
        requires = "utils",
        aliases = c("registry_a", "registry_shared"),
        label = "Registry alpha"
    )))
    expect_identical(pb_list_steps(details = TRUE), before)

    expect_error(
        pb_register_step(
            "registry_alpha",
            function(m) m + 1,
            package = provider,
            kind = "normalization",
            requires = "utils",
            aliases = c("registry_a", "registry_shared"),
            label = "Registry alpha"
        ),
        "already registered with different"
    )
    expect_error(
        pb_register_step(
            "registry_shared",
            step_fun,
            package = "utils"
        ),
        "collides with alias"
    )
    expect_error(
        pb_register_step(
            "registry_beta",
            step_fun,
            package = "utils",
            aliases = "registry_alpha"
        ),
        "collides with canonical"
    )

    expect_false(pb_has_step("registry_beta"))
    expect_identical(
        pb_list_steps("^registry_shared$"),
        "registry_alpha"
    )
    expect_identical(
        pb_list_steps("^Registry alpha$"),
        "registry_alpha"
    )
})

test_that("registry details retain legacy and provider metadata", {
    provider <- "stats"
    on.exit(pb_unregister_steps(provider), add = TRUE)

    step_fun <- function(m, sample_annotation = NULL) m
    pb_register_step(
        "registry_detail",
        step_fun,
        package = provider,
        kind = "correction",
        requires = "utils",
        aliases = c("registry_d1", "registry_d2"),
        label = "Registry detail label"
    )

    details <- pb_list_steps("^registry_d1$", details = TRUE)
    expect_s4_class(details, "DataFrame")
    expect_true(
        all(
            c(
                "step",
                "pkg",
                "env",
                "n_formals",
                "name",
                "package",
                "kind",
                "label",
                "requires",
                "aliases",
                "available"
            ) %in%
                colnames(details)
        )
    )
    expect_identical(as.character(details$step), "registry_detail")
    expect_identical(as.character(details$name), "registry_detail")
    expect_identical(as.character(details$pkg), provider)
    expect_identical(as.character(details$package), provider)
    expect_identical(as.character(details$kind), "correction")
    expect_identical(as.character(details$label), "Registry detail label")
    expect_identical(details$requires[[1L]], "utils")
    expect_identical(
        details$aliases[[1L]],
        c("registry_d1", "registry_d2")
    )
    expect_true(details$available[[1L]])
    expect_identical(details$n_formals[[1L]], 2L)

    expect_identical(
        pb_list_steps("^registry_d1$", available = TRUE),
        "registry_detail"
    )
    expect_identical(
        pb_list_steps("^registry_d1$", available = FALSE),
        character()
    )
})

test_that("aliases resolve canonically with provider identity", {
    provider <- "stats"
    on.exit(pb_unregister_steps(provider), add = TRUE)

    step_fun <- function(m) m
    pb_register_step(
        "registry_resolve",
        step_fun,
        package = provider,
        aliases = "registry_alias"
    )

    expect_true(pb_has_step("registry_alias"))
    expect_true(pb_has_step("registry_alias", available = TRUE))

    resolved <- .pb_resolve_step("registry_alias")
    expect_identical(resolved$input_name, "registry_alias")
    expect_identical(resolved$name, "registry_resolve")
    expect_identical(resolved$package, provider)
    expect_identical(resolved$fun, step_fun)
    expect_true(resolved$available)
    expect_true(resolved$registered)

    expect_error(
        .pb_resolve_step(
            "registry_alias",
            package = "utils",
            require_available = FALSE
        ),
        "not recorded provider 'utils'"
    )
    expect_error(
        .pb_resolve_step(
            "registry_alias",
            package = NA_character_,
            require_available = FALSE
        ),
        "not an ordinary registration"
    )
})

test_that("resolver accepts direct functions without registration", {
    step_fun <- function(m) m

    resolved <- .pb_resolve_step(step_fun)
    expect_identical(resolved$input_name, step_fun)
    expect_identical(resolved$fun, step_fun)
    expect_identical(resolved$name, NA_character_)
    expect_true(resolved$available)
    expect_false(resolved$registered)

    expect_identical(
        .pb_resolve_step(step_fun, package = NA_character_)$fun,
        step_fun
    )
})

test_that("availability reports unloaded providers and requirements", {
    unavailable_provider <- "proBatchRegistryProviderThatDoesNotExist"
    missing_requirement <- "proBatchRegistryRequirementThatDoesNotExist"
    on.exit(pb_unregister_steps(unavailable_provider), add = TRUE)
    on.exit(pb_unregister_steps("stats"), add = TRUE)

    never_run <- function(m) stop("provider function must not run")
    expect_silent(pb_register_step(
        "registry_unloaded_provider",
        never_run,
        package = unavailable_provider
    ))
    expect_true(pb_has_step("registry_unloaded_provider"))
    expect_false(pb_has_step(
        "registry_unloaded_provider",
        available = TRUE
    ))
    expect_identical(
        pb_list_steps(
            "^registry_unloaded_provider$",
            available = FALSE
        ),
        "registry_unloaded_provider"
    )
    expect_error(
        .pb_resolve_step("registry_unloaded_provider"),
        "provider namespace .* is not loaded"
    )

    expect_silent(pb_register_step(
        "registry_missing_requirement",
        never_run,
        package = "stats",
        requires = missing_requirement
    ))
    expect_false(pb_has_step(
        "registry_missing_requirement",
        available = TRUE
    ))
    expect_error(
        .pb_resolve_step("registry_missing_requirement"),
        "requires unavailable package"
    )

    expect_error(
        .pb_resolve_step(
            "registry_provider_replay",
            package = unavailable_provider
        ),
        paste0(
            "recorded from provider '",
            unavailable_provider,
            "' is not registered"
        ),
        fixed = TRUE
    )
    expect_error(
        .pb_resolve_step(
            "registry_ordinary_replay",
            package = NA_character_
        ),
        "recorded as an ordinary registration is not registered",
        fixed = TRUE
    )
})

test_that("availability inspection does not load requirement namespaces", {
    candidate <- setdiff(
        rownames(utils::installed.packages()),
        loadedNamespaces()
    )
    skip_if(!length(candidate), "No installed unloaded namespace is available")
    candidate <- candidate[[1L]]
    provider <- "stats"
    on.exit(pb_unregister_steps(provider), add = TRUE)

    pb_register_step(
        "registry_nonloading_availability",
        identity,
        package = provider,
        requires = candidate
    )
    expect_false(candidate %in% loadedNamespaces())
    expect_true(pb_has_step(
        "registry_nonloading_availability",
        available = TRUE
    ))
    expect_identical(
        pb_list_steps(
            "^registry_nonloading_availability$",
            available = TRUE
        ),
        "registry_nonloading_availability"
    )
    expect_false(candidate %in% loadedNamespaces())
})

test_that("provider unregister removes canonical names and aliases", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)

    pb_register_step(
        "registry_remove_two",
        identity,
        package = provider
    )
    pb_register_step(
        "registry_remove_one",
        identity,
        package = provider,
        aliases = "registry_remove_alias"
    )

    removed <- withVisible(pb_unregister_steps(provider))
    expect_identical(
        removed$value,
        c("registry_remove_one", "registry_remove_two")
    )
    expect_false(removed$visible)
    expect_false(pb_has_step("registry_remove_one"))
    expect_false(pb_has_step("registry_remove_alias"))
    expect_false(pb_has_step("registry_remove_two"))
    expect_identical(
        withVisible(pb_unregister_steps(provider)),
        list(value = character(), visible = FALSE)
    )
})

test_that("selective unregister removes only selected steps and aliases", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)

    for (name in c("registry_select_a", "registry_select_b", "registry_keep")) {
        pb_register_step(
            name,
            identity,
            package = provider,
            aliases = paste0(name, "_alias")
        )
    }

    removed <- withVisible(pb_unregister_steps(
        provider,
        steps = c("registry_select_b", "registry_select_a")
    ))
    expect_identical(removed$value, c("registry_select_a", "registry_select_b"))
    expect_false(removed$visible)
    expect_false(pb_has_step("registry_select_a"))
    expect_false(pb_has_step("registry_select_a_alias"))
    expect_false(pb_has_step("registry_select_b"))
    expect_false(pb_has_step("registry_select_b_alias"))
    expect_true(pb_has_step("registry_keep"))
    expect_true(pb_has_step("registry_keep_alias"))
})

test_that("an empty unregister selection leaves the registry unchanged", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)
    pb_register_step(
        "registry_empty_selection",
        identity,
        package = provider,
        aliases = "registry_empty_selection_alias"
    )
    before <- pb_list_steps(details = TRUE)

    expect_identical(
        withVisible(pb_unregister_steps(provider, steps = character())),
        list(value = character(), visible = FALSE)
    )
    expect_identical(pb_list_steps(details = TRUE), before)
})

test_that("selective unregister ignores absent names and is idempotent", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)
    pb_register_step("registry_retry", identity, package = provider)

    expect_identical(
        pb_unregister_steps(
            provider,
            steps = c("registry_absent", "registry_retry")
        ),
        "registry_retry"
    )
    expect_identical(
        withVisible(pb_unregister_steps(provider, steps = "registry_retry")),
        list(value = character(), visible = FALSE)
    )
    expect_identical(
        withVisible(pb_unregister_steps(provider, steps = "registry_absent")),
        list(value = character(), visible = FALSE)
    )
})

test_that("selective unregister rejects foreign ownership before mutation", {
    provider <- "proBatchExampleProvider"
    other_provider <- "proBatchOtherExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)
    on.exit(pb_unregister_steps(other_provider), add = TRUE)
    pb_register_step(
        "registry_owned",
        identity,
        package = provider,
        aliases = "registry_owned_alias"
    )
    pb_register_step(
        "registry_foreign",
        identity,
        package = other_provider,
        aliases = "registry_foreign_alias"
    )
    before <- pb_list_steps(details = TRUE)

    expect_error(
        pb_unregister_steps(
            provider,
            steps = c("registry_owned", "registry_foreign")
        ),
        "not owned by provider"
    )
    expect_identical(pb_list_steps(details = TRUE), before)
})

test_that("selective unregister rejects aliases before mutation", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)
    pb_register_step(
        "registry_canonical_only",
        identity,
        package = provider,
        aliases = "registry_not_canonical"
    )
    before <- pb_list_steps(details = TRUE)

    expect_error(
        pb_unregister_steps(
            provider,
            steps = c("registry_canonical_only", "registry_not_canonical")
        ),
        "canonical steps, not aliases"
    )
    expect_identical(pb_list_steps(details = TRUE), before)
})

test_that("selective unregister rejects malformed selections without mutation", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)
    pb_register_step(
        "registry_invalid_selection",
        identity,
        package = provider,
        aliases = "registry_invalid_selection_alias"
    )
    before <- pb_list_steps(details = TRUE)
    invalid_selections <- list(
        c("registry_invalid_selection", "registry_invalid_selection"),
        c("registry_invalid_selection", NA_character_),
        1,
        c("registry_invalid_selection", "")
    )

    for (steps in invalid_selections) {
        expect_error(pb_unregister_steps(provider, steps = steps), "`steps`")
        expect_identical(pb_list_steps(details = TRUE), before)
    }
})

test_that("selective refresh preserves collisions for retained registrations", {
    provider <- "proBatchExampleProvider"
    on.exit(pb_unregister_steps(provider), add = TRUE)
    for (name in c("registry_refresh", "registry_retained")) {
        pb_register_step(
            name,
            identity,
            package = provider,
            label = "Original label"
        )
    }
    pb_unregister_steps(provider, steps = "registry_refresh")

    expect_true(isTRUE(pb_register_step(
        "registry_refresh",
        identity,
        package = provider,
        label = "Updated label"
    )))
    details <- pb_list_steps("^registry_refresh$", details = TRUE)
    expect_identical(as.character(details$label), "Updated label")
    before <- pb_list_steps(details = TRUE)

    expect_error(
        pb_register_step(
            "registry_retained",
            identity,
            package = provider,
            label = "Updated label"
        ),
        "already registered with different"
    )
    expect_identical(pb_list_steps(details = TRUE), before)
})

test_that("two-argument registration remains supported", {
    on.exit(pb_unregister_steps("base"), add = TRUE)

    expect_true(isTRUE(pb_register_step("registry_legacy", identity)))
    expect_true(pb_has_step("registry_legacy"))
    expect_true(pb_has_step("registry_legacy", available = TRUE))

    details <- pb_list_steps("^registry_legacy$", details = TRUE)
    expect_identical(as.character(details$package), "base")
})

test_that("registry rejects malformed metadata without mutation", {
    before <- pb_list_steps()

    expect_error(pb_register_step("", identity), "`name`")
    expect_error(pb_register_step("registry_bad", 1), "`fun`")
    expect_error(
        pb_register_step(
            "registry_bad",
            identity,
            aliases = "registry_bad"
        ),
        "must not contain"
    )
    expect_error(
        pb_register_step(
            "registry_bad",
            identity,
            aliases = c("duplicate", "duplicate")
        ),
        "must not contain duplicates"
    )
    expect_error(
        pb_register_step(
            "registry_bad",
            identity,
            requires = NA_character_
        ),
        "`requires`"
    )

    expect_identical(pb_list_steps(), before)
})
