library(shiny)
library(bslib)
library(DT)
library(gt)
library(ggplot2)
library(openxlsx)
library(xgboost)
library(dplyr)

# ── Locate the data directory ─────────────────────────────────────────────────
# Handles running from project root or from scripts/ via runApp()
if (dir.exists("app_data_files")) {
  data_dir <- "app_data_files"
} else if (dir.exists("../app_data_files")) {
  data_dir <- "../app_data_files"
} else {
  stop("Cannot find app_data_files/ directory. Run the app from the project root.")
}

# Make app_data_files/ accessible to the browser for serving images
addResourcePath("imgs", data_dir)

# ── Load model results and training data ─────────────────────────────────────
final_models <- readRDS(file.path(data_dir, "xgb_combos.RDS"))
training_data <- readRDS(file.path(data_dir, "model_data.RDS"))

# Extract XGBoost model objects (row 1 = mets train, row 3 = mets_red train)
mets_xgb_model <- final_models$final_model[[1]]$model
mets_red_xgb_model <- final_models$final_model[[3]]$model

# ── Model configuration ──────────────────────────────────────────────────────
# Each model type has its own variables, metrics, and training sample counts
model_config <- list(
  mets = list(
    label = "All Papaya (Red + Yellow)",
    model = mets_xgb_model,
    variables = strsplit(
      final_models$model_vars[final_models$dataset == "mets"][1], ", "
    )[[1]],
    r_squared = final_models |>
      filter(dataset == "mets", dataset_type == "test") |> pull(r2),
    rmse = final_models |>
      filter(dataset == "mets", dataset_type == "test") |> pull(rmse),
    n_train = final_models |>
      filter(dataset == "mets", dataset_type == "train") |> pull(n),
    n_test = final_models |>
      filter(dataset == "mets", dataset_type == "test") |> pull(n),
    importance = final_models$final_model[[1]]$importance
  ),
  mets_red = list(
    label = "Red-Fleshed Only",
    model = mets_red_xgb_model,
    variables = strsplit(
      final_models$model_vars[final_models$dataset == "mets_red"][1], ", "
    )[[1]],
    r_squared = final_models |>
      filter(dataset == "mets_red", dataset_type == "test") |> pull(r2),
    rmse = final_models |>
      filter(dataset == "mets_red", dataset_type == "test") |> pull(rmse),
    n_train = final_models |>
      filter(dataset == "mets_red", dataset_type == "train") |> pull(n),
    n_test = final_models |>
      filter(dataset == "mets_red", dataset_type == "test") |> pull(n),
    importance = final_models$final_model[[3]]$importance
  )
)

# ── Training data ranges ─────────────────────────────────────────────────────
# Compute min/max/median for each metabolite from the actual training data
compute_ranges <- function(data, variables) {
  lapply(variables, function(var) {
    vals <- data[[var]]
    data.frame(
      Metabolite = var,
      Min = min(vals, na.rm = TRUE),
      Max = max(vals, na.rm = TRUE),
      Median = median(vals, na.rm = TRUE)
    )
  }) |>
    bind_rows() |>
    mutate(
      Units = case_when(
        Metabolite == "Brix" ~ "\u00b0Brix",
        Metabolite == "TA" ~ "% citric acid",
        Metabolite %in% c("Sucrose", "total_sugar") ~ "g/100g FW",
        .default = "ppb (\u00b5g/kg FW)"
      )
    )
}

training_ranges <- list(
  mets = compute_ranges(training_data, model_config$mets$variables),
  mets_red = compute_ranges(
    training_data |> filter(Flesh == "Red"),
    model_config$mets_red$variables
  )
)

# ── Helper: create template with example values from training medians ────────
create_template <- function(variables, ranges) {
  template <- data.frame(Sample_Name = c("Sample_1", "Sample_2"))
  for (var in variables) {
    median_val <- ranges$Median[ranges$Metabolite == var]
    template[[var]] <- round(rep(median_val, 2), 2)
  }
  template
}

# ── Metabolite descriptions (for tooltips / user reference) ──────────────────
metabolite_descriptions <- c(
  Sulcatone = "Ketone volatile (6-methyl-5-hepten-2-one); strong predictor of liking",
  Linalool_oxide = "Oxygenated terpene; linked to aroma intensity in papaya",
  Brix = "Total soluble solids; proxy for sweetness flavour"
)


# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- page_navbar(
  title = div(
    style = "font-size: 2rem; font-weight: bold;",
    "Papaya Flavour Predictor",
    nav_spacer()
  ),
  header = uiOutput("navbar_colour"),
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#FF8C00",
    success = "#498D50",
    warning = "#f4a261",
    danger = "#e63946"
  ) |> bs_add_rules("
    .navbar-collapse {
      flex-grow: 0 !important;
      margin-left: auto !important;
      margin-right: auto !important;
    }
    .navbar-nav .nav-link.active {
      color: #000000 !important;
    }
    .navbar-brand {
      color: white !important;
      text-shadow: -1px -1px 0 #7A3500, 1px -1px 0 #7A3500,
                   -1px  1px 0 #7A3500, 1px  1px 0 #7A3500;
    }
    .navbar {
      background-color: #E56717 !important;
    }
    /* Make every cell in the manual entry table look like an editable field */
    #manual_input_table table.dataTable td {
      border: 1px solid #dee2e6 !important;
      cursor: text;
      min-width: 80px;
    }
    #manual_input_table table.dataTable td:hover {
      background-color: #fff9f0 !important;
    }
  "),

  # ── Tab 1: Home ──────────────────────────────────────────────────────────
  nav_panel(
    "Home",
    card(
      card_header("Welcome to the Papaya Flavour Prediction App"),
      card_body(
        h4("What This App Does"),
        p("This app predicts consumer liking scores for papaya fruit based on metabolite concentrations, using an XGBoost machine learning model trained on sensory panel data from 125 untrained consumers across nine distinct papaya genotypes."),
        h4("Important Limitations"),
        p("Predictions are indicative rather than definitive and may not generalise beyond the genotypes, growing conditions, and analytical methods used in the original study. Metabolite concentrations should be collected using certified reference standards and appropriate internal reference standards following the published methodology — deviations may reduce accuracy. These models do not require matrix correction.
Predicted liking scores reflect the average preferences of the training panel and may not represent other consumer groups or demographics. This tool is intended to support research and decision-making, not replace sensory evaluation."),
        p("For full methodological details, see the published article: ",
          tags$a(
            "Metabolomic modelling of sensory characteristics and consumer liking in papaya fruit.",
            href = "https://doi.org/10.1016/j.foodchem.2026.148323",
            target = "_blank"
          )
        ),
        hr(),

        h4("Step-by-Step Instructions"),
        tags$ol(
          tags$li(strong("Select a Model:"), " Use the sidebar to choose
                  'All Papaya' or 'Red-Fleshed Only'."),
          tags$li(strong("Prepare Your Data:"), " Download the template below
                  and fill in your metabolite concentrations."),
          tags$li(strong("Enter Data:"), " Go to 'Data Input' and either type
                  values manually or upload your completed template."),
          tags$li(strong("Validate:"), " Click 'Validate Data' to check for
                  issues."),
          tags$li(strong("Predict:"), " Go to 'Predictions' and click
                  'Run Prediction'."),
          tags$li(strong("Download:"), " View results and download as Excel.")
        ),
        hr(),

        h4("Download Template"),
        p("Template with correct column names for the selected model:"),
        downloadButton("download_template", "Download Template (.xlsx)",
                       class = "btn-primary"),
        hr(),

        div(
          class = "alert alert-warning",
          role = "alert",
          icon("exclamation-triangle"), strong(" Important:"),
          " Values outside the model training range may reduce prediction
          accuracy. See 'About the Model' for valid ranges."
        )
      )
    )
  ),

  # ── Tab 2: Data Input ────────────────────────────────────────────────────
  nav_panel(
    "Data Input",
    card(
      card_header("Enter Metabolite Data"),
      card_body(
        navset_tab(
          nav_panel(
            "Manual Entry",
            br(),
            actionButton("add_row", "Add Row", icon = icon("plus"),
                         class = "btn-success"),
            actionButton("clear_table", "Clear All Rows", icon = icon("eraser"),
                         class = "btn-warning"),
            actionButton("load_example", "Load Example Data", icon = icon("table"),
                         class = "btn-info"),
            br(), br(),
            # Dropdown + button for removing a specific row by sample name
            div(
              style = "display: flex; align-items: center; gap: 8px; max-width: 420px;",
              div(style = "flex: 1;",
                  selectInput("row_to_remove", NULL, choices = NULL, width = "100%")),
              actionButton("remove_row", "Remove Row", icon = icon("minus"),
                           class = "btn-danger")
            ),
            br(),
            DTOutput("manual_input_table")
          ),
          nav_panel(
            "File Upload",
            br(),
            fileInput("file_upload", "Choose Excel File (.xlsx only)",
                      accept = ".xlsx"),
            uiOutput("upload_status"),
            br(),
            DTOutput("uploaded_data_preview")
          )
        )
      )
    ),
    card(
      card_header("Data Validation"),
      card_body(
        actionButton("validate_data", "Validate Data",
                     icon = icon("clipboard-check"),
                     class = "btn-primary"),
        br(), br(),
        uiOutput("validation_results")
      )
    )
  ),

  # ── Tab 3: Predictions ──────────────────────────────────────────────────
  nav_panel(
    "Predictions",
    card(
      card_header("Generate Predictions"),
      card_body(
        actionButton("run_prediction", "Run Prediction",
                     icon = icon("play-circle"),
                     class = "btn-primary btn-lg"),
        br(), br(),
        uiOutput("prediction_status"),
        br(),
        conditionalPanel(
          condition = "output.predictions_ready",
          downloadButton("download_results", "Download Results (.xlsx)",
                         class = "btn-success"),
          br(), br(),
          h4("Prediction Results"),
          DTOutput("results_table"),
          br(),
          h4("Visual Summary"),
          div(
            style = "width: 90%;",
            plotOutput("results_plot", height = "350px")
          )
        )
      )
    )
  ),

  # ── Tab 4: About the Model ──────────────────────────────────────────────
  nav_panel(
    "About the Model",
    card(
      card_header("About the Model"),
      card_body(
        fill = FALSE,

        h4("Acknowledgements"),
        p("This project was funded by Hort Innovation and the Australian Government, through Hort Frontiers Fund AS19003 and the Griffith University Postgraduate Research Scholarship. All analytical chemistry assays were undertaken with the trusted oversight of Alan White, Griffith University Analytical Facility Manager. All food sensory analysis was managed by Emma Hassal, Queensland Alliance for Agriculture and Food Innovation Nutrition and Food Science Senior Research Assistant and Laboratory Manager."),
        hr(),

        h4("Training Data Ranges"),
        p("Metabolite concentration ranges from the training data.
          Amber rows show where your input values fall outside these ranges."),
        gt_output("training_ranges_table"),
        hr(),

        h4("What is XGBoost?"),
        p("XGBoost (eXtreme Gradient Boosting) is a machine learning algorithm
          that combines many simple decision trees, each correcting the mistakes
          of the previous one. It has learned the relationship between metabolite
          concentrations and consumer liking scores from papaya sensory data."),
        hr(),

        h4("Model Performance"),
        layout_columns(
          col_widths = breakpoints(sm = 12, md = 4),
          uiOutput("model_r2_box"),
          uiOutput("model_rmse_box"),
          uiOutput("model_n_box")
        ),
        hr(),

        h4("Feature Importance"),
        p("How much each metabolite contributes to the prediction:"),
        plotOutput("importance_plot", height = "200px"),
        hr(),

        p(strong("Model Last Updated:"), "February 2026")
      )
    )
  ),

  # ── Tab 5: Flavour Wheel ─────────────────────────────────────────────────
  nav_panel(
    "Flavour Wheel",
    card(
      card_header("Papaya Flavour Wheel"),
      card_body(
        fill = FALSE,
        p("The papaya flavour wheel describes the range of aroma and taste
          descriptors identified in papaya fruit by trained sensory panellists. The inner ring groups the sensory categories; the middle ring are the traits used to score fruit samples; the outer ring descirbes the middle ring traits"),
        tags$div(
          style = "text-align: center;",
          img(src = "imgs/papaya_flavour_wheel.png",
              style = "max-width: 100%; height: auto;")
        )
      )
    )
  ),

  # ── Sidebar ─────────────────────────────────────────────────────────────
  sidebar = sidebar(
    title = "Controls",
    radioButtons(
      "model_choice",
      "Select Model:",
      choices = c(
        "All Papaya (Red + Yellow)" = "mets",
        "Red-Fleshed Only" = "mets_red"
      ),
      selected = "mets_red"
    ),
    hr(),
    actionButton(
      "reset_all",
      "Reset Everything",
      icon = icon("redo"),
      class = "btn-danger btn-block"
    ),
    hr(),
    h5("Progress Tracker"),
    uiOutput("progress_tracker"),
    hr(),
    tags$div(
      style = "text-align: center; color: #888;",
      tags$div(
        style = "display: flex; justify-content: center; align-items: center; gap: 12px; margin-bottom: 6px;",
        img(src = "imgs/Univiersity_logo.png", height = "35px"),
        img(src = "imgs/Funding_logo.png", height = "35px")
      ),
      
      tags$blockquote(
        "Papaya Flavour Research: a project funded by Hort Innovation grant",
        tags$a("Genetics of fruit sensory preferences (AS19003)", 
               href = "https://www.horticulture.com.au/growers/help-your-business-grow/research-reports-publications-fact-sheets-and-more/as19003/", 
               target = "_blank")
      ),
      tags$small(
        tags$a("Joshua Lomax", 
               href = "https://www.linkedin.com/in/joshmarkhillislomax/", 
               target = "_blank"),
        br(),
        tags$a("Rebecca Ford", 
               href = "https://experts.griffith.edu.au/18969-rebecca-ford", 
               target = "_blank"),
        br(),
        tags$a("Ido Bar", href = "http://experts.griffith.edu.au/8327-ido-bar", 
               target = "_blank")
      )
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Reactive values ────────────────────────────────────────────────────
  rv <- reactiveValues(
    data = NULL,
    validated = FALSE,
    validation_messages = NULL,
    predictions = NULL,
    data_source = "manual"
  )

  # ── Helpers to get current model config and ranges ─────────────────────
  current_config <- reactive({
    model_config[[input$model_choice]]
  })

  current_ranges <- reactive({
    training_ranges[[input$model_choice]]
  })

  # ── Dynamically update navbar colour based on model selection ──────────
  output$navbar_colour <- renderUI({
    color <- if (input$model_choice == "mets_red") "#e63946" else "#E56717"
    tags$style(HTML(paste0(".navbar { background-color: ", color, " !important; }")))
  })

  # ── Reset data when model selection changes (start with empty table) ───
  observeEvent(input$model_choice, {
    config <- current_config()
    empty_df <- data.frame(Sample_Name = character(0))
    for (var in config$variables) {
      empty_df[[var]] <- numeric(0)
    }
    rv$data <- empty_df
    rv$validated <- FALSE
    rv$validation_messages <- NULL
    rv$predictions <- NULL
  })

  # ── Progress tracker ───────────────────────────────────────────────────
  output$progress_tracker <- renderUI({
    has_data <- !is.null(rv$data) && nrow(rv$data) > 0
    tags$ul(
      class = "list-unstyled",
      tags$li(
        icon(if (has_data) "check-circle" else "circle",
             class = if (has_data) "text-success" else "text-muted"),
        " Data Entered"
      ),
      tags$li(
        icon(if (rv$validated) "check-circle" else "circle",
             class = if (rv$validated) "text-success" else "text-muted"),
        " Data Validated"
      ),
      tags$li(
        icon(if (!is.null(rv$predictions)) "check-circle" else "circle",
             class = if (!is.null(rv$predictions)) "text-success" else "text-muted"),
        " Predictions Generated"
      )
    )
  })

  # ── Download template ──────────────────────────────────────────────────
  output$download_template <- downloadHandler(
    filename = function() {
      paste0("papaya_template_", input$model_choice, ".xlsx")
    },
    content = function(file) {
      config <- current_config()
      ranges <- current_ranges()
      write.xlsx(create_template(config$variables, ranges), file)
    }
  )

  # ── Manual input table (editable) ──────────────────────────────────────
  output$manual_input_table <- renderDT({
    req(rv$data)
    datatable(
      rv$data,
      editable = list(target = "cell"),
      selection = "none",
      options = list(
        pageLength = nrow(rv$data),
        dom = "t",
        # Display NA cells as blank rather than the text "NA"
        columnDefs = list(list(defaultContent = "", targets = "_all"))
      ),
      rownames = FALSE
    )
  })

  # ── Handle cell edits ──────────────────────────────────────────────────
  observeEvent(input$manual_input_table_cell_edit, {
    info <- input$manual_input_table_cell_edit
    # DT uses 0-based column index; R data.frames use 1-based
    col_idx <- info$col + 1

    if (col_idx == 1) {
      # Sample_Name column stays as text
      rv$data[info$row, col_idx] <- as.character(info$value)
    } else {
      # Metabolite columns converted to numeric
      rv$data[info$row, col_idx] <- suppressWarnings(as.numeric(info$value))
    }

    # Any edit invalidates previous validation and predictions
    rv$validated <- FALSE
    rv$predictions <- NULL
  })

  # ── Add row ────────────────────────────────────────────────────────────
  observeEvent(input$add_row, {
    config <- current_config()
    new_row <- data.frame(
      Sample_Name = paste0("Sample_", nrow(rv$data) + 1)
    )
    for (var in config$variables) {
      new_row[[var]] <- NA_real_
    }
    rv$data <- rbind(rv$data, new_row)
    rv$validated <- FALSE
    rv$predictions <- NULL
  })

  # ── Keep the "remove row" dropdown in sync with current sample names ────
  observe({
    choices <- if (!is.null(rv$data) && nrow(rv$data) > 0) {
      rv$data$Sample_Name
    } else {
      character(0)
    }
    updateSelectInput(session, "row_to_remove", choices = choices)
  })

  # ── Remove the row matching the selected sample name ────────────────────
  observeEvent(input$remove_row, {
    selected_name <- input$row_to_remove
    if (!is.null(selected_name) && selected_name != "" && nrow(rv$data) > 0) {
      # Match on the first row with that sample name (handles any duplicates safely)
      row_idx <- which(rv$data$Sample_Name == selected_name)[1]
      if (!is.na(row_idx)) {
        rv$data <- rv$data[-row_idx, , drop = FALSE]
        rv$validated <- FALSE
        rv$predictions <- NULL
      }
    }
  })

  # ── Load example data from training set ────────────────────────────────
  # Filters to the correct samples for the selected model, then subsets to
  # only the model's metabolite columns. Row names become Sample_Name.
  observeEvent(input$load_example, {
    config <- current_config()

    example_data <- if (input$model_choice == "mets_red") {
      training_data |> filter(Flesh == "Red")
    } else {
      training_data
    }

    rv$data <- data.frame(
      Sample_Name = rownames(example_data),
      example_data[, config$variables, drop = FALSE],
      check.names = FALSE,
      row.names = NULL
    )
    rv$validated <- FALSE
    rv$predictions <- NULL
  })

  # ── Clear all rows (produces an empty table, keeping correct columns) ───
  observeEvent(input$clear_table, {
    config <- current_config()
    empty_df <- data.frame(Sample_Name = character(0))
    for (var in config$variables) {
      empty_df[[var]] <- numeric(0)
    }
    rv$data <- empty_df
    rv$validated <- FALSE
    rv$predictions <- NULL
  })

  # ── File upload ────────────────────────────────────────────────────────
  observeEvent(input$file_upload, {
    req(input$file_upload)

    tryCatch({
      uploaded <- read.xlsx(input$file_upload$datapath)
      rv$data <- uploaded
      rv$data_source <- "upload"
      rv$validated <- FALSE
      rv$predictions <- NULL

      output$upload_status <- renderUI({
        div(
          class = "alert alert-success",
          icon("check-circle"),
          sprintf(" File uploaded: %d rows, %d columns.",
                  nrow(uploaded), ncol(uploaded))
        )
      })
    }, error = function(e) {
      output$upload_status <- renderUI({
        div(
          class = "alert alert-danger",
          icon("exclamation-circle"),
          paste(" Error reading file:", e$message)
        )
      })
    })
  })

  # ── Preview of uploaded data ───────────────────────────────────────────
  output$uploaded_data_preview <- renderDT({
    req(input$file_upload, rv$data)
    if (rv$data_source != "upload") return(NULL)
    datatable(rv$data,
              options = list(pageLength = 10, dom = "tp"),
              rownames = FALSE)
  })

  # ── Validation logic ───────────────────────────────────────────────────
  observeEvent(input$validate_data, {
    req(rv$data)
    config <- current_config()
    errors <- character()

    # 1. Check Sample_Name column exists
    if (!"Sample_Name" %in% names(rv$data)) {
      errors <- c(errors, "Missing required column: 'Sample_Name'")
    }

    # 2. Check all required metabolite columns exist
    missing_cols <- setdiff(config$variables, names(rv$data))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste("Missing metabolite columns:",
                                paste(missing_cols, collapse = ", ")))
    }

    # 3. Only run value-level checks if column structure is correct
    if (length(errors) == 0) {

      # Duplicate sample names
      if (any(duplicated(rv$data$Sample_Name))) {
        dupes <- unique(rv$data$Sample_Name[duplicated(rv$data$Sample_Name)])
        errors <- c(errors, paste("Duplicate sample names:",
                                  paste(dupes, collapse = ", ")))
      }

      # Empty sample names
      empty_names <- which(is.na(rv$data$Sample_Name) | rv$data$Sample_Name == "")
      if (length(empty_names) > 0) {
        errors <- c(errors, paste("Empty Sample_Name at row(s):",
                                  paste(empty_names, collapse = ", ")))
      }

      # Check each metabolite column for missing or non-numeric values
      for (var in config$variables) {
        col_vals <- rv$data[[var]]

        if (!is.numeric(col_vals)) {
          numeric_vals <- suppressWarnings(as.numeric(col_vals))
          non_numeric <- which(is.na(numeric_vals) & !is.na(col_vals))
          if (length(non_numeric) > 0) {
            errors <- c(errors, sprintf(
              "Non-numeric values in '%s' at row(s): %s",
              var, paste(non_numeric, collapse = ", ")
            ))
          }
          col_vals <- numeric_vals
        }

        na_rows <- which(is.na(col_vals))
        if (length(na_rows) > 0) {
          errors <- c(errors, sprintf(
            "Empty values in '%s' at row(s): %s",
            var, paste(na_rows, collapse = ", ")
          ))
        }
      }
    }

    # Store validation result
    if (length(errors) == 0) {
      rv$validated <- TRUE
      rv$validation_messages <- list(
        status = "success",
        messages = "All checks passed! You can now run predictions."
      )
    } else {
      rv$validated <- FALSE
      rv$validation_messages <- list(status = "error", messages = errors)
    }
  })

  # ── Display validation results ─────────────────────────────────────────
  output$validation_results <- renderUI({
    req(rv$validation_messages)

    if (rv$validation_messages$status == "success") {
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(" Validation Passed: "),
        rv$validation_messages$messages
      )
    } else {
      div(
        class = "alert alert-danger",
        icon("exclamation-circle"),
        strong(" Validation Failed:"),
        tags$ul(lapply(rv$validation_messages$messages, tags$li))
      )
    }
  })

  # ── Run prediction using real XGBoost model ────────────────────────────
  observeEvent(input$run_prediction, {
    # Warn if data hasn't been validated yet
    if (!rv$validated) {
      output$prediction_status <- renderUI({
        div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          " Please validate your data first on the 'Data Input' tab."
        )
      })
      return()
    }

    config <- current_config()

    tryCatch({
      # Build prediction matrix with columns in the same order as training
      prediction_matrix <- rv$data |>
        select(all_of(config$variables)) |>
        mutate(across(everything(), as.numeric)) |>
        as.matrix()

      # Create xgb.DMatrix and run prediction
      dmatrix <- xgb.DMatrix(data = prediction_matrix)
      predicted_scores <- predict(config$model, dmatrix)

      # Store results
      rv$predictions <- data.frame(
        Sample_Name = rv$data$Sample_Name,
        Predicted_Liking_Score = round(predicted_scores, 2)
      )

      output$prediction_status <- renderUI({
        div(
          class = "alert alert-success",
          icon("check-circle"),
          sprintf(" Predictions generated for %d samples.",
                  nrow(rv$predictions))
        )
      })
    }, error = function(e) {
      rv$predictions <- NULL
      output$prediction_status <- renderUI({
        div(
          class = "alert alert-danger",
          icon("exclamation-circle"),
          paste(" Prediction error:", e$message)
        )
      })
    })
  })

  # ── Flag for conditionalPanel to show/hide results ─────────────────────
  output$predictions_ready <- reactive({ !is.null(rv$predictions) })
  outputOptions(output, "predictions_ready", suspendWhenHidden = FALSE)

  # ── Results table with red-to-green colour gradient ────────────────────
  output$results_table <- renderDT({
    req(rv$predictions)

    datatable(
      rv$predictions,
      options = list(
        pageLength = nrow(rv$predictions),
        dom = "t",
        order = list(list(1, "desc"))
      ),
      rownames = FALSE
    ) |>
      formatStyle(
        "Predicted_Liking_Score",
        backgroundColor = styleInterval(
          cuts = c(40, 57, 61, 65, 69),
          values = c("#e63946", "#f4845f", "#f4a261",
                     "#e9c46a", "#a7c957", "#52b788")
        )
      )
  })

  # ── Bar chart of predicted scores ──────────────────────────────────────
  output$results_plot <- renderPlot({
    req(rv$predictions)

    rv$predictions |>
      mutate(Sample_Name = reorder(Sample_Name, Predicted_Liking_Score)) |>
      ggplot(aes(x = Sample_Name, y = Predicted_Liking_Score,
                 fill = Predicted_Liking_Score)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "#e63946", high = "#52b788") +
      coord_flip() +
      labs(x = NULL,
           y = "Predicted Consumer Liking Score",
           title = "Predicted Liking Scores by Sample") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })

  # ── Download results as timestamped .xlsx ──────────────────────────────
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("papaya_predictions_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      write.xlsx(rv$predictions, file)
    }
  )

  # ── Training ranges table with amber highlighting for out-of-range ─────
  output$training_ranges_table <- render_gt({
    ranges <- current_ranges() |>
      mutate(across(where(is.numeric), \(x) round(x, 2))) |>
      select(Metabolite, Min, Max, Units)

    # Check if any user input values fall outside training ranges
    config <- current_config()
    if (!is.null(rv$data) && all(config$variables %in% names(rv$data))) {
      out_of_range <- sapply(config$variables, function(var) {
        vals <- suppressWarnings(as.numeric(rv$data[[var]]))
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) return(FALSE)
        range_row <- ranges[ranges$Metabolite == var, ]
        any(vals < range_row$Min | vals > range_row$Max)
      }) |> as.logical()
      ranges <- ranges |> mutate(flag = out_of_range[Metabolite])
    } else {
      ranges <- ranges |> mutate(flag = FALSE)
    }

    ranges |>
      gt() |>
      tab_style(
        style = cell_fill(color = "#f4a261"),
        locations = cells_body(rows = flag)
      ) |>
      cols_hide(flag) |> 
      tab_options(
        column_labels.padding.horizontal = px(20),
        data_row.padding.horizontal = px(20)
      )
  })

  # ── Model performance value boxes ──────────────────────────────────────
  output$model_r2_box <- renderUI({
    config <- current_config()
    value_box(
      title = "R-squared",
      value = sprintf("%.2f", config$r_squared),
      showcase = icon("chart-line"),
      theme = "success",
      p(sprintf("Explains %.0f%% of variance", config$r_squared * 100))
    )
  })

  output$model_rmse_box <- renderUI({
    config <- current_config()
    value_box(
      title = "RMSE",
      value = sprintf("%.2f", config$rmse),
      showcase = icon("bullseye"),
      theme = "primary",
      p("Average prediction error (test set)")
    )
  })

  output$model_n_box <- renderUI({
    config <- current_config()
    value_box(
      title = "Training Samples",
      value = config$n_train,
      showcase = icon("database"),
      theme = "warning",
      p(sprintf("Train: %d | Test: %d", config$n_train, config$n_test))
    )
  })

  # ── Feature importance plot ────────────────────────────────────────────
  output$importance_plot <- renderPlot({
    config <- current_config()

    config$importance |>
      as.data.frame() |>
      mutate(Feature = reorder(Feature, Gain)) |>
      ggplot(aes(x = Feature, y = Gain, fill = Gain)) +
      geom_col(width = 0.6) +
      scale_fill_gradient(low = "#f4a261", high = "#FF8C00") +
      coord_flip() +
      labs(x = NULL,
           y = "Feature Importance (Gain)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            plot.margin = margin(5, 10, 5, 5))
  })

  # ── Reset everything ───────────────────────────────────────────────────
  observeEvent(input$reset_all, {
    config <- current_config()
    empty_df <- data.frame(Sample_Name = character(0))
    for (var in config$variables) {
      empty_df[[var]] <- numeric(0)
    }
    rv$data <- empty_df
    rv$validated <- FALSE
    rv$validation_messages <- NULL
    rv$predictions <- NULL
    output$upload_status <- renderUI(NULL)
    output$prediction_status <- renderUI(NULL)
  })
}


# ── Run the app ──────────────────────────────────────────────────────────────
shinyApp(ui, server)
