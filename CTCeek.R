WWW_DIR  <- file.path("./www")
CTC_DIR  <- file.path("./")
REF_RDS  <- file.path("./pbmcScactcSeurat.withModel.rds")

shiny::addResourcePath("static", WWW_DIR)

# --- 1. Library loading ---
library(shiny)
library(shinycssloaders)
library(shiny.emptystate)
library(Seurat)
library(dplyr)
library(DT)
library(ggplot2)
library(shinyjs)
library(patchwork)
library(promises)
library(future)

# ---- Utility ----
safe_npcs <- function(object, assay = DefaultAssay(object), requested = 50, min_pcs = 2) {
  mat <- tryCatch(GetAssayData(object, assay = assay, slot = "scale.data"), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    nf <- length(VariableFeatures(object, assay = assay))
    nc <- ncol(object)
    max_rank <- max(0, min(nf, nc) - 1L)
  } else {
    max_rank <- max(0, min(nrow(mat), ncol(mat)) - 1L)
  }
  max_allowed <- max(min_pcs, max_rank)
  as.integer(min(requested, max_allowed))
}

# ---- Utility----
pick_ref_label_col <- function(ref_obj) {
  candidates <- c(
    "celltype.l2","predicted.celltype.l2","celltype.l1","predicted.celltype.l1",
    "seurat_annotations","celltype","CellType","cell_type","celltypes","major.celltype",
    "annotation","annotations","label","labels","SubCellType","subclass","class"
  )
  hit <- intersect(candidates, colnames(ref_obj@meta.data))
  if (length(hit) > 0) return(hit[1])
  idvec <- tryCatch(Idents(ref_obj), error = function(e) NULL)
  if (!is.null(idvec)) { ref_obj$ctcscan_idents <- as.character(idvec); return("ctcscan_idents") }
  ref_obj$ctcscan_dummy <- "Reference"; "ctcscan_dummy"
}

# --- 2. UI ---
ui <- fluidPage(
  tags$style(HTML("/* Basic Body and Font Styling */
    html {
      font-size: 16px;
    }
    body {
      background-image: url('background_red.png');
      background-size: cover;
      background-position: center center;
      background-repeat: no-repeat;
      background-attachment: fixed;
      font-family: 'Open Sans', sans-serif;
      color: #333;
      margin: 0;
      padding: 0;
      overflow-x: hidden;
    }

    /* Logo Section */
    #logo-section {
      text-align: center;
      padding: 30px 20px;
      background-color: rgba(245, 218, 218, 0.5);
      margin-bottom: 20px;
      border-bottom: 1px solid transparent;
    }
    #logo-section img {
      max-width: 150%;
      height: 500px;
      object-fit: contain;
    }

    /* Main Content Container (Card-like) */
    .app-card {
      background-color: rgba(255, 255, 255, 1);
      padding: 30px;
      border-radius: 15px;
      box-shadow: 0 8px 16px rgba(0, 0, 0, 0.2);
      margin: 20px auto;
      max-width: 1200px; /* Made card wider for two plots */
      width: 95%;
      box-sizing: border-box;
      animation: fadeIn 1s ease-out;
    }

    /* Animation for fade-in */
    @keyframes fadeIn {
      from { opacity: 0; transform: translateY(20px); }
      to { opacity: 1; transform: translateY(0); }
    }

    /* Section Headings within the card */
    .app-card h2 {
      text-align: center;
      font-weight: 600;
      color: #b30021;
      margin-bottom: 25px;
      padding-bottom: 10px;
      border-bottom: 2px solid #f5dada;
    }
    .app-card h3 {
      font-weight: 500;
      color: #b30021;
      margin-top: 30px;
      margin-bottom: 15px;
      border-bottom: 1px solid #eee;
      padding-bottom: 5px;
    }
    .app-card h4 {
      font-weight: 600;
      color: #b30021;
      margin-top: 25px;
      margin-bottom: 15px;
      border-bottom: 1px solid #f5dada;
      padding-bottom: 5px;
      font-size: 1.3em;
    }

    /* Paragraphs and general text */
    .app-card p {
      font-size: 1rem;
      line-height: 1.6;
      margin-bottom: 15px;
    }
    .app-card p.contact-info {
        font-size: 0.9rem;
        color: #f5dada;
    }

    /* Info/Contact Boxes */
    .info-panel {
      background-color: #f5dada;
      border: 1px solid #f5dada;
      border-radius: 10px;
      padding: 20px;
      margin-bottom: 25px;
      box-shadow: 0 2px 4px rgba(0,0,0,0.05);
      display: flex;
      align-items: flex-start;
      flex-wrap: wrap;
    }
    .info-panel .icon {
      font-size: 2.5em;
      color: #b30021;
      margin-right: 15px;
      flex-shrink: 0;
    }
    .info-panel .text-content {
      flex-grow: 1;
      min-width: 0;
    }
    .info-panel h4.panel-title {
        margin-top: 0;
        margin-bottom: 10px;
        color: #b30021;
        font-weight: 600;
        font-size: 1.2rem;
        border-bottom: none;
        padding-bottom: 0;
    }
    .info-panel p {
        font-size: 0.95rem;
        margin-bottom: 5px;
    }
    .info-panel p strong {
        color: #b30021;
    }
    .info-panel a {
        color: #b30021;
        text-decoration: none;
    }
    .info-panel a:hover {
        text-decoration: underline;
    }

    /* File Input Styling */
    .form-group.shiny-input-container {
      margin-bottom: 25px;
    }
    .btn-file {
      background-color: #b30021;
      color: white;
      border-color: #b30021;
      transition: background-color 0.3s ease;
    }
    .btn-file:hover {
      background-color: #b30021;
      border-color: #b30021;
    }
    .umap-container { border-radius:15px; box-shadow:0 4px 8px rgba(0,0,0,0.1); background:#f8f9fa; padding:15px; margin-top:25px; }
    #footer { text-align:center; padding:10px 20px; font-size:13px; color:#6c757d; background-color: rgba(245,218,218,0.9); margin-top:40px; }
    .footer-logo { max-width:280px; height:auto; margin:5px 20px; display:inline-block; vertical-align:middle; }
  ")),
  div(id = "logo-section", img(src = "static/CTCeek_logo.png", alt = "CTCeek Logo")),
  div(class = "app-card",
      h2("Welcome to CTCeek"),
      fluidRow(
        column(6,
               div(class = "info-panel",
                   tags$div(class = "icon", tags$i(class = "fas fa-microscope")),
                   tags$div(class = "text-content",
                            tags$h4(class = "panel-title", "Software Description"),
                            tags$p("This software automatically annotates CTCs in your scRNA-seq data.")
                   )
               )
        ),
        column(6,
               div(class = "info-panel",
                   tags$div(class = "icon", tags$i(class = "fas fa-envelope")),
                   tags$div(class = "text-content",
                            tags$h4(class = "panel-title", "Contact Us"),
                            tags$p("Prof. Stefano Volinia: ", tags$strong(tags$a(href = "mailto:s.volinia@unife.it", "s.volinia@unife.it"))),
                            tags$p("Pietro Ancona: ", tags$strong(tags$a(href = "mailto:ncnptr@unife.it", "ncnptr@unife.it")))
                   )
               )
        )
      ),
      br(),
      
      h3("1. Upload Your Seurat Object"),
      fileInput("file", "Upload your Seurat object (.rds)", accept = ".rds"),
      uiOutput("empty_state"),
      uiOutput("qc_ui"),
      plotOutput("qc_feature_plot"),
      plotOutput("qc_umi_plot"),
      plotOutput("qc_mt_plot"),
      h3("2. Start Analysis"),
      div(style="text-align:center;", actionButton("Start", "Start analyze", class="btn-lg btn-success")),
      br(),
      
      # ======= 3. Results & Download =======
      h3("3. View Analysis & Download Results"),
      withSpinner(uiOutput("dynamic_table"), type = 8, color = "#007bff"),
      uiOutput("ctc_summary"),
      div(style="text-align:center; margin-top:10px;",
          downloadButton("download_predictions_csv", "Download Predictions (CSV)"),
          downloadButton("download_seurat_obj", "Download Annotated Seurat Object")
      ),
      uiOutput("markers_zip_download_ui"),  # appare dopo il calcolo marker
      br(),
      
      conditionalPanel(
        condition = "output.analysis_complete",
        h3("4. Visualize Data on Reference UMAP"),
        div(id="umap_button_div",
            div(style="text-align:center; margin-top:20px;",
                actionButton("compute_umap","Project onto Reference UMAP", class="btn-primary"),
                helpText("Click to project your data onto the reference UMAP.")
            )),
        shinyjs::hidden(
          div(id="umap_spinner_div", withSpinner(plotOutput("dummy_spinner_plot_1", height="100px"), type=8),
              p("Finding anchors and projecting data onto reference UMAP, please wait..."))),
        shinyjs::hidden(
          div(id="umap_plot_div",
              div(class="umap-container",
                  plotOutput("combined_umap_plot", height=500, width="100%"),
                  downloadButton("download_umap_plot","Download UMAP Plot"))))
      ),
      
      # ======= 5. Marker one-vs-all per cell_type =======
      conditionalPanel(
        condition = "output.analysis_complete",
        h3("5. Find Marker Genes per Cell Type (one-vs-all)"),
        p("Seleziona i ceppi da confrontare contro tutti gli altri (one-vs-all)."),
        uiOutput("group_selector_ui"),
        div(id="markers_button_div",
            div(style="text-align:center; margin-top:10px;",
                actionButton("find_markers_groups","Calculate Selected Markers", class="btn-primary"))),
        shinyjs::hidden(
          div(id="markers_spinner_div",
              withSpinner(plotOutput("dummy_spinner_plot_2", height="90px"), type=8),
              p("Calculating per-group marker genes, please wait..."))),
        uiOutput("per_group_marker_tabs")  # tabs dinamiche, una per ceppo
      )
  ),
  div(id="footer",
      img(src="static/pnrr-logo.jpg", alt="PNRR Logo", class="footer-logo"),
      img(src="static/medtras.png",  alt="MedTras Logo", class="footer-logo"),
      p("Volinia Lab — Department of Translational Medicine, University of Ferrara, Italy."))
)

# --- 3. Server ---
server <- function(input, output, session) {
  options(shiny.maxRequestSize   = 4000 * 1024^2)
  options(future.globals.maxSize = 4000 * 1024^2)
  
  TARGET_CTC_TYPES <- c("CTCepiA", "CTCepiB", "CTCmes")
  
  # --- Ref loading ---
  showNotification("Loading and preparing reference data...", id="ref_load", duration=NULL, type="message")
  if (!file.exists(REF_RDS)) { removeNotification("ref_load"); stop(paste0("File di riferimento mancante: ", REF_RDS)) }
  
  reference_seurat_obj <- readRDS(REF_RDS)
  reference_seurat_obj <- UpdateSeuratObject(reference_seurat_obj)
  DefaultAssay(reference_seurat_obj) <- "SCT"
  
  ref_label_col <- pick_ref_label_col(reference_seurat_obj)
  
  npcs_ref <- safe_npcs(reference_seurat_obj, assay="SCT", requested=50)
  reference_seurat_obj <- RunPCA(reference_seurat_obj, assay="SCT", verbose=FALSE, npcs=npcs_ref)
  
  if (!"umap.scanorama" %in% names(reference_seurat_obj@reductions)) {
    removeNotification("ref_load"); stop("Error: 'umap.scanorama' reduction not found in reference object.")
  }
  
  reference_seurat <- reactiveVal(reference_seurat_obj)
  reference_label_col <- reactiveVal(ref_label_col)
  
  removeNotification("ref_load")
  showNotification("App is ready.", type="message", duration=5)
  
  # --- Reactive & reset ---
  raw_seurat_object        <- reactiveVal(NULL)
  processed_seurat_object  <- reactiveVal(NULL)
  mapped_seurat_object     <- reactiveVal(NULL)
  is_analysis_complete     <- reactiveVal(FALSE)
  analyze_trigger          <- reactiveVal(0)
  selected_gene_for_plot   <- reactiveVal(NULL)
  
  # per group markers and selected
  per_group_markers_rv     <- reactiveVal(NULL)   # named list di data.frame
  selected_groups_rv       <- reactiveVal(NULL)
  
  reset_umap_view <- function(){
    mapped_seurat_object(NULL)
    shinyjs::hide("umap_plot_div");
    shinyjs::hide("umap_spinner_div");
    shinyjs::show("umap_button_div") }
  reset_markers_view <- function(){
    per_group_markers_rv(NULL); selected_groups_rv(NULL)
    shinyjs::hide("markers_spinner_div"); shinyjs::show("markers_button_div")
    output$per_group_marker_tabs <- renderUI({ NULL })
  }
  
  output$empty_state <- renderUI({
    if (is.null(input$file)) {
      div(style="text-align:center; padding:20px 0;",
          shiny.emptystate::empty_state_component(
            title="NO DATA FOUND. Please upload a Seurat object to get started.",
            content=tags$img(src="static/undraw_no-data_ig65.png", style="height:20rem;")
          ))
    }
  })
  
  observeEvent(input$file, {
    req(input$file)
    raw_seurat_object(tryCatch({ readRDS(input$file$datapath) }, error=function(e) NULL))
    is_analysis_complete(FALSE); analyze_trigger(0); reset_umap_view(); reset_markers_view()
  })
  
  observeEvent(input$file, {
    req(input$file)
    
    showNotification("Loading user data and running QC checks...", id="qc_load", duration=NULL, type="message")
    
    # Load user object
    raw_obj <- tryCatch({ readRDS(input$file$datapath) }, error=function(e) NULL)
    
    if (is.null(raw_obj)) {
      removeNotification("qc_load")
      showNotification("Error: cannot read RDS file.", type="error")
      return(NULL)
    }
    
    # Store raw object
    raw_seurat_object(raw_obj)
    
    # Reset states
    is_analysis_complete(FALSE)
    analyze_trigger(0)
    reset_umap_view()
    reset_markers_view()
    
    # --- QC Metrics ---
    # Calculate nFeature_RNA, nCount_RNA, and % mitochondrial
    gene_names <- rownames(raw_obj)
    mt_pattern <- NULL
    if (any(grepl("^MT-", gene_names))) {
      mt_pattern <- "^MT-"
    } else if (any(grepl("^mt-", gene_names))) {
      mt_pattern <- "^mt-"
    } else if (any(grepl("^Mt-", gene_names))) {
      mt_pattern <- "^Mt-"
    } else if (any(grepl("^MT_", gene_names))) {
      mt_pattern <- "^MT_"
    } else if (any(grepl("^mtr", gene_names, ignore.case = TRUE))) {
      mt_pattern <- "(?i)^mtr"  # case-insensitive "mtr"
    }
    
    if (!is.null(mt_pattern)) {
      raw_obj[["percent.mt"]] <- PercentageFeatureSet(raw_obj, pattern = mt_pattern)
      message(paste("Detected mitochondrial genes using pattern:", mt_pattern))
    } else {
      raw_obj[["percent.mt"]] <- 0
      message("No mitochondrial genes detected — check gene naming convention.")
    }
    
    # --- Render QC plots ---
    output$qc_feature_plot <- renderPlot({
      FeatureScatter(raw_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 2) +
        ggtitle("QC: nCount_RNA vs nFeature_RNA")
    })
    
    output$qc_umi_plot <- renderPlot({
      VlnPlot(raw_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.5) +
        ggtitle("QC: feature & UMI distributions")
    })
    
    output$qc_mt_plot <- renderPlot({
      VlnPlot(raw_obj, features = "percent.mt", pt.size = 0.5) +
        ggtitle("QC: % mitochondrial content")
    })
    
    # Optionally show some summary stats
    output$qc_ui <- renderUI({
      n_cells <- ncol(raw_obj)
      n_genes <- nrow(raw_obj)
      avg_genes <- round(mean(raw_obj$nFeature_RNA))
      avg_umis  <- round(mean(raw_obj$nCount_RNA))
      avg_mt    <- round(mean(raw_obj$percent.mt), 2)
      
      tagList(
        tags$h4("Quality Control Summary"),
        tags$p(paste("Cells:", n_cells)),
        tags$p(paste("Genes:", n_genes)),
        tags$p(paste("Avg genes per cell:", avg_genes)),
        tags$p(paste("Avg UMIs per cell:", avg_umis)),
        tags$p(paste("Avg % mitochondrial:", avg_mt, "%"))
      )
    })
    
    removeNotification("qc_load")
    showNotification("QC complete. Review plots before proceeding.", type="message", duration=5)
  })
  
  observeEvent(input$Start, {
    req(raw_seurat_object())
    is_analysis_complete(FALSE); reset_umap_view(); reset_markers_view(); analyze_trigger(analyze_trigger()+1)
  })
  
  output$dummy_spinner_plot_1 <- renderPlot({ NULL })
  output$dummy_spinner_plot_2 <- renderPlot({ NULL })
  
  # --- Mean analysis (SCT + PCA) ---
  observeEvent(analyze_trigger(), {
    if (analyze_trigger() > 0) {
      seurat_obj <- raw_seurat_object(); req(seurat_obj)
      showNotification("Starting Query Preprocessing (SCTransform, PCA)...", id="query_prep_notif", type="message", duration=NULL)
      
      future_promise({
        if (!"SCT" %in% Seurat::Assays(seurat_obj)) {
          seurat_obj <- SCTransform(seurat_obj, assay="RNA", verbose=FALSE)
        }
        DefaultAssay(seurat_obj) <- "SCT"
        
        # --- Logic from First Snippet ---
        max_pcs <- min(ncol(seurat_obj), nrow(seurat_obj)) 
        n_pcs <- min(50, max_pcs - 1)   # More constrained PCA calculation
        if (n_pcs < 2) stop("Too few cells/features to compute PCA. Please provide a larger object.")
        seurat_obj <- RunPCA(seurat_obj, assay="SCT", verbose=FALSE, npcs=n_pcs)
        # --- End Logic from First Snippet ---
        
        seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose=FALSE)
        return(seurat_obj)
      }, seed=TRUE) %...>% {
        processed_seurat_object(.)
        showNotification("Query Preprocessing Complete!", type="message", duration=5)
        is_analysis_complete(TRUE)
        removeNotification("query_prep_notif")
      } %...!% (function(e) {
        showNotification(paste("Preprocessing Error:", e$message), type="error", duration=NULL)
        is_analysis_complete(FALSE)
        removeNotification("query_prep_notif")
      })
    }
  }, ignoreInit=TRUE)
  
  output$analysis_complete <- reactive({ is_analysis_complete() })
  outputOptions(output, "analysis_complete", suspendWhenHidden=FALSE)
  
  
  # --- Reactive trigger for automatic UMAP ---
  umap_trigger <- reactiveVal(0)
  
  # Auto-trigger when analysis completes
  observeEvent(is_analysis_complete(), {
    if (is_analysis_complete()) {
      Sys.sleep(0.5)  # small delay to ensure processed_seurat_object is ready
      umap_trigger(umap_trigger() + 1)
    }
  }, ignoreInit = TRUE)
  
  # Manual button click also triggers
  observeEvent(input$compute_umap, {
    umap_trigger(umap_trigger() + 1)
  })
  
  # The actual UMAP computation
  observeEvent(umap_trigger(), {
    if (umap_trigger() == 0) return()  # Skip initial value
    
    query_seurat <- processed_seurat_object(); req(query_seurat)
    ref_obj  <- reference_seurat(); req(ref_obj)
    ref_col  <- reference_label_col()
    
    shinyjs::hide("umap_button_div"); shinyjs::show("umap_spinner_div")
    showNotification("Starting UMAP projection...", id="umap_proj_notif", type="message", duration=NULL)
    
    future_promise({
      # --- SAFER logic ---
      
      # 1. Ensure common features exist
      common_sct_features <- intersect(VariableFeatures(ref_obj), VariableFeatures(query_seurat))
      if (length(common_sct_features) < 200) {
        stop(paste("Low feature overlap between reference and query:", length(common_sct_features)))
      }
      
      # 2. Safe k-values
      k_score_safe  <- min(30, ncol(query_seurat) - 1, ncol(ref_obj) - 1)
      k_weight_safe <- min(20, ncol(query_seurat) - 1)
      
      # 3. Safe dimension handling
      pcs_ref   <- if ("pca" %in% Reductions(ref_obj)) ncol(Embeddings(ref_obj[["pca"]])) else 30
      pcs_query <- if ("pca" %in% Reductions(query_seurat)) ncol(Embeddings(query_seurat[["pca"]])) else 30
      max_dims  <- min(50, pcs_ref, pcs_query)
      if (max_dims < 5) stop("Too few shared PCs for projection.")
      
      # 4. Build anchors safely
      anchors <- FindTransferAnchors(
        reference = ref_obj,
        query = query_seurat,
        normalization.method = "SCT",
        recompute.residuals = TRUE,
        features = common_sct_features,
        reference.reduction = "pca",
        dims = 1:max_dims,
        k.score = k_score_safe
      )
      
      # 5. Validate UMAP model in reference
      umap_model_name <- if ("umap.scanorama" %in% Reductions(ref_obj)) "umap.scanorama" else "umap"
      if (is.null(ref_obj[[umap_model_name]]@misc$model)) {
        message("⚠️ Reference UMAP model not found, recomputing transferable model...")
        ref_obj <- RunUMAP(ref_obj, reduction = "pca", dims = 1:max_dims, return.model = TRUE)
        umap_model_name <- "umap"
      }
      
      # 6. Map query onto reference
      mapped_query <- MapQuery(
        anchorset = anchors,
        reference = ref_obj,
        query = query_seurat,
        refdata = setNames(list(ref_obj[[ref_col]][,1]), ref_col),
        reference.reduction = "pca",
        reduction.model = umap_model_name,
        transferdata.args = list(k.weight = k_weight_safe)
      )
      
      # 7. Handle predicted label
      possible_names <- c(ref_col, paste0("predicted.", ref_col))
      actual_prediction_col <- intersect(possible_names, colnames(mapped_query@meta.data))
      if (length(actual_prediction_col) == 0) {
        stop(paste("Could not find transferred metadata column. Looked for:",
                   paste(possible_names, collapse = " or ")))
      }
      mapped_query$CTCeek_final_prediction <- mapped_query[[actual_prediction_col[1]]]
      
      mapped_query
    }, seed=TRUE) %...>% {
      removeNotification("umap_proj_notif")
      showNotification("Projection complete!", type="message", duration=4)
      mapped_seurat_object(.) # store final object
      shinyjs::hide("umap_spinner_div"); shinyjs::show("umap_plot_div")
    } %...!% (function(e) {
      removeNotification("umap_proj_notif")
      showNotification(paste("Error during UMAP projection:", e$message), type="error", duration=NULL)
      reset_umap_view()
    })
  })
  
  
  # --- Plot UMAP paired ---
  combined_umap_plot_reactive <- reactive({
    req(processed_seurat_object())
    query_mapped <- mapped_seurat_object()
    req("ref.umap" %in% names(query_mapped@reductions))
    ref_obj <- reference_seurat(); req(ref_obj)
    ref_col <- reference_label_col()
    
    p1 <- DimPlot(ref_obj, reduction="umap.scanorama", group.by=ref_col, label=TRUE, repel=TRUE) +
      NoLegend() + ggtitle("Reference Cells")
    p2 <- DimPlot(query_mapped, reduction="ref.umap", group.by="CTCeek_final_prediction",
                  pt.size=1.5, label=TRUE, repel=TRUE) +
      ggtitle("Query Cells Projected on Reference")
    p1 + p2
  })
  output$combined_umap_plot <- renderPlot({ combined_umap_plot_reactive() })
  output$download_umap_plot <- downloadHandler(
    filename = function() { "CTCeek_UMAP.png" },
    content  = function(file) { ggsave(file, plot = combined_umap_plot_reactive(), width=14, height=7, dpi=300) }
  )
  
  # ======= 3. prediction table + download =======
  output$dynamic_table <- renderUI({
    req(is_analysis_complete(), mapped_seurat_object())
    tagList(h4("Cell Metadata with Predictions"), DTOutput("prediction_table"))
  })
  
  output$prediction_table <- renderDT({
    req(mapped_seurat_object())
    md <- mapped_seurat_object()@meta.data
    if (!"CTCeek_final_prediction" %in% colnames(md)) return(NULL)
    display_df <- data.frame(Cell_Barcode = rownames(md), Predicted_Cell_Type = md$CTCeek_final_prediction)
    datatable(display_df, options = list(pageLength=10, scrollX=TRUE), rownames=FALSE)
  })
  
  output$download_predictions_csv <- downloadHandler(
    filename = function() { "CTCeek_Predictions.csv" },
    content  = function(file) {
      req(mapped_seurat_object())
      md <- mapped_seurat_object()@meta.data
      df <- data.frame(Cell_Barcode = rownames(md), Predicted_Cell_Type = md$CTCeek_final_prediction)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # Info CTC
  output$ctc_summary <- renderUI({
    req(is_analysis_complete(), mapped_seurat_object())
    md <- mapped_seurat_object()@meta.data
    n <- sum(md$CTCeek_final_prediction %in% TARGET_CTC_TYPES, na.rm = TRUE)
    if (n > 0) {
      div(style="font-weight:bold; font-size:1.1rem; padding:10px; margin-top:10px; background:#d4edda; color:#155724; border-radius:8px; text-align:center;",
          HTML(paste0("<strong>", n, " CTC(s) found!</strong>")))
    } else {
      div(style="font-weight:bold; font-size:1.1rem; padding:10px; margin-top:10px; background:#fff3cd; color:#856404; border-radius:8px; text-align:center;",
          "No CTCs detected.")
    }
  })
  
  output$download_seurat_obj <- downloadHandler(
    filename = function() { "CTCeek_Seurat_Object.rds" },
    content  = function(file) {
      req(mapped_seurat_object())
      saveRDS(mapped_seurat_object(), file = file) }
  )
  
  # ======= 5. Marker one-vs-all =======
  
  # UI: selector for group
  output$group_selector_ui <- renderUI({
    req(mapped_seurat_object())
    md <- mapped_seurat_object()@meta.data
    labs <- sort(unique(as.character(md$CTCeek_final_prediction)))
    # metti in cima alcuni gruppi tipici
    priority <- c("CTCmes","CTCepiA","CTCepiB","NK","NK cells","CD4 T","CD8 T","T cells","B cells","Monocytes")
    ordered <- unique(c(priority[priority %in% labs], setdiff(labs, priority)))
    checkboxGroupInput("groups_to_run", NULL, choices = ordered,
                       selected = intersect(ordered, c("CTCmes","CTCepiA","CTCepiB","NK","CD8 T","CD4 T")))
  })
  
  # Marker for selected group (one-vs-all)
  observeEvent(input$find_markers_groups, {
    req(mapped_seurat_object(), input$groups_to_run)
    shinyjs::hide("markers_button_div"); shinyjs::show("markers_spinner_div")
    
    groups <- input$groups_to_run
    selected_groups_rv(groups)
    
    future_promise({
      seurat_obj <- mapped_seurat_object()
      if (!"CTCeek_final_prediction" %in% colnames(seurat_obj@meta.data))
        stop("Prediction metadata ('CTCeek_final_prediction') not found. Ensure UMAP projection is complete.")
      
      Idents(seurat_obj) <- "CTCeek_final_prediction"
      DefaultAssay(seurat_obj) <- "SCT"
      
      out_list <- list()
      for (g in groups) {
        if (!(g %in% levels(Idents(seurat_obj)))) next
        mk <- tryCatch(
          FindMarkers(seurat_obj, ident.1 = g, assay = "SCT", only.pos = TRUE, verbose = FALSE),
          error = function(e) { NULL }
        )
        if (!is.null(mk) && nrow(mk) > 0) {
          mk$gene <- rownames(mk)
          mk <- mk[, c("gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")]
          out_list[[g]] <- mk
        } else {
          out_list[[g]] <- data.frame(
            gene = character(0),
            p_val = numeric(0),
            avg_log2FC = numeric(0),
            pct.1 = numeric(0),
            pct.2 = numeric(0),
            p_val_adj = numeric(0)
          )
        }
      }
      out_list
    }, seed = TRUE) %...>% (function(res_list) {
      per_group_markers_rv(res_list)
      shinyjs::hide("markers_spinner_div"); shinyjs::show("markers_button_div")
      
      # Tabs for each group
      output$per_group_marker_tabs <- renderUI({
        req(per_group_markers_rv())
        tabs <- lapply(names(per_group_markers_rv()), function(g) {
          tabPanel(
            title = g,
            fluidRow(
              column(12,
                     DTOutput(outputId = paste0("tbl_", g)),
                     downloadButton(outputId = paste0("dl_", g), label = paste0("Download ", g, " markers"))
              )
            )
          )
        })
        do.call(tabsetPanel, c(id = "marker_tabs", tabs))
      })
      
      # Table render + download 
      lapply(names(res_list), function(g) {
        local({
          grp <- g
          df  <- res_list[[grp]]
          output[[paste0("tbl_", grp)]] <- renderDT({
            datatable(df, rownames = FALSE, selection = 'single',
                      options = list(pageLength = 10, scrollX = TRUE))
          })
          output[[paste0("dl_", grp)]] <- downloadHandler(
            filename = function() paste0("Markers_", gsub("[^A-Za-z0-9_]+","_", grp), ".csv"),
            content  = function(file) write.csv(df, file, row.names = FALSE)
          )
        })
      })
      
      # global download ZIP 
      output$markers_zip_download_ui <- renderUI({
        req(per_group_markers_rv())
        div(style="text-align:center; margin-top:10px;",
            downloadButton("download_all_markers_zip", "Download ALL marker tables (ZIP)"))
      })
      
    }) %...!% (function(e) {
      shinyjs::hide("markers_spinner_div"); shinyjs::show("markers_button_div")
      showNotification(paste("Marker Calculation Error:", e$message), type="error", duration=NULL)
    })
  })
  
  # ZIP file for markers
  output$download_all_markers_zip <- downloadHandler(
    filename = function() { "CTCeek_AllMarkers.zip" },
    content = function(file) {
      req(per_group_markers_rv())
      td <- tempdir()
      fpaths <- c()
      for (nm in names(per_group_markers_rv())) {
        fn <- file.path(td, paste0("Markers_", gsub("[^A-Za-z0-9_]+","_", nm), ".csv"))
        write.csv(per_group_markers_rv()[[nm]], fn, row.names = FALSE)
        fpaths <- c(fpaths, fn)
      }
      oldwd <- getwd(); setwd(td)
      utils::zip(zipfile = file, files = basename(fpaths))
      setwd(oldwd)
    }
  )
  
}

# --- 4. Run App ---
shinyApp(ui = ui, server = server)

