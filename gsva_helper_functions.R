gsva_wrapper_classII <- function(output_dir, ranks_input, reference_ranks, metadata, genesets, algorithm, analysis_name, tool_colour){

    annot_colouring <- data.frame(row.names = colnames(ranks_input)) %>%
        mutate(annotation_group = metadata$annotation[match(colnames(ranks_input), metadata[, c("filename", "annotation")]$filename)]) %>% 
        dplyr::arrange(annotation_group)

    annotation_list <- sort(unique(metadata$annotation))

    # initialise final df
    gsva_output_full <- data.frame(GOI_set = sort(names(genesets)))
    invisible(lapply(c('ES_matrix', 'ES_plots'), function(x) dir.create(file.path(output_dir, x), recursive = TRUE)))

    # run algorithm across annotation group, e.g. sample from cancer type against all other cancer types
    for (annotation_group in annotation_list){

        annotation_group_og <- annotation_group
        # simpler/shorter name saving for later (also avoid file saving issue with spaces)
        if (grepl("\\(", annotation_group) & grepl("\\)", annotation_group)) {
            annotation_group <- sub(".*\\((.*)\\).*", "\\1", annotation_group)
        }
        annotation_group <- gsub(' ', '_', gsub('/', '-', gsub(' / ', '-', annotation_group)))

        print(annotation_group)
        gsva_output <- data.frame(GOI_set = names(genesets))

        metadata_of_annotation <- metadata %>% 
            filter(annotation == annotation_group_og)
        ranks_input_filt <- ranks_input[, intersect(colnames(ranks_input), metadata_of_annotation$filename), drop = FALSE]

        print(paste0('-----Running ', toupper(algorithm), '------'))
        # each sample ran individually against the rest
        for (i in 1:length(names(ranks_input_filt))){
            sample <- names(ranks_input_filt)[[i]]
            sample_df <- ranks_input_filt %>% select(sample)
            sample_df <- sample_df[!is.na(sample_df), , drop = F]
            input_df <- merge(sample_df, reference_ranks, by = "row.names", all = TRUE)
            input_df <- column_to_rownames(input_df, var = 'Row.names')
            input_df[] <- lapply(input_df, function(x) as.numeric(as.character(x)))
            
            # format gsva input
            if (toupper(algorithm) == 'GSVA'){
                gsvaparam_object <- gsvaParam(as.matrix(input_df), genesets,
                                                kcdf = 'auto', sparse = FALSE,
                                                minSize = 10, use = 'na.rm')
            } else if (toupper(algorithm) == 'PLAGE'){
                gsvaparam_object <- plageParam(as.matrix(input_df), genesets, minSize = 10)
            } else if (toupper(algorithm) == 'ZSCORE'){
                gsvaparam_object <- zscoreParam(as.matrix(input_df), genesets, minSize = 10)
            } else if (toupper(algorithm) == 'SSGSEA'){
                gsvaparam_object <- ssgseaParam(as.matrix(input_df), genesets, normalize = TRUE, minSize = 10)
            } else {
                message = "ERROR: Invalid algorithm specified, must be one of [GSVA, PLAGE, ZSCORE, SSGSEA]"
                stop(message)
            }
            # calculate ES with gene ranking
            gsva_results <- as.data.frame(gsva(gsvaparam_object, verbose = F))
            # take score only of sample of interest
            sample_output <- as.data.frame(gsva_results)[, 1, drop = F]
            sample_output['GOI_set'] <- row.names(sample_output)
            gsva_output <- merge(gsva_output, sample_output, by = "GOI_set", all = TRUE)
        }
        gsva_output <- column_to_rownames(gsva_output, var = 'GOI_set')
        gsva_output <- gsva_output[order(rownames(gsva_output)), , drop = FALSE] # genesets in alphabetical order
        plot_GSEApheatmap_wNAs(gsva_output,
            png_name = file.path(output_dir, 'ES_plots', paste0(annotation_group, ".png")),
            plot_title = paste0(annotation_group_og, " - gsva_", toupper(algorithm)),
            tool_colour = tool_colour)
        gsva_output <- rownames_to_column(gsva_output, var = 'GOI_set')
        write.table(gsva_output, file.path(output_dir, 'ES_matrix', paste0(annotation_group, ".tsv")), sep ='\t', quote = FALSE, row.names = FALSE)
        gsva_output_full <- merge(gsva_output_full, gsva_output, by = "GOI_set", all = TRUE)
    }
    gsva_output_full <- column_to_rownames(gsva_output_full, var = 'GOI_set')
    gsva_output_full <- gsva_output_full[order(rownames(gsva_output_full)), , drop = FALSE]
    gsva_output_full <- gsva_output_full %>% # if .x suffix at end of sample names
                    rename_with(~ str_replace(.x, "\\.x$", ""), ends_with(".x"))
    write.table(rownames_to_column(gsva_output_full, var = 'GOI_Set'), file.path(output_dir, paste0(analysis_name, "-fullNES.tsv")), sep ='\t', quote = FALSE, row.names = FALSE)
    # no legend
    plot_GSEApheatmap_wNAs(gsva_output_full,
        png_name = file.path(output_dir, 'ES_plots', paste0(analysis_name, ".png")),
        plot_title = paste0(analysis_name, " - gsva_", toupper(algorithm)),
        tool_colour = tool_colour, wannotation = annot_colouring)
    # with legend
    plot_GSEApheatmap_wNAs(gsva_output_full,
        png_name = file.path(output_dir, 'ES_plots', paste0(analysis_name, "_annot.png")),
        plot_title = paste0(analysis_name, " - gsva_", toupper(algorithm)),
        tool_colour = tool_colour, wannotation = annot_colouring, wlegend = TRUE)
}

gsva_wrapper_classI_III <- function(ranks_data, metadata, genesets, analysis_name, output_folder, input_type, tool_colour, algorithm){

    annot_colouring <- data.frame(row.names = colnames(ranks_data)) %>%
        mutate(annotation = metadata$annotation[match(colnames(ranks_data), metadata[, c("filename", "annotation")]$filename)]) %>% 
        dplyr::arrange(annotation)

    ranks_data_ordered <- ranks_data %>% 
        dplyr::select(rownames(annot_colouring)) %>%
        mutate(across(where(is.integer), as.numeric))
    
    invisible(lapply(c('ES_matrix', 'ES_plots'), function(x) dir.create(file.path(output_folder, x), recursive = TRUE)))

    if (toupper(algorithm) == 'GSVA'){
        gsvaparam_object <- gsvaParam(as.matrix(ranks_data_ordered), genesets,
                                        kcdf = 'auto', sparse = FALSE,
                                        minSize = 10, use = 'na.rm')
    } else if (toupper(algorithm) == 'PLAGE'){
        gsvaparam_object <- plageParam(as.matrix(ranks_data_ordered), genesets, minSize = 10)
    } else if (toupper(algorithm) == 'ZSCORE'){
        gsvaparam_object <- zscoreParam(as.matrix(ranks_data_ordered), genesets, minSize = 10)
    } else if (toupper(algorithm) == 'SSGSEA'){
        gsvaparam_object <- ssgseaParam(as.matrix(ranks_data_ordered), genesets, normalize = TRUE, minSize = 10, use = 'na.rm')
    } else {
        message = "ERROR: Invalid algorithm specified, must be one of [GSVA, PLAGE, ZSCORE, SSGSEA]"
        stop(message)
    }
    ssgsea_results <- as.data.frame(gsva(gsvaparam_object, verbose = T))
    ssgsea_output <- ssgsea_results[order(rownames(ssgsea_results)), , drop = FALSE] # genesets in alphabetical order
    write.table(rownames_to_column(ssgsea_output, var = 'GOI_Set'), file.path(output_folder, paste0(analysis_name, "-fullNES.tsv")), sep ='\t', quote = FALSE, row.names = FALSE)

    plot_GSEApheatmap_wNAs(ssgsea_output,
        file.path(output_folder, 'ES_plots', paste0(analysis_name, ".png")),
        paste0(analysis_name, " - ssgsea_", input_type),
        tool_colour, wannotation = annot_colouring)
    plot_GSEApheatmap_wNAs(ssgsea_output,
        file.path(output_folder, 'ES_plots', paste0(analysis_name, "_annot.png")),
        paste0(analysis_name, " - ssgsea_", input_type),
        tool_colour, wannotation = annot_colouring, wlegend = TRUE)

    # iterate throgh annotation groups for subplots/individual ES matrices
    annotation_list <- unique(metadata$annotation)
    for (annotation_group in annotation_list){

        annotation_group_og <- annotation_group
        # simpler/shorter name saving for later (also avoid file saving issue with spaces)
        if (grepl("\\(", annotation_group) & grepl("\\)", annotation_group)) {
            annotation_group <- sub(".*\\((.*)\\).*", "\\1", annotation_group)
        }
        annotation_group <- gsub(' ', '_', gsub('/', '-', gsub(' / ', '-', annotation_group)))
        print(annotation_group)

        metadata_of_annotation <- metadata %>% 
            filter(annotation == annotation_group_og)
        ssgsea_output_annot <- ssgsea_output[, colnames(ssgsea_output) %in% metadata_of_annotation$filename, drop = FALSE]
        
        write.table(rownames_to_column(ssgsea_output_annot, var = 'GOI_Set'), file.path(output_folder, 'ES_matrix', paste0(annotation_group, ".tsv")), sep = '\t', quote = FALSE, row.names = FALSE)
        
        plot_GSEApheatmap_wNAs(ssgsea_output_annot,
            file.path(output_folder, 'ES_plots', paste0(annotation_group, '.png')),
            paste0(annotation_group_og, " - ssgsea_", input_type),
            tool_colour)
    }
}


plot_GSEApheatmap_wNAs <- function(ES_matrix, png_name, plot_title, tool_colour, wannotation = NA, wlegend = FALSE){
    # function to plot all GSEA pheatmaps the same way, handling NAs in ES matrix outputs
    ES_matrix[ES_matrix == "---"] <- NA
    for (i in seq_along(ES_matrix)) {
        ES_matrix[[i]] <- suppressWarnings(as.numeric(ES_matrix[[i]]))
    }
    
    ES_matrix <- as.matrix(ES_matrix)
    
    # --- Scale columns ignoring NAs ---
    scale_columns <- function(mat) {
        out <- mat
        
        for (j in seq_len(ncol(mat))) {
            col <- mat[, j]
            
            # if all values are NA → leave column untouched
            if (all(is.na(col))) next
            
            m <- mean(col, na.rm = TRUE)
            s <- sd(col, na.rm = TRUE)
            
            # If no variation: leave column as-is
            if (is.na(s) || s == 0) {
                next
            } else {
                out[, j] <- (col - m) / s
            }
        }
        out
    }
    
    scaled_ES <- scale_columns(ES_matrix)
    myColor <- colorRampPalette(c("white", "white", tool_colour))(101)

    # --- Create guaranteed-unique breaks ---
    data_min <- min(scaled_ES, na.rm = TRUE)
    data_max <- max(scaled_ES, na.rm = TRUE)


    # If all values equal, expand artificial range
    if (data_min == data_max) {
        data_min <- data_min - 0.1
        data_max <- data_max + 0.1
    }

    breaks <- seq(data_min, data_max, length.out = length(myColor) + 1)

    png(png_name)
    pheatmap::pheatmap(scaled_ES,
                    show_rownames = T, show_colnames = F,
                    treeheight_row = 0, treeheight_col = 0,
                    cluster_cols = F, cluster_rows = F, scale = 'none',
                    color = myColor, breaks = breaks,
                    main = plot_title, fontsize_row = 5,
                    annotation_col = wannotation,
                    annotation_legend = wlegend)
    dev.off()
}

