#!/usr/local/bin/R

suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE))

args = commandArgs(trailingOnly=TRUE)

file.list <- file.info(list.files(getwd(), full.names = T))
latest.file <- rownames(file.list)[which.max(file.list$mtime)]

mat <- as.matrix(data.table::fread(latest.file, header=FALSE))

if(ncol(mat)==8 && nrow(mat)==8){
  colnames(mat) <- paste0("x",seq(1,ncol(mat)))
  rownames(mat) <- paste0("x",seq(1,nrow(mat)))
  mat[lower.tri(mat)] <- NA

  df <- mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("dim1") %>%
  tidyr::pivot_longer(-c(dim1), names_to = "dim2", values_to = "zscore")

  df$dim1 <- factor(df$dim1, levels=rev(sort(unique(df$dim1))))
  diag <- df$dim1 == df$dim2

  mat.plot <- ggplot(df, aes(x=dim2, y=dim1, fill=zscore)) +
    geom_tile(color = "black", size = 0.5) +
    geom_tile(data = df[is.na(df$zscore), ], aes(x=dim2, y=dim1, fill=zscore), color = "gray85", size = 0.5) +
    geom_tile(data = df[diag==TRUE,], aes(x=dim2, y=dim1, fill=zscore), color = "black", size = 0.5) +
    geom_text(data = df %>% dplyr::filter(zscore != 0),
              aes(label=round(zscore, 3)),
              size = 4) +
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = 0,
                         na.value = "gray85",
                         name = "z-score",
                         limits = c(floor(min(df$zscore)), ceiling(max(df$zscore))),
                         guide = guide_colorbar(frame.colour = "black",
                                                frame.linewidth = 1,
                                                ticks.colour = "black",
                                                ticks.linewidth = 1,
                                                barwidth = 1,
                                                barheight = 21.8,
                                                title.position = "top")) +
    theme_bw() +
    theme(panel.border = element_rect(colour="black", size=1),
          axis.title = element_blank(),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12)) +
    scale_x_discrete(expand = c(0, 0),
                     labels = c("\u03bc","h","+","-","\u03c0","A","P","G")) +
    scale_y_discrete(expand = c(0, 0),
                     labels = rev(c("\u03bc","h","+","-","\u03c0","A","P","G")))

  ggsave(filename = paste0(args[1],"_zscore_matrix_plot.tiff"),
         plot = mat.plot,
         device = "tiff",
         width = 6,
         height = 5,
         units = "in")
} else{
  mat <- as.matrix(data.table::fread(latest.file, header=TRUE))
  rownames(mat) <- seq(1,nrow(mat))

  if(any(stringr::str_length(unique(colnames(mat))) > 7)){
    warning("Trucating sequence names. The length of each sequence name should be <= 5 characters.")
    colnames(mat) <- stringr::str_sub(colnames(mat), 1, 7)
    }

  df <- mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("features") %>%
    dplyr::mutate(features = as.numeric(features)) %>%
    tidyr::pivot_longer(-c(features), names_to = "sequences", values_to = "zscore") %>%
    dplyr::arrange(sequences)

    df$sequences <- factor(df$sequences, levels=rev(colnames(mat)))
    df$features <- factor(df$features, levels=sort(unique(df$features)))

    aa.labs <- c(paste0("\u03bc","\u03bc"),
                 paste0("\u03bc","h"),
                 paste0("\u03bc","+"),
                 paste0("\u03bc","-"),
                 paste0("\u03bc","\u03c0"),
                 paste0("\u03bc","A"),
                 paste0("\u03bc","P"),
                 paste0("\u03bc","G"),
                 paste0("h","h"),
                 paste0("h","+"),
                 paste0("h","-"),
                 paste0("h","\u03c0"),
                 paste0("h","A"),
                 paste0("h","P"),
                 paste0("h","G"),
                 paste0("+","+"),
                 paste0("+","-"),
                 paste0("+","\u03c0"),
                 paste0("+","A"),
                 paste0("+","P"),
                 paste0("+","G"),
                 paste0("-","-"),
                 paste0("-","\u03c0"),
                 paste0("-","A"),
                 paste0("-","P"),
                 paste0("-","G"),
                 paste0("\u03c0","\u03c0"),
                 paste0("\u03c0","A"),
                 paste0("\u03c0","P"),
                 paste0("\u03c0","G"),
                 paste0("A","A"),
                 paste0("A","P"),
                 paste0("A","G"),
                 paste0("P","P"),
                 paste0("P","G"),
                 paste0("G","G"))

    if(nrow(mat) >= ncol(mat)){
      mat.plot <- ggplot(df, aes(x=features, y=sequences, fill=zscore)) +
        geom_tile(color = "black", size = 0.5) +
        geom_text(data = df %>% dplyr::filter(zscore != 0),
                  aes(label=round(zscore, 3)),
                  size = 3) +
        scale_fill_gradient2(low = "blue",
                             mid = "white",
                             high = "red",
                             midpoint = 0,
                             na.value = "gray85",
                             name = "z-score",
                             limits = c(floor(min(df$zscore)), ceiling(max(df$zscore))),
                             guide = guide_colorbar(frame.colour = "black",
                                                    frame.linewidth = 1,
                                                    ticks.colour = "black",
                                                    ticks.linewidth = 1,
                                                    barwidth = 25,
                                                    barheight = 1,
                                                    title.position = "left",
                                                    title.vjust = 0.95,
                                                    )) +
        theme_bw() +
        theme(panel.border = element_rect(colour="black", size=1),
              axis.title = element_blank(),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              legend.position = "bottom",
              legend.direction = "horizontal") +
        scale_x_discrete(expand = c(0, 0),
                         position="top",
                         labels = aa.labs) +
        scale_y_discrete(expand = c(0, 0))

        ggsave(filename = paste0(args[1],"_zscore_matrix_plot.tiff"),
               plot = mat.plot,
               device = "tiff",
               width = 15,
               height = (0.5*ncol(mat))+0.9,
               units = "in")
    }
    if(nrow(mat) < ncol(mat)){
      warning("Plotting huge sequence datasets is not yet implemented. Ignoring the -plot flag. Custom plots can be made with the written NARDINI output.")
    }
}
