# if (!requireNamespace("types", quietly = TRUE)) {
#     install.packages("types", repos = "https://cloud.r-project.org/")
# }

# if (!requireNamespace("UpSetR", quietly = TRUE)) {
#     install.packages("UpSetR", repos = "https://cloud.r-project.org/")
# }

# if (!requireNamespace("lemon", quietly = TRUE)) {
#     install.packages("lemon", repos = "https://cloud.r-project.org/")
# }

# if (!requireNamespace("dplyr", quietly = TRUE)) {
#     install.packages("dplyr", repos = "https://cloud.r-project.org/")
# }

# if (!requireNamespace("cmapR", quietly = TRUE)) {
#     if (!require("BiocManager", quietly = TRUE)) {
#         install.packages("BiocManager", repos = "https://cloud.r-project.org/")
#     }

#     BiocManager::install("cmapR")
# }
# if (!requireNamespace("data.table", quietly = TRUE)) {
#     install.packages("data.table", repos = "https://cloud.r-project.org/")
# }
# if (!requireNamespace("kableExtra", quietly = TRUE)) {
#     install.packages("kableExtra", repos = "https://cloud.r-project.org/")
#     # tinytex::reinstall_tinytex(repository = "illinois")
#     # tinytex::tlmgr_install('multirow')
# }


library(glue)
library(stringr)
library(dplyr)
library(UpSetR)
library(knitr)
library(lemon)
library(cmapR)
library(data.table)
library(kableExtra)

get_random_dataframes <- function(unique_region_type = `?`(character), path){
      # get directory list of subdirectories of the current
      # working directory
      directories <- list.dirs(path = path, full.names = TRUE, recursive = FALSE)
      
      randomfiles <- list()
      background <- list()
      randomSetnames <- c()
      for (i in directories) {
        if (!(identical(list.files(i, pattern = paste0(
          "*", unique_region_type,
          "_(10|[1-9]|1[1-9]|20)_random_intersect_counts.bed"
        )),
        character(0)))) {
          randomfiles <- c(randomfiles, paste0(i, "/", list.files(i,
                          pattern = paste0(
                            "*", unique_region_type,
                            "_(10|[1-9]|1[1-9]|20)_random_intersect_counts.bed"
                          )
          )))
          randomSetnames <- c(randomSetnames, list.files(i,
                               pattern = paste0("*", 
                                                unique_region_type, 
                                                "_(10|[1-9]|1[1-9]|20)_random_intersect_counts.bed")
          ))
        }
      }
      
      randomnames <- stringr::str_replace(randomSetnames,
                                          pattern = "_random_intersect_counts.bed",
                                          replacement = ""
      )
      samples <- stringr::str_replace(randomnames,
                                      pattern = paste0("_", unique_region_type, "_\\d+"),
                                      replacement = ""
      )
      
      # Read in the random dataframes
      randomdataframes <- lapply(randomfiles,
                                 fread,
                                 header = FALSE,
                                 sep = "\t"
      )
      
      # give each random frame the corresponding sample as a name
      names(randomdataframes) <- samples
      
      return(randomdataframes)
  
}


calculate_background_rank_counts <- function(unique_region_type = `?`(character), path) {
    # read in random dataframes
    randomdataframes <- get_random_dataframes(unique_region_type, path)
    samples <- names(randomdataframes)
    
    
    # get the relative count per random region and rank accordingly
    ranked_random_frames <- list()
    for (i in seq_along(randomdataframes)) {
      randomframe <- randomdataframes[[i]]

      colnames(randomframe) <- c(
        "chr_background", "start",
        "stop", "name", "phase", "strand", "chr_ribo", "start_ribo",
        "stop_ribo", "name_ribo", "phase_ribo", "strand_ribo",
        "nr_bp_overlap"
      )
      
      randomframe$new_name <- paste(randomframe$name,
                                    randomframe$chr_background,
                                    randomframe$start, randomframe$stop,
                                    sep = "_"
      )
      
      result <- randomframe %>%
        mutate(len = stop - start) %>%
        filter(len > 0) %>% 
        group_by(new_name, name_ribo) %>%
        summarize(
          total_bp_overlap = sum(nr_bp_overlap, na.rm = TRUE),
          count_combinations = n(),
          len = first(len),
          chr_background = first(chr_background),
          start = first(start),
          stop = first(stop)
        ) %>%
        # filter(total_bp_overlap > 9) %>%  # Filter based on total bp overlap
        # this is incorrect, as otherwise all 0 counts are disregarded
        group_by(new_name) %>%
        summarize(
          distinct_ribo_count = n_distinct(
            name_ribo[total_bp_overlap > 9]
          ),
          len = first(len),
          chr_background = first(chr_background),
          start = first(start),
          stop = first(stop)
        ) %>%
        mutate(
          num_reads = as.numeric(distinct_ribo_count),
          len = as.numeric(len),
          relative_count = distinct_ribo_count / len
        ) %>%
        arrange(desc(relative_count))
      
      
      if (names(randomdataframes)[[i]] %in% names(ranked_random_frames)){
        ranked_random_frames[[names(randomdataframes)[[i]]]] <- c(ranked_random_frames[[names(randomdataframes)[[i]]]], result)
      }else{
        ranked_random_frames[[names(randomdataframes)[[i]]]] <- c(result)
      }
    }
    
    # loop through the samples and average the random counts per rank
    random_ranked_dfs <- list()
    for (sample in unique(samples)){
        relative_count_background_list <- ranked_random_frames[[sample]][names(ranked_random_frames[[sample]]) %in% "relative_count"]
        background_rel_count_matrix <- do.call(cbind, relative_count_background_list)
        mean_relative_background_counts <- rowMeans(background_rel_count_matrix)
        random_region_names <- paste0("RR", 1:length(mean_relative_background_counts))
        random_region_type <- rep("RR", length(mean_relative_background_counts))
        
        
        random_ranked_df <- data.frame(
          new_name = random_region_names,
          name = random_region_names,
          relative_count = mean_relative_background_counts,
          region_type = random_region_type,
          distinct_ribo_count = rep(NA, length(mean_relative_background_counts)),
          len = rep(NA, length(mean_relative_background_counts)),
          chr_unique = rep(NA, length(mean_relative_background_counts)),
          start = rep(NA, length(mean_relative_background_counts)),
          stop = rep(NA, length(mean_relative_background_counts)),
          genomic_region = rep(NA, length(mean_relative_background_counts)),
          num_reads = rep(NA, length(mean_relative_background_counts))
        )
        
        random_ranked_dfs[[sample]] <- random_ranked_df
    }
    return(random_ranked_dfs)
}





get_intersect_files <- function(unique_region_type = `?`(character), path) {
    directories <- list.dirs(path = path, full.names = TRUE, recursive = TRUE)
    files <- list() # Initialize an empty list
    basenames <- list() # Initialize an empty list
    for (i in directories) {
        # Use list.files with the correct pattern to find
        # matching files
        matching_files <- list.files(i, pattern = paste0(
            "*",
            unique_region_type, "_intersect_counts_sorted.bed"
        ), full.names = TRUE)

        # If matching files are found, append to the lists
        if (length(matching_files) > 0) {
            files <- append(files, matching_files)
            basenames <- append(basenames, basename(matching_files))
            # Use basename to get file names only
        }
    }


    # Read the files into a list of data frames
    dataframes <- lapply(files, fread, header = FALSE, sep = "\t")

    # Assign names to the data frames based on the basenames
    names(dataframes) <- stringr::str_replace(basenames,
        pattern = paste0("_", unique_region_type, "_intersect_counts_sorted.bed"),
        replacement = ""
    )

    return(dataframes)
}


get_relative_ur_intersect_counts <- function(dataframes) {
    frames_preprocessed <- list()
    genomic_frames <- list()
    unique_names_per_sample <- list()

    for (frame in dataframes) {
        colnames(frame) <- c(
            "chr_unique", "start", "stop", "name",
            "phase", "strand", "chr_ribo", "start_ribo", "stop_ribo",
            "name_ribo", "phase_ribo", "strand_ribo", "nr_bp_overlap"
        )


        frame$new_name <- paste(frame$name, frame$chr_unique,
            frame$start, frame$stop,
            sep = "_"
        )

        UR_lengths <- frame %>%
            mutate(len = stop - start) %>%
            group_by(name) %>%
            distinct(start, len) %>%
            summarize(len = sum(len))

        frame <- frame %>% left_join(UR_lengths, by = "name")

        processed_frame <- frame %>%
            filter(nr_bp_overlap > 0) %>%
            group_by(name, name_ribo) %>%
            summarize(
                total_bp_overlap = sum(nr_bp_overlap, na.rm = TRUE),
                count_combinations = n(),
                len = first(len),
                chr_unique = first(chr_unique),
                start = min(start),
                stop = max(stop),
                new_name = paste(new_name, collapse = "_")
            )

        processed_frame <- processed_frame %>%
            filter(total_bp_overlap > 9) %>% # Filter based on total bp overlap
            group_by(name) %>% # Group by name to count distinct name_ribo per name
            summarize(
                distinct_ribo_count = n_distinct(name_ribo),
                len = first(len),
                chr_unique = first(chr_unique),
                start = min(start),
                stop = max(stop),
                new_name = paste(strsplit(first(new_name), "_")[[1]][1], start, stop, sep = "_"),
                genomic_region = paste(chr_unique, start, stop, sep = "_"),
                region_type = "UR"
            ) %>%
            mutate(
                num_reads = as.numeric(distinct_ribo_count),
                len = as.numeric(len),
                relative_count = distinct_ribo_count / len
            ) %>%
            arrange(desc(relative_count))

        genomic_frame <- processed_frame %>%
            # do not count the same genomic region twice
            # group_by genomic region and just take the first ORF information
            group_by(genomic_region) %>%
            summarize(
                distinct_ribo_count = first(distinct_ribo_count),
                len = first(len),
                chr_unique = first(chr_unique),
                start = first(start),
                stop = first(stop),
                new_name = first(new_name),
                name = first(name),
                num_reads = first(num_reads),
                relative_count = first(relative_count),
                region_type = first(region_type)
            ) %>%
            arrange(desc(relative_count))

        # drop the genomic region as anyway present in the new_name column
        genomic_frame <- genomic_frame %>% select(-1)
        colnames(genomic_frame) <- c(
            "distinct_ribo_count",
            "len",
            "chr_unique",
            "start",
            "stop",
            "new_name",
            "name",
            "num_reads",
            "relative_count",
            "region_type"
        )

        frames_preprocessed[[length(frames_preprocessed) +
            1]] <- processed_frame
        genomic_frames[[length(genomic_frames) +
            1]] <- genomic_frame
    }
    
    
    intersect_frames_list <- list()
    intersect_frames_list$frames_processed <- frames_preprocessed
    intersect_frames_list$genomic_frames <- genomic_frames
    
    return(intersect_frames_list)
    }
    
    
    
    
    
    
calculate_fdr_urs <- function(intersect_frames_list, random_ranked_dfs){
    ur_rr_concat_dfs <- list()
    ur_dfs <- list()
    significant_dfs <- list()
    unique_regions <- list()
    upsetlist <- list()
    relevantregionscount <- list()
    genomic_ur_dfs <- list()
    genomic_unique_regions <- list()
    relevantregionscount_genomic <- list()
    upsetlist_genomic <- list()
    
    for (i in seq_along(intersect_frames_list$frames_processed)){
      random_df_colnames_order <- random_ranked_dfs[[i]][colnames(intersect_frames_list$frames_processed[[i]])]
      nr_random_regions <- length(rownames(random_df_colnames_order))
      ur_rr_concat_df <- rbind(intersect_frames_list$frames_processed[[i]], random_df_colnames_order)
      ur_rr_concat_df <- ur_rr_concat_df %>%
        arrange(desc(relative_count))
      
      
      ur_rr_concat_df <- ur_rr_concat_df %>% mutate(cumsum_rr = cumsum(region_type == "RR"))
      
      ur_rr_concat_df$p_value <- (ur_rr_concat_df$cumsum_rr + 1)/(nr_random_regions + 1)
      
      ur_df <- ur_rr_concat_df[ur_rr_concat_df$region_type == "UR",]
      
      # nr_random_regions is the number of random regions originally sampled = nr of unique regions considered
      ur_df$q_value <- p.adjust(ur_df$p_value, method = "BH", n = nr_random_regions)
      
      # for significance also require at least 2 reads
      ur_df$significant <- ifelse(ur_df$q_value < 0.05 & ur_df$num_reads > 2, 1, 0)
      
      ur_rr_concat_dfs[[i]] <- ur_rr_concat_df
      ur_dfs[[i]] <- ur_df
      
      significant_df <- ur_df[ur_df$significant == 1,]
      significant_dfs[[i]] <- significant_df
      
      
      relevantregionscount <- c(relevantregionscount, nrow(significant_df))
      unique_regions <- significant_df$name
      upsetlist <- c(upsetlist, list(unique_regions))
      
      
      ##########################################################################
      # Genomic frame handling
      genomic_frame <- intersect_frames_list$genomic_frames[[i]]
      genomic_frame <- merge(genomic_frame, ur_df[, c("name", "q_value", "p_value", "significant")],
                         by = "name", 
                         all.x = TRUE) 
      genomic_frame$ID <- genomic_frame$name
      genomic_frame_significant <- genomic_frame[genomic_frame$significant == 1,]
      relevantregionscount_genomic <- c(relevantregionscount_genomic, nrow(genomic_frame_significant))
      genomic_unique_regions <- genomic_frame_significant$name
      upsetlist_genomic <- c(upsetlist_genomic, list(genomic_unique_regions))
      genomic_ur_dfs[[i]] <- genomic_frame
      
    }
    
    counts_and_names_list <- list()
    counts_and_names_list$counts_above_t_relative_length <- relevantregionscount
    counts_and_names_list$unique_names <- upsetlist
    counts_and_names_list$ur_rr_concat_dfs <- ur_rr_concat_dfs
    counts_and_names_list$ur_dfs <- ur_dfs
    counts_and_names_list$genomic_ur_dfs <- genomic_ur_dfs
    names(counts_and_names_list$genomic_ur_dfs) <- names(counts_and_names_list$ur_dfs)
    
    
    counts_and_names_list$counts_above_t_relative_length_genomic <- relevantregionscount_genomic
    counts_and_names_list$unique_names_genomic <- upsetlist_genomic
    return(counts_and_names_list)
}  
    
    
    


print_unique_region_above_t <- function(
    relevantregionscount,
    dataframes) {
    printframe <- data.frame(x = relevantregionscount)
    printframe <- as.data.frame(t(printframe))
    rownames(printframe) <- names(dataframes)
    colnames(printframe) <- paste0("Number of unique regions with relative count >= threshold")
    kable(printframe,
        caption = "Regions above the threshold",
        escape = TRUE
    ) %>%
        kable_styling(font_size = 8)
}


get_top_5_unique_regions <- function(dataframes_list_genomic) {
  # upsetlist <- list()
  j <- 1
  # names <- names(dataframes_list_genomic)
  for (unique_region_frame in dataframes_list_genomic) {
    unique_region_frame <- unique_region_frame %>%
      filter(significant == 1) %>%
      arrange(desc(relative_count))
    
    if (dim(unique_region_frame)[1] > 0) {
      # unique_regions <- unique_region_frame$ID
      # upsetlist <- c(upsetlist, list(unique_regions))
      
      unique_region_frame <- unique_region_frame[, c("ID", "len", "chr_unique", "start", "stop", "num_reads", "relative_count", "q_value")]
      
      unique_region_frame$ID <-
        gsub("\\|", "-", unique_region_frame$ID)
      
      
      print(kable(unique_region_frame[1:5, ],
            caption = names(dataframes_list_genomic)[j],
            row.names = FALSE, escape = TRUE
      ) %>%
        kable_styling(font_size = 8))
    } 
    # else {
      # names <- names[!names %in% c(names[j])]
    # }
    j <- j + 1
  }
  # names(upsetlist) <- names
}



get_start_positions <- function(pos_string) {
    single_ORFs <- unlist(str_split(pos_string, ","))
    single_positions <- unlist(str_split(single_ORFs, "-"))
    ORF_starts <- sort(as.integer(single_positions[seq(1, length(single_positions), by = 2)]))
    return(ORF_starts)
}

get_print_ORF_ranks <- function(SO_pipe_path, frames) {
    SO_results <- read.delim(file.path(SO_pipe_path, "UniqueProteinORFPairs.txt"))
    SO_results$ORF_starts <- sapply(SO_results$OrfPos, get_start_positions)
    head(SO_results)

    nr_URs <- c()
    nr_first <- c()
    nr_second <- c()
    perc_first <- c()
    perc_second <- c()
    return_frames <- list()
    for (frame in frames) {
        # filter for significant URs
        # frame <- frame[frame$significant == 1, ]
        frame$OrfTransID <- str_split_fixed(frame$name, "\\|", 4)[, 2]
        frame$OrfTransID <- str_split_fixed(frame$OrfTransID, ":", 2)[, 1]
        frame$ORF_start <- as.integer(str_split_fixed(frame$name, ":", 4)[, 3])
        merged_df <- merge(frame, SO_results, by = "OrfTransID")
        merged_df$ORF_rank <- mapply(function(x, y) match(x, y), merged_df$ORF_start, merged_df$ORF_starts)
        # print(head(merged_df[c('ORF_start', 'ORF_starts', 'OrfPos', 'ORF_rank')]))

        ORF_ranks_dict <- setNames(merged_df$ORF_rank, merged_df$OrfTransID)

        # Map values from df1 to df2 based on the shared column
        frame$ORF_ranks <- ORF_ranks_dict[as.character(frame$OrfTransID)]

        nr_URs <- c(length(rownames(frame)), nr_URs)
        nr_first <- c(length(rownames(merged_df[merged_df$ORF_rank == 1, ])), nr_first)
        nr_second <- c(length(rownames(merged_df[merged_df$ORF_rank > 1, ])), nr_second)
        perc_first <- c(length(rownames(merged_df[merged_df$ORF_rank == 1, ])) / length(rownames(frame)), perc_first)
        perc_second <- c(length(rownames(merged_df[merged_df$ORF_rank > 1, ])) / length(rownames(frame)), perc_second)
        return_frames[[length(return_frames) + 1]] <- frame
    }

    NR_ORFs_df <- data.frame(nr_URs = nr_URs, nr_first = nr_first, nr_second = nr_second, perc_first = perc_first, perc_second = perc_second)
    rownames(NR_ORFs_df) <- names(counts_and_names_list$frames_processed)
    kable(NR_ORFs_df, caption = "ORF positions of validated URs") %>%
        kable_styling(font_size = 8)
    return(return_frames)
}


write_csv <- function(dataframes_genomic, unique_names_per_sample_genomic, dataframes, unique_names_per_sample, outdir) {
    frame_number <- 1
    for (frame in dataframes_genomic) {
        
        write_frame <- frame[, c(
            "ID",
            "len",
            "chr_unique",
            "start",
            "stop",
            "new_name",
            "num_reads",
            "relative_count",
            "p_value",
            "q_value",
            "significant"
        )]

        write.csv(write_frame, file.path(outdir, paste(names(dataframes_genomic)[frame_number],
            "genomicunique_regions.csv",
            sep = "_"
        )))
        frame_number <- frame_number + 1
    }

    frame_number <- 1
    for (frame in dataframes) {
        # colnames(frame) <- c()
        # frame$significant <- ifelse(frame$name %in% unique_names_per_sample[[frame_number]],
        #     1, 0
        # )

        write_frame <- frame[, c(
            "name",
            "len",
            "chr_unique",
            "start",
            "stop",
            "new_name",
            "num_reads",
            "relative_count",
            "significant",
            "p_value",
            "q_value",
            "ORF_ranks"
        )]

        write.csv(write_frame, file.path(outdir, paste(names(dataframes)[frame_number],
            "unique_regions.csv",
            sep = "_"
        )))
        frame_number <- frame_number + 1
    }
}


Split_ORFs_validation <- function(dataframes, unique_names_per_sample, path) {
    # Assign column names

    frame_number <- 1
    for (unique_region_frame in dataframes) {
        unique_region_frame <- unique_region_frame[, c("name", "distinct_ribo_count", "len", "chr_unique", 
                                                       "start", "stop", "new_name", "num_reads", "p_value",
                                                       "q_value", "relative_count", "significant")]

        # print(unique_region_frame)
        unique_region_frame <- unique_region_frame %>%
            filter(significant == 1)

        if (nrow(unique_region_frame) > 0) {
            unique_region_frame$transcript <- sapply(strsplit(unique_region_frame$name, ":"), `[`, 1)
            unique_region_frame$orf <- sapply(strsplit(unique_region_frame$name, ":"), `[`, 2)

            SO_2_uniques_validated <- unique_region_frame %>%
                group_by(transcript) %>%
                summarise(
                  across(everything(), first),
                  ribo_cov_orfs = n_distinct(orf)) %>%
                filter(ribo_cov_orfs >= 2) %>%
                ungroup() %>%
                arrange(transcript)
            
            if (nrow(SO_2_uniques_validated) > 0) {
              write.csv(SO_2_uniques_validated[, c(
                  "transcript",
                  "orf",
                  "distinct_ribo_count",
                  "len",
                  "chr_unique",
                  "start",
                  "stop",
                  "new_name",
                  "num_reads",
                  "p_value",
                  "q_value",
                  "relative_count"
              )], file.path(path, paste0(names(dataframes)[frame_number], "_two_regions_validated_on_transcript.csv")))
              frame_number <- frame_number + 1
            }
        }
  }
}

        
