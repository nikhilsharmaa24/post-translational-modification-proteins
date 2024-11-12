rm(list = ls())
library(tidyverse)

read_n_skip_fun <- function(folder_address_of_csv_files) {
  main_path <<- folder_address_of_csv_files
  
  # Create the folder for ratio
  ratio_folder <<- file.path(main_path, "ratio")
  dir.create(ratio_folder, showWarnings = FALSE)
  
  # reading files
  file_list <- list.files(main_path, pattern = "*.csv", full.names = TRUE)
  
  skip_lines <- sapply(file_list, function(x) grep("prot_hit_num", readr::read_lines(x))[1] - 1, USE.NAMES = FALSE)
  
  list(csv_location = file_list, no_lines_to_skip = skip_lines)
}

get_non_zero_positions <- function(column) {
  positions <- sapply(strsplit(column, ""), function(x) which(as.numeric(x) != 0))
  positions
}

add_peptide_positions <- function(data, prot_seq_col, pep_seq_col, pep_var_mod_pos_col) {
  positions <- mapply(gregexpr, data[[pep_seq_col]], data[[prot_seq_col]])
  
  start_of_pep_in_prot <- sapply(positions, function(x) x[1])
  end_of_pep_in_prot <- sapply(positions, function(x) x[1] + attr(x, "match.length") - 1)
  
  data$start_of_pep_in_prot <- start_of_pep_in_prot
  data$end_of_pep_in_prot <- end_of_pep_in_prot
  
  non_zero_positions <- get_non_zero_positions(data[[pep_var_mod_pos_col]])
  data$position_modified <- Map(function(start, non_zero) {
    if (length(non_zero) > 0) {
      start + non_zero - 1  # Adjust for 0-based indexing
    } else {
      start
    }
  }, data$start_of_pep_in_prot, non_zero_positions)
  
  non_zero_values <- mapply(function(pos, mod_pos) {
    non_zero_vals <- strsplit(gsub("0", "", mod_pos), "")[[1]]
    if (length(non_zero_vals) > 0) {
      non_zero_vals <- gsub("1", "R", non_zero_vals)
      non_zero_vals <- gsub("2", "C", non_zero_vals)
      non_zero_vals <- gsub("3", "K", non_zero_vals)
      non_zero_vals <- gsub("4", "M", non_zero_vals)
      non_zero_vals <- gsub("5", "P", non_zero_vals)
      non_zero_vals <- gsub("6", "T", non_zero_vals)
      non_zero_vals
    } else {
      NA
    }
  }, data$position_modified, data[[pep_var_mod_pos_col]])
  
  data$non_zero_values <- non_zero_values
  
  return(data)
}

calculate_modification_summary <- function(data) {
  count_R <- count_C <- count_K <- count_M <- count_P <- count_T <- numeric(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    non_zero_vals <- data$non_zero_values[[i]]
    if (!is.null(non_zero_vals)) {
      counts <- table(unlist(non_zero_vals))
      count_R[i] <- counts["R"]
      count_C[i] <- counts["C"]
      count_K[i] <- counts["K"]
      count_M[i] <- counts["M"]
      count_P[i] <- counts["P"]
      count_T[i] <- counts["T"]
    }
  }
  
  total_count_R <- sum(count_R, na.rm = T)
  total_count_C <- sum(count_C, na.rm = T)
  total_count_K <- sum(count_K, na.rm = T)
  total_count_M <- sum(count_M, na.rm = T)
  total_count_P <- sum(count_P, na.rm = T)
  total_count_T <- sum(count_T, na.rm = T)
  
  result <- data.frame(total_modified_R = total_count_R,
                       total_modified_C = total_count_C,
                       total_modified_K = total_count_K,
                       total_modified_M = total_count_M,
                       total_modified_P = total_count_P,
                       total_modified_T = total_count_T)
  
  return(result)
}

calculate_summary <- function(data) {
  protein_seq <- data$prot_seq[1]
  total_count <- table(unlist(strsplit(protein_seq, "")))
  
  result <- data.frame(total_count_R = total_count["R"],
                       total_count_C = total_count["C"],
                       total_count_K = total_count["K"],
                       total_count_M = total_count["M"],
                       total_count_P = total_count["P"],
                       total_count_T = total_count["T"]
  )
  
  return(result)
}

raw_to_ratio <- function(Folder_address) {
  list1 <- read_n_skip_fun(Folder_address)
  
  for (i in seq_along(list1$csv_location)) {
    data_list <- list()  # Create an empty list
    data <- read_csv(list1$csv_location[i], skip = list1$no_lines_to_skip[i], col_types = cols())
    # Extract the name from the file path
    file_name <- tools::file_path_sans_ext(basename(list1$csv_location[i]))
    # Assign the data to the list using the extracted name
    data_list[[file_name]] <- data
    
    processed_data <- lapply(data_list, function(dat) {
      dat %>%
        select(prot_hit_num, prot_acc, pep_seq, prot_seq, pep_var_mod, pep_var_mod_pos) %>%
        na.omit() %>%
        distinct() %>%
        mutate(
          end_of_pep_in_prot = str_locate(pattern = pep_seq, prot_seq)[, 2],
          pep_var_mod_pos = sub('.*\\.(.*?)\\..*', '\\1', pep_var_mod_pos),
          pep_var_mod_pos = gsub('2', 0, pep_var_mod_pos),
          pep_var_mod_pos = gsub('4', 0, pep_var_mod_pos)
        ) %>%
        filter(!is.na(end_of_pep_in_prot) & !grepl("^0*$", pep_var_mod_pos))
    })
    
    processed_data <- lapply(processed_data, function(dat) {
      dat %>%
        mutate(non_zero_positions = get_non_zero_positions(pep_var_mod_pos))
    })
    
    
    processed_data <- lapply(processed_data, function(dat) {
      add_peptide_positions(dat, "prot_seq", "pep_seq", "pep_var_mod_pos")
    })
    
    split_processed_data <- lapply(processed_data, function(dat) {
      # Process the current dataset
      processed_dat <- add_peptide_positions(dat, "prot_seq", "pep_seq", "pep_var_mod_pos")
      
      # Split the processed dataset into sublists based on prot_acc
      split_dat <- split(processed_dat, processed_dat$prot_acc)
      
      return(split_dat)
    })
    
    # Apply calculate_modification_summary function on each sublist
    summary_list <- lapply(split_processed_data, function(sublist) {
      lapply(sublist, calculate_modification_summary)
    })
    
    # Apply calculate_summary function on each sublist
    summary_list <- lapply(split_processed_data, function(sublist) {
      lapply(sublist, function(data) {
        modification_summary <- calculate_modification_summary(data)
        overall_summary <- calculate_summary(data)
        
        bind_cols(modification_summary, overall_summary)
      })
    })
    
    lapply(summary_list, function(sublist) {
      do.call(rbind, sublist)
    }) -> final_result
    
    lapply(final_result, function(data) {
      data$ratio_R <- data$total_modified_R / data$total_count_R
      data$ratio_C <- data$total_modified_C / data$total_count_C
      data$ratio_K <- data$total_modified_K / data$total_count_K
      data$ratio_M <- data$total_modified_M / data$total_count_M
      data$ratio_P <- data$total_modified_P / data$total_count_P
      data$ratio_T <- data$total_modified_T / data$total_count_T
      return(data)
    }) -> final_result
    
    # Iterate over each dataframe in final_result and save as CSV
    lapply(names(final_result), function(name) {
      # Generate the filename using the list name and the ratio folder
      filename <- file.path(ratio_folder, paste0("result_", name, ".csv"))
      
      # Save the dataframe as CSV
      write.csv(final_result[[name]], file = filename, row.names = TRUE)
    })
    rm(processed_data, split_processed_data, data_list, data)
  }
}

raw_to_ratio('./demo data/')


data_list <- list()  # Create an empty list
data <- read_csv(list1$csv_location[1], skip = list1$no_lines_to_skip[1], col_types = cols())
# Extract the name from the file path
file_name <- tools::file_path_sans_ext(basename(list1$csv_location[1]))
# Assign the data to the list using the extracted name
data_list[[file_name]] <- data

processed_data <- lapply(data_list, function(dat) {
  dat %>%
    select(prot_hit_num, prot_acc, pep_seq, prot_seq, pep_var_mod, pep_var_mod_pos) %>%
    na.omit() %>%
    distinct() %>%
    mutate(
      end_of_pep_in_prot = str_locate(pattern = pep_seq, prot_seq)[, 2],
      pep_var_mod_pos = sub('.*\\.(.*?)\\..*', '\\1', pep_var_mod_pos),
      pep_var_mod_pos = gsub('2', 0, pep_var_mod_pos),
      pep_var_mod_pos = gsub('4', 0, pep_var_mod_pos)
    ) %>%
    filter(!is.na(end_of_pep_in_prot) & !grepl("^0*$", pep_var_mod_pos))
})

processed_data <- lapply(processed_data, function(dat) {
  dat %>%
    mutate(non_zero_positions = get_non_zero_positions(pep_var_mod_pos))
})


processed_data <- lapply(processed_data, function(dat) {
  add_peptide_positions(dat, "prot_seq", "pep_seq", "pep_var_mod_pos")
})

split_processed_data <- lapply(processed_data, function(dat) {
  # Process the current dataset
  processed_dat <- add_peptide_positions(dat, "prot_seq", "pep_seq", "pep_var_mod_pos")
  
  # Split the processed dataset into sublists based on prot_acc
  split_dat <- split(processed_dat, processed_dat$prot_acc)
  
  return(split_dat)
})

# Apply calculate_modification_summary function on each sublist
summary_list <- lapply(split_processed_data, function(sublist) {
  lapply(sublist, calculate_modification_summary)
})

# Apply calculate_summary function on each sublist
summary_list <- lapply(split_processed_data, function(sublist) {
  lapply(sublist, function(data) {
    modification_summary <- calculate_modification_summary(data)
    overall_summary <- calculate_summary(data)
    
    bind_cols(modification_summary, overall_summary)
  })
})

lapply(summary_list, function(sublist) {
  do.call(rbind, sublist)
}) -> final_result

lapply(final_result, function(data) {
  data$ratio_R <- data$total_modified_R / data$total_count_R
  data$ratio_C <- data$total_modified_C / data$total_count_C
  data$ratio_K <- data$total_modified_K / data$total_count_K
  data$ratio_M <- data$total_modified_M / data$total_count_M
  data$ratio_P <- data$total_modified_P / data$total_count_P
  data$ratio_T <- data$total_modified_T / data$total_count_T
  return(data)
}) -> final_result

valid_column_names <- names(final_result[[1]])[nchar(names(final_result[[1]])) > 0]
# Sort valid column names based on the last letter
sorted_column_names <- sort(valid_column_names, key = function(x) substr(x, nchar(x), nchar(x)))


# Set the output directory path
output_dir <- "D:/Nikhil/demo data for  package/"

# Iterate over each dataframe in final_result and save as CSV
lapply(names(final_result), function(name) {
  # Generate the filename using the list name
  filename <- paste0(output_dir, "result_", name, ".csv")
  print(filename)
  # Save the dataframe as CSV
  write.csv(final_result[[name]], file = filename, row.names = T)
})
