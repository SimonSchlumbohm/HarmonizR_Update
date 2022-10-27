#' Main function
#'
#' This function executes the entire harmonizR program and executes all other
#' functions found in this package. Therefore, this is the only function in
#' need of calling.
#'
#' @param data_as_input Path to input data. Please read the SOP for the correct
#' input format found on https://github.com/SimonSchlumbohm/HarmonizR.
#' Additionally, the input can be a data.frame with proper row- and column
#' names.
#' @param description_as_input Path to input description. Please read the SOP
#' for the correct input format found on
#' https://github.com/SimonSchlumbohm/HarmonizR. Additionally, the input can be
#' a data.frame with three columns total.
#' @param algorithm Optional. Pass either "ComBat" or "limma" to select the
#' preferred adjustment method. Defaults to ComBat.
#' @param ComBat_mode Optional. Pass a number between 1 and 4 to select the
#' desired ComBat parameters. Can only be set when ComBat is used. For
#' information on the meaning of the numbers, please view the SOP found on
#' https://github.com/SimonSchlumbohm/HarmonizR. Defaults to 1.
#' @param plot Optional. Takes either "samplemeans" for sample specific means,
#' "featuremeans" for feature specific means or "CV" for the coeffivient of
#' variation as input and creates before/after plots for the given data.
#' Defaults to FALSE -> Turned off.
#' @param sort XXX
#' @param block XXX
#' @param output_file Optional. Takes a string as input for the .tsv file name.
#' This can also be a path. Defaults to "cured_data", hence yielding a
#' "cured_data.tsv" file.
#' @param prints Optional. Toggles the amount of stuff printed out by the
#' HarmonizR algorithm during execution.
#' @return A .tsv file by default called cured_data.tsv will be written out as
#' a result.
#' As a return, the harmonizR function will yield the batch effect adjusted
#' data.frame.
#' @export


# ----- MAIN (FUNCTION CALLS) -------------------------------------------------
# Contents: harmonizR()


# ----- harmonizR(
#         data_as_input,
#         description_as_input,
#         algorithm,
#         ComBat_mode,
#         plot,
#         sort,
#         block,
#         output_file,
#         prints
#       ) ---------------------------------------------------------------------
# This function will be called by the User and will execute the entire program
harmonizR <- function(data_as_input = NULL,
                      description_as_input = NULL,
                      ...,
                      algorithm = "ComBat",
                      ComBat_mode = 1,
                      plot = FALSE,
                      sort = FALSE,
                      block = NULL,
                      output_file = "cured_data",
                      prints = "standard") {

  # Check whether data and description are present
  if (is.null(data_as_input) && is.null(description_as_input)) {
    stop("No parameters given. Usage: harmonizR(\"path/to/data\", \"path/to/description\")")
  } else if (is.null(description_as_input)) {
    stop("Not enough parameters. Usage: main(\"path/to/data\", \"path/to/description\")")
  }

  # Check algorithm input
  if (algorithm != "ComBat" && algorithm != "limma") {
    print("Please set the algorithm parameter to either ComBat or limma to choose the prefered adjustment. Parameter is now being set to default (ComBat).")
    algorithm <- "ComBat"
  }

  # Check ComBat_mode input
  if (ComBat_mode < 1 || ComBat_mode > 4) {
    print("Invalid ComBat Mode chosen. Select 1, 2, 3 or 4. Parameter is now being set to default (1).")
    ComBat_mode <- 1
  }

  # Check plot input
  if (plot != FALSE &&
    plot != "samplemeans" &&
    plot != "featuremeans" &&
    plot != "CV") {
    print("Please set the plot parameter to either samplemeans, featuremeans, CV or FALSE.")
    plot <- FALSE
  }
  
  # Check sort input TODO! SORT LOGIC HERE!
  #if (sort != FALSE && sort != TRUE){
  #  print("Please set the sort parameter to either TRUE or FALSE. Parameter is now set to default (FALSE)")
  #  sort <- FALSE
  #}
  #if (sort == TRUE && is.null(block)){
  #  print("Please set the block parameter alongside the sort parameter since sorting is pointless without blocking.")
  #  sort <- FALSE
  #}

  # Check output_file input
  if (is.character(output_file) != TRUE) {
    print("Please only pass a string via the output_file parameter.")
    output_file <- "cured_data"
  }
  
  # Check prints input
  if (prints != "standard" && prints != "mute" && prints != "all") {
    print("Please set the prints parameter to either 'mute', 'standard' or 'all'. Parameter is now being set to default ('standard').")
    prints <- "standard"
  }
  
  # PRINT-OUT
  if (prints == "all") {
    # This line checks the OS
    print(paste0("Current operating system: ", toupper(.Platform$OS.type)))
  }


  # ----- READING SECTION -----------------------------------------------------
  # During reading, R already removes completely empty rows in the data.frame
  # PRINT-OUT
  if (prints == "standard" || prints == "all") {
    print("Reading the files...")
  }

  # Logic for the Perseus plugin
  if (is.character(data_as_input)) {
    # Read in the data
    main_data <- read_main_data(data_as_input)
  } else {
    # Read in the data
    main_data <- data_as_input
  }

  # Logic for the Perseus plugin
  if (is.character(description_as_input)) {
    # Read in the batch-descriptions
    batch_data <- read_description(description_as_input)
  } else {
    # Read in the batch-descriptions
    batch_data <- description_as_input
  }

  # Information about which row of batch_data is in which batch
  batch_list <- fetch_batch_overview(batch_data)
  
  # A backup later to be used for visualize since I never re-sort the
  # description file
  backup_batch_list <- batch_list
  
  number_features_before <- dim(main_data)[1]

  # Remove duplicates
  main_data <- unique(main_data)
  
  number_features_after <- dim(main_data)[1]
  
  if (prints == "all") {
    duplicates <- number_features_before - number_features_after
    print(paste("Duplicates found and removed:", duplicates))
  }
  
  order_to_go_by <- FALSE
  
  
  # ----- SORTING SECTION -----------------------------------------------------
  # This is toggleable by the user
  if (sort == TRUE || sort == "sparsity_sort" || sort == "seriation_sort" || sort == "jaccard_sort") {
    
    # PRINT-OUT
    if (prints == "standard" || prints == "all") {
      print(paste("Sorting with", sort, "..."))
    }
    
    order_to_go_by <- FALSE
    
    
    # - - - BINARY REDUCTION SUBSECTION - - - - - - - - - - - - - - - - - - - -
    # Needed for seriation- and jaccard-sorting
    if (sort == "seriation_sort" || sort == "jaccard_sort"){
      # Here, binary matrix reduction is used to get "binary_df"
      bin_input <- main_data
      bin_input[!is.na(bin_input)] <- 1
      bin_input[is.na(bin_input)] <- 0
      needed_val <- 1
      
      if (ComBat_mode == 1 ||
          ComBat_mode == 3 ||
          algorithm == "limma") {
        needed_val <- 2
      }
    
      # Function call
      binary_df <- binary_matrix_reduction(bin_input, batch_list, needed_val)
      #print(binary_df)
    
      if (sort == "seriation_sort") {
        # Usage of seriation sorting 
        # IMPORTANT! HERE I NEED THE SERIATION-LIBRARY. 
        # "seriate" AND "get_order" COME FROM IT!
        seriation_result <- seriation::seriate(binary_df, margin = 2)
        #print(seriation_result)
        order_to_go_by <- seriation::get_order(seriation_result)
      }
      
      if (sort == "jaccard_sort"){
        order_to_go_by <- jaccard(binary_df)
      }
    }
    
    
    # - - - ACTUAL SORTING SUBSECTION - - - - - - - - - - - - - - - - - - - - -
    # Get original names for later
    saved_colnames <- colnames(main_data)
    
    # 'sorting()' returns a list containing sorted main_data and sorted 
    # batch_data. It takes "order_to_go_by" as an argument. This is FALSE if
    # "sparsity_sort" is used and a list of length equal to the amount of
    # batches if another sorting approach is chosen
    main_desc <- sorting(
      main_data, 
      batch_list, 
      batch_data, 
      order_to_go_by, 
      prints
    )
    
    # Note: data.frame is important since elsewise it will be treated as a 
    # list element
    # Update main_data
    main_data <- data.frame(main_desc[1], check.names = FALSE)
    # Update batch_data
    batch_data <- data.frame(main_desc[2])
    # Update batch_list
    batch_list <- fetch_batch_overview(batch_data)
  }
  
  
  # ----- BLOCKING SECTION ----------------------------------------------------
  # Create empty "block_list" to be filled by 'blocking()' function
  block_list <- c()
  
  # Set "block_list" the same as "batch_list" in case blocking is not used
  if (is.null(block)) {
    block_list <- batch_list
  } 
  # Else create "block_list" based on user input
  else {
    number_batches <- tail(batch_list, n=1)
    # For example: valid input for "block" with 5 batch data would be 2, 3 or 4
    if (is.double(block) && block < number_batches && block > 1) {
      
      # PRINT-OUT
      if (prints == "standard" || prints == "all") {
        print("Blocking...")
      }
      
      calculated_block_list <- blocking(batch_list, block)
      # Update "block_list" with the calculated result
      block_list <- calculated_block_list
    
    } else {
      print("Please enter a valid number to be blocked by. Continuing without blocking.")
      block_list <- batch_list
      block <- NULL
    }
    
    # PRINT-OUT
    if (prints == "all") {
      print("The given batch listings for all samples:")
      print(batch_list)
      print(paste("Blocked (block = ",block ,"):", sep = ""))
      print(block_list)
    }
  }
  
  
  # ----- SPOTTING SECTION ----------------------------------------------------
  # PRINT-OUT
  if (prints == "standard" || prints == "all") {
    print("Preparing...")
  }

  # Create the "affiliation_list", the list of vectors
  if (ComBat_mode == 1 ||
    ComBat_mode == 3 ||
    algorithm == "limma") {
    affiliation_list <- spotting_missing_values(
      main_data, 
      batch_list, 
      block_list, 
      2, 
      prints
    )
  }
  if (ComBat_mode == 2 || ComBat_mode == 4) {
    affiliation_list <- spotting_missing_values(
      main_data, 
      batch_list, 
      block_list, 
      1, 
      prints
    )
  }
  
  orig_copy <- affiliation_list
  
  
  # ----- REMOVING UNIQUE COMBINATIONS SECTION --------------------------------
  # Removing all unique combinations from "affiliation_list". Currently 
  # untouchable by the user
  new_affiliation_list <- unique_removal(affiliation_list)
  
  # Update "affiliation_list"
  affiliation_list <- new_affiliation_list

# VERBOSITY #
############################################################### From here  
  
  killcount <- 0
  for (element in affiliation_list) {
    if (length(element) == 0){
      killcount <- killcount + 1
    }
  }
  print("Features with insufficient data to be considered:")
  print(killcount)
  
  # Keep for evaluation maybe for the runs without unique removal?
  new_orig <- list()
  for (element in orig_copy) {
    new_orig <- append(new_orig, toString(element))
  }
  #print(length(unique(new_orig)))
  
  
  string_affiliation_list <- list()
  for (element in affiliation_list) {
    string_affiliation_list <- append(string_affiliation_list, toString(element))
  }
  print("Sub-dataframes produced in total:")
  print(length(unique(string_affiliation_list)))
  
  # Error-check. Likely unnecessary now
  #for (element in string_affiliation_list){
  #  if (sum(string_affiliation_list == element) < 2){
  #    if (nchar(element) > 1){
  #      print("Error, the following element only appears once!")
  #      print(element)
  #    }
  #  }
  #}
  
############################################################### To here  
  

  # ----- SPLITTING SECTION ---------------------------------------------------
  # PRINT-OUT
  if (prints == "standard" || prints == "all") {
    print(paste("Splitting the data using", algorithm, "adjustment..."))
  }

  # Split up the "main_data" dataframe into the sub-dataframes
  cured_subdfs <- splitting(
    affiliation_list,
    main_data,
    batch_data,
    block_list,
    algorithm,
    ComBat_mode,
    block,
    prints
  )
  

  # ----- REBUILD SECTION -----------------------------------------------------
  # PRINT-OUT
  if (prints == "standard" || prints == "all") {
    print("Rebuilding...")
  }
  
  # Build the result file by rebuilding the cured matrix from all viable
  # sub-dataframes
  cured <- rebuild(cured_subdfs)
  
  
  # ----- RESORT SECTION ------------------------------------------------------
  # Resorting both main_data and cured
  if (sort == TRUE || sort == "sparsity_sort" || sort == "seriation_sort" || sort == "jaccard_sort") {
    
    # PRINT-OUT
    if (prints == "standard" || prints == "all") {
      print("Sorting back to normal...")
    }
    
    main_data <- main_data[ , saved_colnames]
    cured <- cured[ , saved_colnames]
  }
  
  
  # ----- WRITE-OUT SECTION ---------------------------------------------------
  # The cured_data.tsv file is written
  outfilename <- paste(output_file, "tsv", sep = ".")
  # 'unlink()' makes sure that a file with that ecaxt name gets deleted prior
  # to avoid conflicts/errors
  unlink(outfilename)
  write.table(cured, outfilename, sep = "\t", col.names = NA)


  # ----- VISUALIZATION SECTION -----------------------------------------------
  # Use the original "backup_batch_list" for the visualize functions
  batch_list <- backup_batch_list
  
  # Save a version of the logged cured for the harmonizR() return (since 
  # 'visual()', 'visual2()' and 'visual3' change it)
  original_cured <- cured
  
  # Reverse the log for visualizing
  main_data <- 2^main_data
  cured <- 2^cured

  if (plot == "featuremeans") {
    # PRINT-OUT
    if (prints == "standard" || prints == "all") {
      print("Visualizing feature means...")
    }

    original <- visual(main_data, batch_list)
    corrected <- visual(cured, batch_list)

    par(mfrow = c(1, 2))
    boxplot(original, main = "Original", las = 2)
    boxplot(corrected, main = "Corrected", las = 2)
  } else if (plot == "samplemeans") {
    # PRINT-OUT
    if (prints == "standard" || prints == "all") {
      print("Visualizing sample means...")
    }
    
    original <- visual2(main_data, batch_list)
    corrected <- visual2(cured, batch_list)

    lmts <- range(original, corrected)

    par(mfrow = c(1, 2))
    boxplot(original, main = "Original", las = 2, ylim = lmts)
    boxplot(corrected, main = "Corrected", las = 2, ylim = lmts)
  } else if (plot == "CV") {
    # PRINT-OUT
    if (prints == "standard" || prints == "all") {
      print("Visualizing CV...")
    }
    
    original <- visual3(main_data, batch_list)
    corrected <- visual3(cured, batch_list)

    lmts <- range(original, corrected)

    par(mfrow = c(1, 2))
    boxplot(original, main = "Original", las = 2, ylim = lmts)
    boxplot(corrected, main = "Corrected", las = 2, ylim = lmts)
  }


  # ----- END -----------------------------------------------------------------
  # PRINT-OUT
  if (prints == "standard" || prints == "all") {
    print("Termination.")
  }
  
  return(original_cured)
}