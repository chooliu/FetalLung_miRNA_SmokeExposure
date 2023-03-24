# ==============================================================================
# 02_table1.R
# summarize demographics by IUS
# ==============================================================================




# table 1 helper fxns ----------------------------------------------------------

summaryCountPercent <-
  function(x, values, count_NAs = T,
         digits = 1,
         NA_percent_in_brackets = F,
         NA_count_in_brackets = F,
         NA_count_in_brackets_alt = T,
         fuzzy = F, inverse = F) {
  
  n_NA <- sum(is.na(x))
  ratio_miss <- n_NA/length(x)*100
  
  if (n_NA != 0) {
    warning(paste0(n_NA,
                   " values are missing/NA in input vector"))
  }
  
  if (n_NA == length(x)) {
    warning("all values missing!")
    return("-")
  }
  
  # exact match (default)
  n <- sum(x %in% values, na.rm = T)
  tot <- ifelse(count_NAs, length(x), length(x) - n_NA)
  
  # fuzzy matching
  if (fuzzy) {
    search_term <-
      paste(values, collapse = "|")
    n <- grepl(search_term, x, ignore.case = T) %>% sum
  }
  
  # inverse
  if (inverse) {
    n <- tot - n
  }
  
  output <- paste0(n, " (", formatC(n / tot * 100, format = "f", digits = digits),
                   "%)")# (", n, "/", tot, ")")
  
  if (NA_count_in_brackets & n_NA != 0) {
    output <-
      paste0(output, " [", n_NA, "]")
  }
  
  if (NA_percent_in_brackets & n_NA != 0) {
    output <-
      ratio_miss %>%
      formatC(., format = "f", digits = digits) %>%
      paste0(output, " [", ., "%]")
  }
  
  return(output)
}



summaryMedianIQR <-
  function(x, digits = 1, na.rm = F,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = F,
           NA_count_in_brackets_alt = T,
           symbol = NULL) {
    
    if ( !(typeof(x) %in% c("double", "integer")) ) {
      warning("input vector not numeric")
      x <- as.numeric(x)
    }
    
    tot <- length(x)
    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA/length(x)*100
    
    if (n_NA != 0) {
      warning(paste0(n_NA,
                     " values are missing/NA in input vector"))
    }
    
    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }
    
    output <- quantile(
      as.numeric(x),
      c(0.25, 0.5, 0.75), na.rm = na.rm
    )
    
    output <- formatC(output, format = "f", digits = digits)
    
    output <- ifelse(is.null(symbol),
                     paste0( output[2], " (", output[1], ", ", output[3], ")"),
                     paste0( output[2], symbol, " (", output[1], symbol, ", ",
                             output[3], symbol, ")"))
    
    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }
    
    if (NA_count_in_brackets_alt & !NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [n = ", tot - n_NA, "]")
    }
    
    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }
    
    
    return(output)
    
  }






# make table -------------------------------------------------------------------

Table1 <-
  metadata_miRNA %>%
    mutate(IUS = if_else(SmokeBinary == 1, "IUS-exposed", "unexposed")) %>%
    group_by(IUS) %>%
    summarize(`N` = n() %>% as.character,
              `Age (dpc)` = summaryMedianIQR(Age),
              `Placental Cotinine (ng/g)` = summaryMedianIQR(Cotinine),
              `Sex (Male)` = summaryCountPercent(Sex, "TRUE")) %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    pivot_wider(names_from = IUS)


metadata_miRNA %>%
  group_by(SmokeBinary, Sex) %>%
  tally()

