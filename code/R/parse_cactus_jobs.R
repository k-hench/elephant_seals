instructions <- readLines("../results/cactus/cactus_instructions.sh")

round_idx <- 0L
round_name <- "preprocessor"
job_idx <- 0

script_dir <- 'sh/cactus'
is_title <- \(str){grepl("Preprocessor|Round|merging",str)}

line_type <- \(str){
  if(is_title(str)){
    return("title")
  } else if (substr(str, 1, 1) == "#") {
    return("comment")
  } else if (str == "") {
    return( "empty" )
  } else {
    return ("code")
  }
}

inventory <- c()

j <- 0L; d <- 0L
for( i in seq_along(instructions)){
  ln_cur <- instructions[i]
  ln_type <- line_type(ln_cur)
  
  if(ln_type == "title") {
    round_idx <- round_idx + 1L
    j <- 1L
    d <- d + 1L
  } else if (ln_type == "empty") {
    j <- j + 1L
  } else if (ln_type == "code") {
    inventory[d] <- j
  }
}

write.table(data.frame(round = seq_along(inventory),
                       n_jobs = inventory), 
            file = "../results/cactus/job_inventory.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

n_rounds_digits <- floor(log10(length(inventory))) + 1L
inventory_digits <- floor(log10(inventory)) + 1L

round_idx <- 0L
round_name <- "preprocessor"
job_idx <- 0
job_table <- data.frame(round_idx = c(),
                        round = c(),
                        job_idx = c(),
                        job = c())
code_nr <- 0L

for( i in seq_along(instructions)){
  ln_cur <- instructions[i]
  ln_type <- line_type(ln_cur)
  
  if(ln_type == "title") {
    round_idx <- round_idx + 1L
    job_idx <- 1L
    dir.create(paste0(script_dir, "/round_", 
                      sprintf(paste0("%0", n_rounds_digits,".0f"), round_idx)),
               recursive = TRUE)
  } else if (ln_type == "empty") {
    job_idx <- job_idx + 1L
    code_nr <- 0L
  } else if (ln_type == "code") {
    code_nr <- code_nr + 1L 
    if(code_nr == 1L){
      job_table <- rbind(
        job_table,
        data.frame(round_idx = round_idx,
                   round = paste0("round_", sprintf(paste0("%0", n_rounds_digits,".0f"), round_idx)),
                   job_idx = job_idx,
                   job = paste0("job_", sprintf(paste0("%0", inventory_digits[round_idx],".0f"), job_idx),".sh"))
        )
    }
    write(ln_cur,
          file = paste0(
            script_dir, "/round_", sprintf(paste0("%0", n_rounds_digits,".0f"), round_idx), "/job_",
            sprintf(paste0("%0", inventory_digits[round_idx],".0f"), job_idx),".sh"),
          append = TRUE)
  }
}

write.table(job_table, 
            file = "../results/cactus/job_list.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

