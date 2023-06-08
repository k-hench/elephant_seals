base_dir <- "sh/cactus/"
rounds <- dir("sh/cactus/", pattern = "round")

steps_table <- data.frame(round_idx = c(0),
                              round = c("round_0"),
                              job_idx = c(0),
                              job = c("job_0.sh"),
                              step_idx = c(0),
                              step = c("step_0"))

final_step_table <- data.frame(round_idx = c(0),
                              round = c("round_0"),
                              job_idx = c(0),
                              job = c("job_0.sh"),
                              step_idx = c(0),
                              step = c("step_0"))

for(round_idx in seq_along(rounds)){
  round_base <- paste0(base_dir, rounds[[round_idx]])
  jobs <- dir(round_base, pattern = "sh")
  n_jobs <- length(jobs)
  n_jobs_digits <- floor(log10(n_jobs)) + 1L

  for(job_idx in seq_along(jobs)){
    steps <- readLines(paste0(round_base,"/", jobs[[job_idx]]))
    jobs_dir <- paste0(round_base,"/job_",
                       sprintf(paste0("%0", n_jobs_digits,".0f"), job_idx))
    
    dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
    
    n_steps <- length(steps)
    n_steps_digits <- floor(log10(n_steps)) + 1L
    
    for(step_idx in seq_along(steps)){
      s_cur <- step_idx
      write(steps[[step_idx]],
            file = paste0(
              jobs_dir, "/step_",
              sprintf(paste0("%0", n_steps_digits,".0f"), step_idx),".sh"))
    
    steps_table <- rbind(
      steps_table,
      data.frame(round_idx = round_idx,
                 round = c(rounds[[round_idx]]),
                 job_idx = job_idx,
                 job = c(jobs[[job_idx]]),
                 step_idx = s_cur,
                 step = c(paste0("step_", sprintf(paste0("%0", n_steps_digits,".0f"), s_cur)))
                 ))
    }

    final_step_table <- rbind(
      final_step_table,
      data.frame(round_idx = round_idx,
                 round = c(rounds[[round_idx]]),
                 job_idx = job_idx,
                 job = c(jobs[[job_idx]]),
                 step_idx = s_cur,
                 step = c(paste0("step_", sprintf(paste0("%0", n_steps_digits,".0f"), s_cur)))
                 ))
  }
}

write.table(steps_table, 
            file = "../results/cactus/job_list.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(final_step_table, 
            file = "../results/cactus/job_final_steps.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)