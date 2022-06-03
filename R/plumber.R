library(plumber)
library(here)

pr(here("R", "plumber_fgsea.R")) %>%
  pr_run(port=8000)
