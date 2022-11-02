library(plumber)
library(here)

pr(here("R", "plumber_fgsea.R")) %>%
  pr_run(host="0.0.0.0", port=8000)

