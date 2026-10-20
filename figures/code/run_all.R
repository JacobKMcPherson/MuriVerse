# Run the full figure pipeline: Rscript figures/code/run_all.R
source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))

for (script in c("04_data_muriverse.R",
                 "04b_data_muriverse_strict.R",
                 "05_data_halflife_model.R", "fig2.R", "fig3.R",
                 "06_data_polymyxin.R", "fig4.R",
                 "02_data_preprocessing.R", "supp_fig1.R", "supp_fig2.R")) {
  message("-> ", script)
  source(file.path(DIR_CODE, script))
}

sessionInfo() |> capture.output() |> writeLines(file.path(DIR_TABS, "sessionInfo.txt"))
