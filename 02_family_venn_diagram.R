script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}
source(file.path(script_dir, "_shared.R"), encoding = "UTF-8")

ensure_packages(c("eulerr", "grid"))

project_dir <- project_dir_from_script(script_dir)
data_dir <- data_dir_from_env(project_dir)
output_dir <- output_dir_from_env(project_dir, "family_venn")

cmbl_file <- Sys.getenv("CMBL_FAMILY_FILE", unset = file.path(data_dir, "CMBL_family_data.csv"))
ncbi_file <- Sys.getenv("NCBI_FAMILY_FILE", unset = file.path(data_dir, "NCBI_family_data.csv"))
morphology_file <- Sys.getenv(
  "MORPHOLOGY_FAMILY_FILE",
  unset = file.path(data_dir, "zhejiang_morphology_family_data.csv")
)

read_family_column <- function(path, family_column = "family") {
  data <- read_csv_checked(path, header = TRUE, check.names = FALSE)

  if (!family_column %in% names(data)) {
    stop("Column not found in ", path, ": ", family_column, call. = FALSE)
  }

  clean_character_set(data[[family_column]])
}

family_sets <- list(
  CMBL = read_family_column(cmbl_file),
  NCBI = read_family_column(ncbi_file),
  Morphology = read_family_column(morphology_file)
)

fit <- eulerr::euler(family_sets)

sink(file.path(output_dir, "family_level_euler_fit.txt"))
print(fit)
sink()

pdf(file.path(output_dir, "family_level_area_proportional_venn.pdf"), width = 7, height = 6)
plot(
  fit,
  fills = list(
    fill = c("#E68D3D", "#E26472", "#6270B7"),
    alpha = 0.55
  ),
  edges = list(
    col = "grey25",
    lwd = 1.2
  ),
  labels = list(
    font = 2,
    fontsize = 15,
    col = "grey10"
  ),
  quantities = list(
    fontsize = 14,
    font = 2,
    col = "grey10"
  ),
  legend = FALSE,
  main = NULL
)
dev.off()

print(fit)

