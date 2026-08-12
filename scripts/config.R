# Run scripts from the repository root. Environment variables may override
# any default directory without changing the analysis code.
project_dir <- normalizePath(
  Sys.getenv("RNA_LNP_PROJECT_DIR", unset = getwd()),
  mustWork = FALSE
)
input_dir <- normalizePath(
  Sys.getenv("RNA_LNP_INPUT_DIR", unset = file.path(project_dir, "data", "inputs")),
  mustWork = FALSE
)
primary_dir <- normalizePath(
  Sys.getenv("RNA_LNP_PRIMARY_DIR", unset = file.path(project_dir, "data", "primary")),
  mustWork = FALSE
)
core_output_dir <- normalizePath(
  Sys.getenv("RNA_LNP_CORE_OUTPUT_DIR", unset = file.path(project_dir, "results", "core")),
  mustWork = FALSE
)
revision_output_dir <- normalizePath(
  Sys.getenv("RNA_LNP_REVISION_OUTPUT_DIR", unset = file.path(project_dir, "results", "revision")),
  mustWork = FALSE
)
figure_output_dir <- normalizePath(
  Sys.getenv("RNA_LNP_REVISION_FIGURE_DIR", unset = file.path(project_dir, "results", "figures")),
  mustWork = FALSE
)

first_existing <- function(...) {
  candidates <- unlist(list(...), use.names = FALSE)
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("Required input not found: ", paste(candidates, collapse = ", "))
  }
  found[[1]]
}
