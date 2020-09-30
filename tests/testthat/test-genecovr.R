pkg_path <- dirname(file.path(system.file(package = "genecovr")))
create_file <- function(fn) system.file("extdata", fn, package = "genecovr")

data <- as.data.frame(rbind(
    c("nonpol", create_file("transcripts2nonpolished.psl"),
      create_file("nonpolished.fai"), create_file("transcripts.fai")),
    c("pol", create_file("transcripts2polished.psl"),
      create_file("polished.fai"), create_file("transcripts.fai"))
))
colnames(data) <- NULL

tmp <- tempdir(check = TRUE)
withr::local_dir(tmp)
withr::local_file(write.csv(data, file = "assemblies.csv", row.names = FALSE))
withr::local_temp_libpaths()
devtools::install(pkg_path, quick = TRUE, upgrade = "never")
withr::local_path(file.path(.libPaths()[1], "genecovr", "bin"))

test_that("genecovr runs as expected", {
    system(paste(paste0("R_LIBS_USER=", .libPaths()[1]),
                 "genecovr", "assemblies.csv"))
    expect_equal(length(list.files(pattern = "*.png")), 13)
    expect_equal(length(grep("Rplots", list.files(pattern = "*.pdf"),
                             invert = TRUE), 13))
    expect_equal(length(list.files(pattern = "*.csv.gz")), 4)
})

withr::defer(unlink(dir, recursive = TRUE))
