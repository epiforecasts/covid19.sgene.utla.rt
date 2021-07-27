library(osfr)

project <- osf_retrieve_node("gjskf")
files <- osf_ls_files(project)
osf_download(files, path = ".", recurse = TRUE, progress = TRUE,
             conflicts = "overwrite")
