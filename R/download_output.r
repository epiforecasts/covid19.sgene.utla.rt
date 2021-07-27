library(osfr)

project <- osf_retrieve_node("gjskf")

osf_download(project, path = "output", recurse = TRUE, progress = TRUE,
             conflicts = "overwrite")
