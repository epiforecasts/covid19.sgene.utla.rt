library(osfr)

# authorise uploading using PAT
#osf_auth(<PAT here>)

# project reference
project <- osf_retrieve_node("gjskf")

#upload processed results
osf_upload(project, path = "output", recurse = TRUE, progress = TRUE,
           conflicts = "overwrite")
