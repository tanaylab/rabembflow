# script for downloading and initializing the single-cell RNA-seq database for the metacell package

if(!dir.exists("output")) {
  dir.create("output")
} 

download.file("https://rabemb.s3.eu-west-1.amazonaws.com/oc_mm_scdb_embexe.tar.gz","oc_mm_scdb_embexe.tar.gz")

system("tar -xvf oc_mm_scdb_embexe.tar.gz output/scrna_db")

file.remove("oc_mm_scdb_embexe.tar.gz")