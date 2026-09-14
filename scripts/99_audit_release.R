# Read-only release audit.
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset=getwd()), winslash="/", mustWork=TRUE)
manifest_path <- file.path(root,"documentation","FILE_MANIFEST_SHA256.csv")
manifest <- read.csv(manifest_path, check.names=FALSE)
if(anyDuplicated(manifest$file) || any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",manifest$file))) stop("Unsafe or duplicated manifest paths.",call.=FALSE)
sha <- function(p) digest::digest(file=p,algo="sha256",serialize=FALSE)
actual <- vapply(file.path(root,manifest$file),sha,character(1))
if(any(actual != manifest$sha256)) stop("Manifest mismatch: ",paste(manifest$file[actual != manifest$sha256],collapse=", "),call.=FALSE)
files <- list.files(root,recursive=TRUE,full.names=TRUE,all.files=TRUE); files <- files[!dir.exists(files)]
relative <- substring(files,nchar(root)+2L)
relative <- relative[!startsWith(relative,".git/")]
if(!setequal(relative,c(manifest$file,"documentation/FILE_MANIFEST_SHA256.csv"))) stop("Unlisted or missing package file.",call.=FALSE)
if(any(grepl("\\.(rds|rdata|xpt|parquet|gz)$",files,ignore.case=TRUE))) stop("A data-like binary file is present.",call.=FALSE)
message("Release audit passed: ",nrow(manifest)," files matched SHA-256; no data-like binary file present.")

