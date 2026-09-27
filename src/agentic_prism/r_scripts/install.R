# Optional pinned packages only; this script does not install the R interpreter.
a<-commandArgs(TRUE);lib<-normalizePath(a[1],mustWork=FALSE)
dir.create(lib,recursive=TRUE,showWarnings=FALSE);.libPaths(c(lib,.libPaths()))
pins<-read.csv(a[2],stringsAsFactors=FALSE)
repo<-'https://cloud.r-project.org'
av<-available.packages(repos=repo)
for(i in seq_len(nrow(pins))) {
 p<-pins$package[i];v<-pins$version[i]
 if(requireNamespace(p,quietly=TRUE) && as.character(packageVersion(p))==v) next
 # Binary repository metadata must advertise the pinned version; otherwise use
 # the exact source release, then its archive URL. No unpinned dependency solve.
 if(p %in% rownames(av) && av[p,'Version']==v) {
   install.packages(p,lib=lib,repos=repo,dependencies=FALSE)
 } else {
   dest<-tempfile(fileext='.tar.gz')
   urls<-c(paste0(repo,'/src/contrib/',p,'_',v,'.tar.gz'),paste0(repo,'/src/contrib/Archive/',p,'/',p,'_',v,'.tar.gz'))
   ok<-FALSE
   for(url in urls) if(isTRUE(tryCatch({download.file(url,dest,quiet=TRUE);TRUE},error=function(e)FALSE))) {ok<-TRUE;break}
   if(!ok) stop(paste('Cannot download pinned',p,v))
   install.packages(dest,lib=lib,repos=NULL,type='source',dependencies=FALSE)
 }
 if(!requireNamespace(p,quietly=TRUE) || as.character(packageVersion(p))!=v) stop(paste('Pinned installation failed:',p,v))
}
cat('Optional R package pins verified\n')
