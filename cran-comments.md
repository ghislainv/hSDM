# hSDM v1.4.4

We have remove the cause of the following error ("%" in .bib file):

Check Details

Version: 1.4.3
Check: re-building of vignette outputs
Result: ERROR
    Error(s) in re-building vignettes:
    --- re-building ‘hSDM.Rmd’ using rmarkdown
    --- finished re-building ‘hSDM.Rmd’
    
    --- re-building ‘publications.Rmd’ using rmarkdown
    Error reading bibliography file bib/biblio-publications.bib:
    (line 50, column 3):
    unexpected '%'
    expecting white space
    Error: processing vignette 'publications.Rmd' failed with diagnostics:
    pandoc document conversion failed with error 25
    --- failed re-building ‘publications.Rmd’
    
    SUMMARY: processing the following file failed:
     ‘publications.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
Flavors: r-devel-linux-x86_64-fedora-gcc, r-release-windows-x86_64, r-oldrel-windows-x86_64 

# hSDM v1.4.3

We have added \alias{hSDM-package} to file 'hSDM/man/hSDM-package.Rd' as suggested by Kurt Hornik by email on 18/08/2023.

We have removed the links to the following URLs which were invalid:

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>’
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1139/cjfas-2018-0148
      From: inst/doc/publications.html
      Status: 403
      Message: Forbidden
    URL: https://doi.org/10.1670/14-075
      From: inst/doc/publications.html
      Status: Error
      Message: code d’erreur libcurl 60:
        	SSL certificate problem: unable to get local issuer certificate
    URL: https://hal.science/tel-02519161
      From: inst/doc/publications.html
      Status: Error
      Message: code d’erreur libcurl 28:
        	Operation timed out after 60001 milliseconds with 0 bytes received
