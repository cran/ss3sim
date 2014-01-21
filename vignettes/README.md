This folder contains the vignette document for the ss3sim R package.

The knitr source file `ss3sim-vignette.Rnw` is compiled by R when installing the package.
The resulting PDF file is then copied to `inst/doc/`.

You can download a PDF version of the vignette [here](https://dl.dropboxusercontent.com/u/254940/ss3sim-vignette.pdf).

You can build the PDF vignette yourself by running the following command in R:
```{r}
# install.packages("knitr") # if needed
knitr::knit2pdf("ss3sim-vignette.Rnw")
```
