# Recording and rebuilding your environment

This exercise will teach you how to use an `.renv` lock file to record your analysis dependencies


## Create an Renv lock file
Open an R session in the terminaln or and Rstudio session (see Jamie's demo) and try installing a couple of packages with `renv`

```R
## Install two mini packages
renv::install("glue")
renv::install("here")
## Create a snapshot that records your new added packages
renv::snapshot()
```

## Reconstruct the environment from the lockfile
```bash
cd exercises/02_rebuilding_your_R_environment
docker build -t example_renv . 
```
## Rebuild your environment with 
```bash
docker run -it example_renv R
```
This should open an R terminal where in theory you could 









