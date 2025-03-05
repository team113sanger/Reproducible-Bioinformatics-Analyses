# Recording and rebuilding your environment

This exercise will teach you how to use an `.renv` lock file to record your analsysis dependencies

## Prerequisites

- Docker installed on your machine
- Basic knowledge of R and its packages

## Create an Renv lock file
Open an R session or and Rstudio session and try installing a couple of packages

```R
## Install two mini packages
renv::install("glue")
renv::install("here")
## Create a snapshot that records your new added packages
renv::snapshot()
```

## Reconstruct the environment from the lockfile
```bash
cd exercises/02_rebuilding_your_environment
docker build -t example_renv . 
```
## Rebuild your environment with 
```bash
docker run -it example_renv R
```
This should open an R terminal where in theory you could 

```
docker run --rm -p 8789:8789 example_renv sudo sudo rstudio-server start
```








