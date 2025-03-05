# Docker 101: Extending a Rocker Image

This exercise will teach you how to extend a Docker image for R by adding additional packages to an existing Rocker container.

## Prerequisites

- Docker installed on your machine
- Basic knowledge of R and its packages


## The Base Dockerfile

We're starting with a basic Dockerfile that uses the `rocker/tidyverse:4.2.0` image and adds a few common R packages:

```dockerfile
FROM rocker/tidyverse:4.2.0
# Install additional packages from CRAN
RUN install2.r --error \
 rmarkdown \
 shiny \
 knitr
```

## Your Task

Extend the provided Dockerfile to include two additional R packages:
1. `glue` - A string interpolation package for R
2. `here` - A package that helps with file path handling in projects

## Instructions

1. Create a new file named `Dockerfile` in an empty directory
2. Copy the base Dockerfile content into your file
3. Modify the Dockerfile to add the required packages
4. Build the Docker image
5. Test your image by running R and verifying the packages are installed

## Hints

- The `install2.r` command is provided by the Rocker image and is a convenient way to install R packages
- Make sure to maintain the existing structure of the Dockerfile
- You can add the packages to the existing `install2.r` command rather than creating a new one

## Building Your Docker Image

Once you've modified the Dockerfile, build your image using:

```bash
docker build -t my-r-environment .
```

## Testing Your Image

Test your image with:

```bash
docker run --rm -it my-r-environment R -e "library(glue); library(here); cat('Packages loaded successfully!\n')"
```

## Solution

<details>
<summary>Click to reveal the solution</summary>

```dockerfile
FROM rocker/tidyverse:4.2.0
# Install additional packages from CRAN
RUN install2.r --error \
 rmarkdown \
 shiny \
 knitr \
 glue \
 here
```

</details>

## Next Steps

After completing this exercise, you might want to:
- Add more complex dependencies to your Docker image
- Create a custom R script and run it inside the container
- Mount a volume to share data between your host and the container
- Learn how to push your image to Docker Hub or a private registry

Happy Dockerizing!