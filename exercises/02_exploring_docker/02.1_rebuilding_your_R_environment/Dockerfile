FROM rocker/tidyverse:4.2.0

# Install additional packages from CRAN
RUN install2.r --error renv

RUN R -e "renv::restore()"

