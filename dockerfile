FROM mambaorg/micromamba AS tarpon
USER root
RUN apt-get update && apt-get install -y procps python3-pip cmake r-base libcurl4-openssl-dev r-cran-lme4
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/*.yaml /
RUN micromamba install -y -n base --file /tarpon.yaml && \
        micromamba clean --all --yes
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('viridis')"
RUN Rscript -e "install.packages('ggpointdensity')"
RUN Rscript -e "install.packages('ggpubr')"
