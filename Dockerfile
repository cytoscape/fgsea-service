FROM rstudio/plumber

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('fgsea')"
RUN R -e "install.packages('readr')"
RUN R -e "install.packages('here')"

WORKDIR /app
COPY ./R ./R
COPY ./data ./R/data

CMD ["/app/R/plumber_fgsea.R"]
