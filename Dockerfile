FROM rstudio/plumber:v1.2.0

RUN R -e "install.packages('remotes')"

# packages versions are associated with versions of bioconductor
RUN R -e "remotes::install_bioc('3.20/fgsea')" # installs 1.32.2
RUN R -e "remotes::install_bioc('3.20/limma')" # installs 3.62.2
RUN R -e "remotes::install_bioc('3.20/edgeR')" # installs 4.4.1
RUN R -e "install.packages('readr', version = '2.1.5')"
RUN R -e "install.packages('here', version = '1.0.1')"

WORKDIR /app
COPY ./R ./R
COPY ./data ./R/data

CMD ["/app/R/plumber_fgsea.R"]
