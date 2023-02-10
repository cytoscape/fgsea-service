# fgsea-service

- To create docker container that hosts the fgsea-service:
  ```
  docker build -t fgsea .
  docker run --rm -p 8000:8000 fgsea
  ```

- To create docker container for development:
  ```
  docker build -t mike_rstudio -f Dockerfile_bioconductor .
  docker run -v /home/mkucera/fgsea-service:/home/rstudio/fgsea-service -e USERID=$UID -e PASSWORD=password -p 8888:8787 mike_rstudio
  ```
