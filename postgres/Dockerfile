FROM postgres:15.1
RUN apt-get update && apt-get install -y wget && \
    wget -O miniforge.sh https://github.com/conda-forge/miniforge/releases/download/22.9.0-2/Miniforge3-22.9.0-2-Linux-$(arch).sh && \
    chmod a+x miniforge.sh && \
    ./miniforge.sh -b && \
    chmod a+x /root && \
    chmod a+x /root/miniforge3
COPY environment.yml .
RUN . /root/miniforge3/bin/activate && \
     conda env create -f environment.yml
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod a+x /usr/local/bin/docker-entrypoint.sh
