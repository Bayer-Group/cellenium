FROM postgres:15.1
RUN apt-get update && apt-get install -y wget unzip lz4 && \
    wget -O awscliv2.zip "https://awscli.amazonaws.com/awscli-exe-linux-$(arch).zip" && \
    unzip awscliv2.zip && \
    aws/install && \
    rm -r awscliv2.zip aws && \
    wget -O miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(arch).sh && \
    chmod a+x miniforge.sh && \
    ./miniforge.sh -b && \
    chmod a+x /root && \
    chmod a+x /root/miniforge3
COPY environment.yml .
RUN . /root/miniforge3/bin/activate && \
     conda env create -f environment.yml
COPY docker-entrypoint.sh /usr/local/bin/
RUN ln -s /bin/mkdir /usr/bin/mkdir
RUN bash -c "source /root/miniforge3/bin/activate postgres_python && git clone https://github.com/vibhorkum/pg_background.git"
RUN bash -c "source /root/miniforge3/bin/activate postgres_python && cd pg_background && make"
RUN bash -c "source /root/miniforge3/bin/activate postgres_python && cd pg_background && make install"
COPY s3_backup.sh .
RUN chmod a+x /usr/local/bin/docker-entrypoint.sh
