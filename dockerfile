FROM mambaorg/micromamba as teloseq
USER root
RUN apt-get update && apt-get install -y procps python3-pip
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/*.yaml /

FROM teloseq AS seqkit
RUN micromamba install -y -n base --file /seqkit.yaml && \
        micromamba clean --all --yes

#FROM augustus AS eggnog_mapper
#RUN micromamba install -y -n base --file /eggnog_mapper.yaml && \
#        micromamba clean --all --yes
#RUN mkdir /opt/conda/lib/python3.10/site-packages/data
#SHELL ["micromamba", "run", "-n", "base", "/bin/bash", "-c"]
#RUN download_eggnog_data.py -y
