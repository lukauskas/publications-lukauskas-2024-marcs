FROM python:3.6-stretch

RUN echo "deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y --force-yes cmake \
                                      git \
                                      graphviz \
                                      libblas-dev \
                                      libgsl-dev \
                                      libhdf5-dev \
                                      libudunits2-dev \
                                      liblapack-dev \
                                      openjdk-8-jdk \
                                      python-dev \
                                      r-base \
                                      ttf-freefont \
    && apt-get -yq autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install some needed R packages
RUN Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("BiocManager")' \
     && Rscript -e 'BiocManager::install(c("limma", "statmod", "AnnotationDbi", "fgsea", "reactome.db", "impute", "GO.db", "preprocessCore", "enrichplot", "ReactomePA", "clusterProfiler", "meshes", "org.Hs.eg.db", "GOSemSim", "minet"))' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("precrec")' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("ggplot2")' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("ggrepel")' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("corpcor")' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("fdrtool")' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("ppcor")' \
     && Rscript -e 'chooseCRANmirror(graphics=FALSE, ind=1); install.packages("matrixcalc")'

# Cache Reactome and GO Data so it's faster when we need to run it
RUN Rscript -e 'library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); enrichGO(c(4312, 8318, 10874, 55143, 55388, 991), OrgDb=org.Hs.eg.db, ont="MF"); enrichPathway(c(4312, 8318, 10874, 55143, 55388, 991));'

# python setup
ADD requirements.txt /requirements.txt
RUN pip install --upgrade pip \
    && pip install --upgrade -r /requirements.txt \
    && rm /requirements.txt \
    && rm -rf ~/.cache/pip \
    && rm -rf /src

# Put forceatlas jar here too
RUN mkdir /jars/ \
    && wget -q https://github.com/lukauskas/gephi-toolkit-forceatlas2-standalone/releases/download/v0.0.5/gephi-toolkit-forceatlas-0.0.5-jar-with-dependencies.jar -O /jars/gephi-forceatlas2.jar
