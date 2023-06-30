FROM lukauskas/marcsenvironment

# Install Helvetica fonts
# If you don't have Helvetica, comment-out the next two lines, and the echo 'font.sans-serif ... line below
COPY fonts/*.ttf  /usr/local/share/fonts/
RUN fc-cache

RUN mkdir -p /root/.config/matplotlib \
    && echo 'backend : Agg' >> /root/.config/matplotlib/matplotlibrc \
    && echo 'font.sans-serif : helvetica' >> /root/.config/matplotlib/matplotlibrc

# Set up environment

ENV SNAPANALYSIS_RAW_DATA="/data/raw" \
    SNAPANALYSIS_EXTERNAL_DATA="/data/external" \
    SNAPANALYSIS_INTERIM_DATA="/data/interim" \
    SNAPANALYSIS_OUTPUT="/output" \
    SNAPANALYSIS_LOG_CONFIG="/conf/logging.yaml"

RUN mkdir -p $SNAPANALYSIS_RAW_DATA \
    && mkdir -p $SNAPANALYSIS_EXTERNAL_DATA \
    && mkdir -p $SNAPANALYSIS_INTERIM_DATA \
    && mkdir -p $SNAPANALYSIS_OUTPUT \
    && mkdir -p /log

# Copy raw and external data to their place
COPY data/raw /data/raw
COPY data/downloaded /data/external
COPY logging.docker.yaml /conf/logging.yaml

#### Main analysis #######################
ARG BUILD_DATE
ARG COMMIT

LABEL org.label-schema.vcs-ref="$COMMIT" \
      org.label-schema.build-date="$BUILD_DATE"

RUN mkdir /build && echo "DATE=$BUILD_DATE\nCOMMIT=$COMMIT" >> /build/info.txt

COPY src /snap
RUN  pip install --no-deps -e /snap/

# Build current data
RUN cd \
    # Reads the output from MaxQuant
    && python -m snapanalysis.preprocessing.raw.extract \
    # Cleans up the identifiers, etc.
    && python -m snapanalysis.preprocessing.cleanup.main \
    # Some metadata
    && python -m snapanalysis.preprocessing.pulldown_metadata \
    # Imputation & othersi
    && python -m snapanalysis.models.enrichment.generate \
    # Metadata about proteins (Domains, complexes, etc.)
    && python -m snapanalysis.preprocessing.protein_metadata \
    # Read biogrid dataset
    && python -m snapanalysis.external.interactions.biogrid \
    # Train network
    && python -m snapanalysis.models.network.training \
    && python -m snapanalysis.models.network.drawall \
    && python -m snapanalysis.models.ptm_response.main

#### Figures #######################

COPY extra-figures-notebooks /notebooks

RUN \
    # Metadata
    jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/metadata/Figure Nucleosome Map.ipynb" \
    # Network table
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/network/Table Network.ipynb" \
    # Preprocessing
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/preprocessing/Figure Preprocessing Workflow.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/preprocessing/Figure Imputation Stats.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/preprocessing/Table Pull-Down Data.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/preprocessing/Table Pull-Down Heatmap.ipynb" \
    # PTM response figures ()
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Complexes Barplots.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=2400 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Complexes Heatmaps.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=2400 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Method.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Network Projection.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Pairwise.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Pairwise Gap Plot.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Volcanoes.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Table Complexes.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Table Proteins.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Table Response Clustering.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/ptm-response/Figure Total Counts.ipynb" \
    # Scatterplots
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/scatterplots/Figure Annotated Scatterplots.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/scatterplots/Figure Scatterplots Highlights.ipynb" \
    && jupyter nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute --output /dev/null "/notebooks/scatterplots/Figure Scatterplot Grids.ipynb" 


# Based on neat script here: https://www.r-bloggers.com/list-of-user-installed-r-packages-and-their-versions/
RUN Rscript -e 'p <- as.data.frame(installed.packages()[,c(1,3:4)]); rownames(p) <- NULL; p <- p[is.na(p$Priority),1:2,drop=FALSE]; print(p, row.names=FALSE)' > $SNAPANALYSIS_OUTPUT/r-packages.txt \
    && pip freeze > $SNAPANALYSIS_OUTPUT/python-packages.txt \
    && cp /log/snap.log $SNAPANALYSIS_OUTPUT/marcs.log

# Generate data for the web
COPY ./for-web /for-web
RUN python /for-web/gene_select.py \
    && python /for-web/predictors.py \
    && python /for-web/collect_downloads.py \
    && mkdir -p /output/web/api/precompiled \
    && python /for-web/api/extract_gene_meta.py \
    && python /for-web/api/extract_matrix.py \
    && python /for-web/api/extract_predictors.py \
    && python /for-web/api/extract_annotations.py \
    && python /for-web/api/extract_networks.py \
    && python /for-web/api/extract_ptm_response.py \
    && python /for-web/api/extract_build_version.py
