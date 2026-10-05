FROM mambaorg/micromamba:2.3.2

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

# WDL runners can replace ENTRYPOINT. Keep the environment active through PATH.
ENV PATH=/opt/conda/bin:/usr/local/nvidia/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/lib/python3.10/site-packages/tensorflow:/usr/local/nvidia/lib:/usr/local/nvidia/lib64 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    VARIANT_SCORER_DIR=/opt/variant-scorer/src \
    TF_FORCE_GPU_ALLOW_GROWTH=true \
    PYTHONUNBUFFERED=1
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN mkdir -p /opt/chrombpnet /opt/variant-scorer \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/chrombpnet /opt/variant-scorer
USER $MAMBA_USER
RUN git -C /opt/variant-scorer init \
    && git -C /opt/variant-scorer remote add origin https://github.com/kundajelab/variant-scorer.git \
    && git -C /opt/variant-scorer fetch --depth 1 origin 0e1e34199e63112aa618748bb79a206fc491300a \
    && git -C /opt/variant-scorer checkout --detach FETCH_HEAD \
    && test "$(git -C /opt/variant-scorer rev-parse HEAD)" = 0e1e34199e63112aa618748bb79a206fc491300a \
    && rm -rf /opt/variant-scorer/.git
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /opt/chrombpnet/scripts
RUN python /opt/chrombpnet/scripts/link_cuda.py \
    && python -m pip check
WORKDIR /tmp
