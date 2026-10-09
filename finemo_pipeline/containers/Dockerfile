FROM mambaorg/micromamba:2.3.2
ENV PIP_NO_CACHE_DIR=1
COPY --chown=$MAMBA_USER:$MAMBA_USER finemo_pipeline/containers/environment.yml /tmp/finemo-environment.yml
RUN micromamba install --yes --name base --file /tmp/finemo-environment.yml \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:/usr/local/nvidia/bin:$PATH \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN mkdir -p /opt/finemo_pipeline \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/finemo_pipeline
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER finemo_pipeline/scripts /opt/finemo_pipeline/scripts
RUN python -m pip check
WORKDIR /tmp
RUN --mount=type=bind,source=tests,target=/tests \
    python /tests/smoke_finemo_image.py
