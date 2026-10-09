FROM mambaorg/micromamba:2.3.2

# CUDA wheels are large. Do not retain their downloaded archives in the image.
ENV PIP_NO_CACHE_DIR=1

COPY --chown=$MAMBA_USER:$MAMBA_USER motif_pipeline/containers/environment.yml /tmp/motif-environment.yml
RUN micromamba install --yes --name base --file /tmp/motif-environment.yml \
    && micromamba clean --all --yes

# Keep CUDA and the Python environment available when Cromwell replaces ENTRYPOINT.
ENV PATH=/opt/conda/bin:/usr/local/nvidia/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/lib/python3.10/site-packages/tensorflow:/usr/local/nvidia/lib:/usr/local/nvidia/lib64 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    TF_FORCE_GPU_ALLOW_GROWTH=true \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/opt/chrombpnet-upstream \
    MPLBACKEND=Agg \
    CHROMBPNET_COMMIT=09938fdb4397ec0006510e5251e48920a505d4de
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN mkdir -p /opt/motif_pipeline /opt/chrombpnet-upstream \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/motif_pipeline /opt/chrombpnet-upstream
USER $MAMBA_USER
# Retain the pinned source and license. Import only its attribution helpers;
# installing its old training requirements would replace the pinned TF environment.
RUN git -C /opt/chrombpnet-upstream init \
    && git -C /opt/chrombpnet-upstream remote add origin https://github.com/kundajelab/chrombpnet.git \
    && git -C /opt/chrombpnet-upstream fetch --depth 1 origin "$CHROMBPNET_COMMIT" \
    && git -C /opt/chrombpnet-upstream checkout --detach FETCH_HEAD \
    && test "$(git -C /opt/chrombpnet-upstream rev-parse HEAD)" = "$CHROMBPNET_COMMIT" \
    && rm -rf /opt/chrombpnet-upstream/.git
COPY --chown=$MAMBA_USER:$MAMBA_USER motif_pipeline/scripts /opt/motif_pipeline/scripts
RUN python /opt/motif_pipeline/scripts/expose_cuda.py \
    && python -m pip check
WORKDIR /tmp

# GitHub runs the real CPU smoke test before image export. The read-only
# BuildKit mount does not add test files or generated models to the image.
RUN --mount=type=bind,source=tests,target=/tests \
    python /tests/smoke_motif_image.py
