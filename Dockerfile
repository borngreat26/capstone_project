FROM mambaorg/micromamba:2.0.5
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -n capstone -f /tmp/environment.yml && micromamba clean --all --yes
ENV ENV_NAME=capstone MPLBACKEND=Agg MPLCONFIGDIR=/tmp/matplotlib
WORKDIR /workspace
# Explicit allowlist: full datasets, local executables and private notebooks are never copied.
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ src/
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ scripts/
COPY --chown=$MAMBA_USER:$MAMBA_USER tests/ tests/
COPY --chown=$MAMBA_USER:$MAMBA_USER data/demo/ data/demo/
CMD ["python", "scripts/run_demo.py", "--output", "/tmp/demo"]
