# Use official mambaforge image
FROM condaforge/mambaforge:23.3.1-0

# Avoid interactive APT prompts (e.g. tzdata)
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install retry script with exponential backoff
RUN printf '#!/bin/sh\n\
for i in 1 2 3 4 5; do\n\
  echo "Attempt $i: $@"\n\
  $@ && exit 0\n\
  sleep $((i*2))\n\
done\n\
echo "Failed after 5 attempts: $@"\n\
exit 1' > /usr/local/bin/retry && \
    chmod +x /usr/local/bin/retry

# 2. Securely update APT and install system keyring + tzdata preconfig
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gnupg \
        ca-certificates \
        debian-archive-keyring \
        tzdata && \
    ln -fs /usr/share/zoneinfo/Etc/UTC /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    rm -rf /var/lib/apt/lists/*

# 3. Install system dependencies
RUN retry apt-get update -qy && \
    retry apt-get install -y --no-install-recommends \
        software-properties-common \
        procps \
        libgl1-mesa-glx \
        libglib2.0-0 \
        build-essential \
        zlib1g-dev \
        libssl-dev \
        openssl \
        locales \
        gfortran \
        libgfortran5 && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8

# 4. Copy and create environment
COPY environment.yml .
RUN mamba env create -f environment.yml && \
    mamba clean -afy

# 5. Install Perl modules
RUN conda run -n rmap-1.0 cpanm -n \
    LWP::Protocol::https \
    LWP::Simple \
    Module::Build && \
    conda clean -afy

# 6. Set up working directory and permissions
WORKDIR /app
COPY . /app
RUN chmod +x /app/rMAP && chmod -R a+rX /app

# 7. Runtime configuration
ARG CONDA_ENV=rmap-1.0
ENV PERL5LIB="/opt/conda/envs/${CONDA_ENV}/lib/perl5/site_perl/5.32.1:/opt/conda/envs/${CONDA_ENV}/lib/perl5/5.32.1" \
    PATH="/opt/conda/envs/${CONDA_ENV}/bin:$PATH" \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# 8. Configure matplotlib
RUN mkdir -p ${MPLCONFIGDIR} && chmod -R 777 ${MPLCONFIGDIR}

EXPOSE 8080

# 🟢 Run the Bash pipeline script
ENTRYPOINT ["/bin/bash", "/app/rMAP"]