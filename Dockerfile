FROM debian:stretch-slim
LABEL maintainer="BcBio Maintainers (https://github.com/bcbio/bcbio-nextgen)"

# Setup a base system
RUN apt-get update && \
    apt-get install -y wget python git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/tmp/*

# bcbio-nextgen installation
RUN mkdir -p /tmp/bcbio-nextgen-install && \
    cd /tmp/bcbio-nextgen-install && \
    wget --no-check-certificate \
      https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py && \
    python bcbio_nextgen_install.py \
      /usr/local/share/bcbio-nextgen \
      --tooldir=/usr/local \
      --isolate \
      --minimize-disk \
      --nodata \
      -u development && \
# setup paths and cleanup
    echo 'export PATH=/usr/local/bin:$PATH' >> /etc/profile.d/bcbio.sh && \
    /usr/local/share/bcbio-nextgen/anaconda/bin/conda clean --yes --all && \
    rm -rf /tmp/bcbio-nextgen-install

# add user run script
RUN wget --no-check-certificate -O createsetuser \
      https://raw.githubusercontent.com/bcbio/bcbio-nextgen-vm/master/scripts/createsetuser && \
    chmod a+x createsetuser && mv createsetuser /sbin

# Create directories and symlinks for data
RUN mkdir -p /mnt/biodata && \
    mv /usr/local/share/bcbio-nextgen/galaxy/bcbio_system.yaml /usr/local/share/bcbio-nextgen/config && \
    rmdir /usr/local/share/bcbio-nextgen/galaxy && \
    ln -s /mnt/biodata/galaxy /usr/local/share/bcbio-nextgen/galaxy && \
    ln -s /mnt/biodata/gemini_data /usr/local/share/bcbio-nextgen/gemini_data && \
    ln -s /mnt/biodata/genomes /usr/local/share/bcbio-nextgen/genomes && \
    ln -s /mnt/biodata/liftOver /usr/local/share/bcbio-nextgen/liftOver && \
    chmod -R a+rwX /usr/local/share/bcbio-nextgen && \
# Ensure permissions are set for update in place by arbitrary users
    find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
    find /usr/local -perm /u+w -execdir chmod a+w {} \;
