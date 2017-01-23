FROM stackbrew/ubuntu:14.04
MAINTAINER Brad Chapman "https://github.com/chapmanb"

# v1.0.0a -- https://github.com/chapmanb/bcbio-nextgen/commit/33a9ed5

# Setup a base system 
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 xz-utils pigz && \
# Support inclusion in Arvados pipelines
    apt-get install -y --no-install-recommends libcurl4-gnutls-dev mbuffer python2.7-dev python-virtualenv && \
# Not added by default to save space
#   apt-get install -y libglu1-mesa && \

# bcbio-nextgen installation
    mkdir -p /tmp/bcbio-nextgen-install && cd /tmp/bcbio-nextgen-install && \
    wget --no-check-certificate \
      https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py && \
    python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen \
      --isolate --nodata -u development --tooldir=/usr/local && \
    git config --global url.https://github.com/.insteadOf git://github.com/ && \
    /usr/local/share/bcbio-nextgen/anaconda/bin/bcbio_nextgen.py upgrade --isolate --tooldir=/usr/local --tools && \
    /usr/local/share/bcbio-nextgen/anaconda/bin/bcbio_nextgen.py upgrade --isolate -u development --tools && \

# setup paths
    echo 'export PATH=/usr/local/bin:$PATH' >> /etc/profile.d/bcbio.sh && \

# add user run script
    wget --no-check-certificate -O createsetuser \
      https://raw.github.com/chapmanb/bcbio-nextgen-vm/master/scripts/createsetuser && \
    chmod a+x createsetuser && mv createsetuser /sbin && \

# clean filesystem
    cd /usr/local && \ 
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/tmp/* && \
    /usr/local/share/bcbio-nextgen/anaconda/bin/conda clean --yes --tarballs && \
    rm -rf /usr/local/share/bcbio-nextgen/anaconda/pkgs/qt* && \
    rm -rf /usr/local/.git && \
    rm -rf /.cpanm && \
    rm -rf /tmp/bcbio-nextgen-install && \

# Create directories and symlinks for data 
    mkdir -p /mnt/biodata && \
    mkdir -p /tmp/bcbio-nextgen && \
    mv /usr/local/share/bcbio-nextgen/galaxy/bcbio_system.yaml /usr/local/share/bcbio-nextgen/config && \
    rmdir /usr/local/share/bcbio-nextgen/galaxy && \
    ln -s /mnt/biodata/galaxy /usr/local/share/bcbio-nextgen/galaxy && \
    ln -s /mnt/biodata/gemini_data /usr/local/share/bcbio-nextgen/gemini_data && \
    ln -s /mnt/biodata/genomes /usr/local/share/bcbio-nextgen/genomes && \
    ln -s /mnt/biodata/liftOver /usr/local/share/bcbio-nextgen/liftOver && \
    chmod a+rwx /usr/local/share/bcbio-nextgen && \
    chmod a+rwx /usr/local/share/bcbio-nextgen/config && \
    chmod a+rwx /usr/local/share/bcbio-nextgen/config/*.yaml && \

# Ensure permissions are set for update in place by arbitrary users
    find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
    find /usr/local -perm /u+w -execdir chmod a+w {} \;
