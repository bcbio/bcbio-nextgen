FROM stackbrew/ubuntu:13.10
MAINTAINER Brad Chapman

# Setup a base system 
RUN mkdir -p /tmp/bcbio-nextgen-install
RUN cd /tmp/bcbio-nextgen-install
RUN apt-get update
RUN apt-get install -y build-essential zlib1g-dev wget python-setuptools git

# Fake a fuse install; openjdk pulls this in 
# https://github.com/dotcloud/docker/issues/514
# https://gist.github.com/henrik-muehe/6155333
RUN apt-get install libfuse2
RUN apt-get download fuse
RUN dpkg-deb -x fuse_* .
RUN dpkg-deb -e fuse_*
RUN rm fuse_*.deb
RUN echo -en '#!/bin/bash\nexit 0\n' > DEBIAN/postinst
RUN dpkg-deb -b . /fuse.deb
RUN dpkg -i /fuse.deb

# bcbio-nextgen installation
RUN wget --no-check-certificate https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
RUN python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local --nodata --nosudo
RUN bcbio_nextgen.py upgrade --isolate
RUN mkdir -p /mnt/biodata
RUN ln -s /mnt/biodata/genomes /usr/local/share/bcbio-nextgen/genomes
RUN ln -s /mnt/biodata/liftOver /usr/local/share/bcbio-nextgen/liftOver
RUN mkdir -p /tmp/bcbio-nextgen

# gemini installation
RUN wget --no-check-certificate https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
RUN python gemini_install.py /usr/local /usr/local/share/gemini --nosudo --notools --nodata
RUN ln -s /mnt/biodata/gemini_data /usr/local/share/gemini/data
# XXX add 'data' to bcbio-nextgen toolplus defaults: /usr/local/share/bcbio-nextgen/config/install-params.yaml
    
EXPOSE 8085
ENTRYPOINT ["/usr/local/bin/bcbio_nextgen.py", "server", "--port", "8085"]