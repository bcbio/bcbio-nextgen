FROM stackbrew/ubuntu:13.10
MAINTAINER Brad Chapman "https://github.com/chapmanb"

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
RUN python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local --toolplus data --nodata --nosudo -u development
RUN bcbio_nextgen.py upgrade --isolate
RUN mkdir -p /mnt/biodata
RUN mv /usr/local/share/bcbio-nextgen/galaxy/bcbio_system.yaml /usr/local/share/bcbio-nextgen/config
RUN rmdir /usr/local/share/bcbio-nextgen/galaxy
RUN ln -s /mnt/biodata/galaxy /usr/local/share/bcbio-nextgen/galaxy
RUN ln -s /mnt/biodata/gemini_data /usr/local/share/bcbio-nextgen/gemini_data
RUN ln -s /mnt/biodata/genomes /usr/local/share/bcbio-nextgen/genomes
RUN ln -s /mnt/biodata/liftOver /usr/local/share/bcbio-nextgen/liftOver
RUN chmod a+rwx /usr/local/share/bcbio-nextgen
RUN chmod a+rwx /usr/local/share/bcbio-nextgen/config
RUN chmod a+rwx /usr/local/share/bcbio-nextgen/config/*.yaml
RUN chmod a+rwx /usr/local/share/bcbio-nextgen/gemini-config.yaml
RUN mkdir -p /tmp/bcbio-nextgen
RUN echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> /etc/bash.bashrc
RUN echo 'export PERL5LIB=/usr/local/lib/perl5:/usr/local/lib/perl5/site_perl:${PERL5LIB}' >> /etc/bash.bashrc

# Clean up space
RUN rm -rf $(brew --cache)
RUN rm -rf /.cpanm

EXPOSE 8085
ENTRYPOINT ["/usr/local/bin/bcbio_nextgen.py", "server", "--port", "8085", "--biodata_dir", "/mnt/biodata"]