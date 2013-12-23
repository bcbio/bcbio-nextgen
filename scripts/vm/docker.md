Notes on using docker for bcbio-nextgen installation.

## General

Install docker http://docs.docker.io/en/latest/installation/ubuntulinux/

    sudo sh -c "curl https://get.docker.io/gpg | apt-key add -"
    sudo sh -c "echo deb http://get.docker.io/ubuntu docker main > /etc/apt/sources.list.d/docker.list"
    sudo apt-get update
    sudo apt-get install lxc-docker

Setup docker group: http://docs.docker.io/en/latest/use/basics/#dockergroup
Have to do `newgrp docker` to refresh groups until log in and out.

    sudo groupadd docker
    sudo gpasswd -a ${USERNAME} docker
    sudo service docker restart
    newgrp docker

## Creating images from scratch

Start up docker

    DID=$(docker run -d -i -t -p 8085:8085 stackbrew/ubuntu:13.10 /bin/bash)
    docker attach $DID

Install bcbio-nextgen via instructions in Dockerfile. Then commit:

    docker commit $DID chapmanb/bcbio-nextgen-devel

or build directly:

    docker build -t chapmanb/bcbio-nextgen-devel .

## Update images to index

    DID=$(docker run -d -i -t -p 8085:8085 chapmanb/bcbio-nextgen-devel /bin/bash)
    docker attach $DID
    docker commit $DID chapmanb/bcbio-nextgen-devel
    docker push chapmanb/bcbio-nextgen-devel
