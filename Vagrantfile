# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "bento/ubuntu-18.04"

  config.vm.hostname = "bcbio"

  config.vm.provider :virtualbox do |vb|
    vb.memory = 4096
    vb.cpus = 2
  end

  # to make any additional data on the host available inside the VM
  # (for example: reference genomes, pipeline inputs, etc)
  # set BCBIO_DATA_DIR environment variable on the host to a directory that contains the data
  if ENV["BCBIO_DATA_DIR"]
    config.vm.synced_folder ENV["BCBIO_DATA_DIR"], "/data"
  end

  config.vm.provision :shell, :path => "scripts/vagrant.sh"

end
