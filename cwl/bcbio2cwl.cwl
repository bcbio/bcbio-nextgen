#!/usr/bin/env cwl-runner

class: Workflow

inputs:
  - id: "#system_config"
    type:
      - "null"
      - File
    description: bcbio system configuration file. Can be null to use the default.

  - id: "#run_config"
    type: File
    description: bcbio run configuration file in YAML format.

hints:
  - class: DockerRequirement
    dockerImport: https://s3.amazonaws.com/bcbio_nextgen/bcbio-nextgen-docker-image.gz
    dockerImageId: chapmanb/bcbio-nextgen-devel

outputs:
  - id: "#world"
    type: File
    source: "#organize_samples.world"

steps:
  - id: "#bcbioprep"
    run: {import: bcbioprep-tool.cwl}
    inputs:
      - {id: "#bcbioprep.system_config", source: "#system_config"}
    outputs:
      - {id: "#bcbioprep.system_config_prep"}

  - id: "#organize_samples"
    run: {import: organize_samples-tool.cwl}
    inputs:
      - {id: "#organize_samples.run_config", source: "#run_config"}
      - {id: "#organize_samples.system_config_prep", source: "#bcbioprep.system_config_prep"}
    outputs:
      - { id: "#organize_samples.world" }