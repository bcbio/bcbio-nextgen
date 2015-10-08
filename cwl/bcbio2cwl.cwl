#!/usr/bin/env cwl-runner

class: Workflow

inputs:
  - id: "#system_config"
    default: null
    type: File
    description: bcbio system configuration file. Can be null to use the default.
  - id: "#run_config"
    type: File
    description: bcbio run configuration file in YAML format.

requirements:
  - class: EnvVarRequirement
    envDef:
      - envName: "MPLCONFIGDIR"
        envValue: '.'

hints:
  - class: DockerRequirement
    dockerImport: https://s3.amazonaws.com/bcbio_nextgen/bcbio-nextgen-docker-image.gz
    dockerImageId: chapmanb/bcbio-nextgen-devel

outputs:
  - id: "#workflow"
    type: File
    source: "#create_cwl.workflow"

steps:
  - id: "#prep_system"
    run: {import: prep_system-tool.cwl}
    inputs:
      - {id: "#prep_system.run_config", source: "#run_config"}
      - {id: "#prep_system.system_config", source: "#system_config"}
    outputs:
      - {id: "#prep_system.system_config_prep"}

  - id: "#organize_samples"
    run: {import: organize_samples-tool.cwl}
    inputs:
      - {id: "#organize_samples.run_config", source: "#run_config"}
      - {id: "#organize_samples.system_config_prep", source: "#prep_system.system_config_prep"}
    outputs:
      - { id: "#organize_samples.world" }

  - id: "#create_cwl"
    run: {import: create_cwl-tool.cwl}
    inputs:
      - {id: "#create_cwl.world", source: "#organize_samples.world"}
      - {id: "#create_cwl.run_config", source: "#run_config"}
    outputs:
      - {id: "#create_cwl.workflow"}
