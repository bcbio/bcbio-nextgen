class: CommandLineTool

description: "Parse a bcbio system file to provide a system configuration file and directories"

requirements:
  - class: CreateFileRequirement
    fileDef:
      - filename: bcbio_system.yaml
        fileContent:
          engine: "cwl:JsonPointer"
          script: "job/system_config"

inputs:
  - id: "#system_config"
    type: File
    description:
      bcbio system YAML file.

outputs:
  - id: "#system_config_prep"
    type: File
    outputBinding:
      glob: "bcbio_system.yaml"

baseCommand: ["bcbio_nextgen.py", "runfn", "bcbioprep"]

arguments:
  - "bcbio_system.yaml"