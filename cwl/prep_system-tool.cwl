class: CommandLineTool

description: "Parse a bcbio system file to provide a system configuration file and directories"

inputs:
  - id: "#run_config"
    type: File
    description: bcbio run information YAML file
      bcbio system YAML file.
    inputBinding:
      position: 0
  - id: "#system_config"
    default: null
    type: File
    description:
      bcbio system YAML file.
    inputBinding:
      position: 1

outputs:
  - id: "#system_config_prep"
    type: File
    outputBinding:
      glob: "system_config_prep.yaml"

baseCommand: ["bcbio_nextgen.py", "runfn", "prep_system", "--raw", "--out", "system_config_prep.yaml"]