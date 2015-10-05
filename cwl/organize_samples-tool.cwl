class: CommandLineTool

description: "Convert sample input file and configuration into a world object with items to run."

inputs:
  - id: "#system_config_prep"
    type: File
    description: bcbio system YAML file and directories, prepped.
    inputBinding:
      position: 0
  - id: "#run_config"
    type: File
    description: bcbio run configuration YAML file

outputs:
  - id: "#world"
    type: File
    outputBinding:
      glob: "world.yaml"

baseCommand: ["/usr/local/bin/bcbio_nextgen.py", "runfn", "organize_samples", "--out", "world.yaml"]
