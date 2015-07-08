class: CommandLineTool

description: "Convert sample input file and configuration into a world object with items to run."

inputs:
  - id: "#system_config_prep"
    type: File
    description: bcbio system YAML file and directories, prepped.
    InputBinding:
      position: 1
  - id: "#run_config"
    type: File
    description: bcbio run configuration YAML file
    InputBinding:
      position: 2

outputs:
  - id: "#world"
    type: File
    outputBinding:
      glob: "*organize_samples*.yaml"

baseCommand: ["bcbio_vm.py", "runfn", "organize_samples"]