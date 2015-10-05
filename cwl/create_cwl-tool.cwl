class: CommandLineTool

description: "Convert sample input file and configuration into a world object with items to run."

inputs:
  - id: "#world"
    type: File
    description: World object with sample and work
    inputBinding:
      position: 0
  - id: "#run_config"
    type: File
    description: bcbio run configuration YAML file
    inputBinding:
      position: 1

outputs:
  - id: "#workflow"
    type: File
    outputBinding:
      glob: "*-workflow/*.cwl"

baseCommand: ["/usr/local/bin/bcbio_nextgen.py", "runfn", "create_cwl"]
