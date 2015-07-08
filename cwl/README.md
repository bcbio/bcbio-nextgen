## bcbio with the Common Workflow Language (CWL)

In progress work to run bcbio using a tool that supports the
[Common Workflow Language (CWL)][0]. The initial implementation parses a bcbio
system YAML file and run info YAML file to produce a CWL workflow and run that
handles multiple samples.

Usage;
```
cwl-runner bcbio2cwl.cwl bcbio2cwl-testinput.json
```

[0]: https://github.com/common-workflow-language/common-workflow-language
