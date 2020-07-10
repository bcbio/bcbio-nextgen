# HLA typing

`hlacaller` -- Perform identification of highly polymorphic HLAs with human build 38 (hg38).
The recommended option is `optitype`, using the [OptiType](https://github.com/FRED-2/OptiType) caller.
Also supports using the [bwa HLA typing implementation](https://github.com/lh3/bwa/blob/master/README-alt.md#hla-typing) with `bwakit`
For this to run, `aligner` must be set to `bwa`.
