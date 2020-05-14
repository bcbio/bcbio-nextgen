# Disambiguation

* `disambiguate` For mixed or explant samples, provide a list of `genome_build`
identifiers to check and remove from alignment.
Currently supports cleaning a single organism. For example, with `genome_build: hg19`
and `disambiguate: [mm10]`, it will align to hg19 and mm10, run disambiguation and discard
reads confidently aligned to mm10 and not hg19.
Affects fusion detection when `star` is chosen as the aligner. Aligner must be set to a non false value for this to run.
