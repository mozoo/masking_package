# masking_package
A custom, simple suite of bash and R scripts to remove phylogenetic noise using different softwares.

Multiple sequence alignment is a key step in a wide range of biological studies.
Many different approaches and softwares have been developed to this purpose; however,
despite the striking improvement in alignment accuracy in recent years, regions of
uncertainty may still hamper a correct retrieval of positional homologies.

Several metrics were suggested to estimate overall alignment quality. Measuring
site-wise alignment accuracy (i.e. at each individual position), however, is often
mandatory, especially in a phylogenetic context, where poorly aligned sites may
contribute noise to the tree inference. Those poorly aligned sites should be detected
and discarded: this task is typically known as trimming or masking.

Masking can be carried out through many different methods and tools: the present package,
which is called masking_package, was written to call five different masking tools,
compare their results and produce a final consensus accounting for all of them.

masking_package was written in bash and R languages to perform and compare masking
from different softwares.

It is also possible to provide separate alignments (e.g., separate genes) in the same
analysis: they will be separately aligned, masked, and evaluated; then, they will be
concatenated in the final resulting dataset.
