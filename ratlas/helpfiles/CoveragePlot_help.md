### CoveragePlot of ATAC signals for gene or genomic coordinates

You may search for ATAC-seq signals by 1) gene name by selecting the `CoveragePlot` plot option and/or 2) by genomic coordinates by selecting the `CoveragePlot_genome` option.

For either choice, you may split the tracks by metadata available in the dataset by selecting a group from `Choose all or split by group`. Due to the amount of data that may be displayed, we suggest to select specific cell clusters of interest to be visualized (by choosing clusters of interest in `CoveragePlot plot option: choose all or show specific clusters`) to limit the number of tracks shown and speed up plotting, however this is not a requirement.

A few notes linked to plotting genomic coordinates:

* __IMPORTANT__: The current dataset was mapped against the rat __rn6__ / __Rnor_6.0__ assembly from Ensembl (v.95, see archive from Ensembl [here](http://jan2019.archive.ensembl.org/Rattus_norvegicus/Info/Index)). All genomic coordinate inputs should be relative to this reference assembly and __not__ newer versions (e.g.:  mRatBN7.2).
* The format for the gene coordinates in the dataset should be `chromossome_number-start-end` (e.g: `1-172934-175664`), however an input that includes `:` by the chromosome name and commas in the start/end positions (`1:172,934-175,664`) is also acceptable. In these cases the app will automatically convert the input to the expected input format and you should see a message at the bottom right corner indicating this automatic change.
* Adding `chr` to the chromossome number is also acceptable and the app will automatically remove it as well (e.g: `chr1-172934-175664`  or `chr1:172,934-175,664`) for compatibility with the annotation format of the datasets.
* Other input formats that are not outlined above are not valid inputs.