pbmpi
=====

partitioned phylobayes mpi

Options unique to the partitioned version:

`-p <partition-file>` partitioning scheme file in [PartitionFinder](http://www.robertlanfear.com/partitionfinder/) format

`-linkgam` link gamma distribution shape parameters across +G partitions (by default, each +G partition has its own shape parameter) 

`-unlinkgtr` unlink GTR exchangeabilities across GTR partitions (by default all GTR partitions share exchangeabilities)

`-ratemult` use independent partition-specific substitution rate multipliers (by default all partitions have the same relative substitution rate)

Note:

When using `-cat` or `-catgtr` with a partition scheme, any input scheme for partitioning of equilibirum frequencies is ignored, and instead profiles are models homogeneously across all sites according the the CAT mixture.
