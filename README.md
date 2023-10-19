The following corrections were added to the README file:
### Fastq files filtering function
#### run_fastq_filter
This function calculates sequences GC content, length, and quality and filters fastq sequences according to counted values. As input, it takes the following parameters:
* **input_path** - path to input fastq file. For example, `input_path = 'example_fastq.fastq'`
* **gc_bounds** - boundaries for filtering according to sequence GC content. Two numbers indicate the lower and upper limits (inclusively) of sequence filtering. If the user specifies only one value to filter, it will be considered the upper bound.
* **length_bounds** - boundaries for filtering according to sequence length. Default - (0, 2**32). Two numbers indicate the lower and upper limits (inclusively) of sequence filtering. If the user specifies only one value to filter, it will be considered the upper bound.
* **quality_threshold** - boundaries for filtering according to sequence quality. Default - 0.
* **output_filename** - the name of the output file. Write 'None' if you don't want to specify filename. The output fastq file is written in the separate folder `fastq_filtrator_resuls`


