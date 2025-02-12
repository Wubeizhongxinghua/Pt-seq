# Pt-seq analysis pipeline

Pt seq is a high-throughput library construction method specifically designed to detect high-sensitivity single base resolution Pt drag binding sites on the genome in various Pt drag treated biological samples.

## Preperation

1. Clone the repository

Clone this repository into your service and `cd` into the directory.

After clonning, your directory structure is like this:
```
.
├── 0.process_platin_q0_unique.snake.py
├── 1.from_q0_to_q1to10_and_40.snake.py
├── 2.downstream_pipeline
│   ├── main_downs.py
│   └── snk_script
│       ├── ...
├── cal_covered_genome.py
├── filter_reads_withSbfI.py
├── README.md
└── tools
    └── 0.all_signal.py
```

2. Prepare the raw data

Inside the repository direction, create a directory called `data`, and put your Pt-seq `.fastq.gz` into the `data/{project_name}/` directory. In this example, we name the project `example_data`. After that, record the name of samples into `sample_id/{project_name}.txt` file.

Directory structure:

```
.
├── 0.process_platin_q0_unique.snake.py
├── 1.from_q0_to_q1to10_and_40.snake.py
├── 2.downstream_pipeline
│   ├── main_downs.py
│   └── snk_script
│       ├── ...
├── cal_covered_genome.py
├── data
│   └── example_data
│       ├── sample1_1.fq.gz
│       ├── sample1_2.fq.gz
│       ├── sample2_1.fq.gz
│       └── sample2_2.fq.gz
├── sample_id
│   └── example_data.txt
├── ...
```

`example_data.txt`:
```
sample1
sample2
```

**NOTE: The sample name should not contain character underline "\_" !!**

3. Set the reference genome

Edit the main pipeline file `0.process_platin_q0_unique.snake.py`

You only need to modify these variables based on your needs:
- `genome_dir`: The directory (absolute address recommended) where the reference genome `fasta` file and the relavent bowtie2 index locate.
- `genome_fasta_file`: Reference genome `fasta` file basename.
- `genome_index_basename`: Prefix of bowtie2 index.
- `datasets`: This is a list, you need to add your project name inside. In this example, we modify the datasets to `datasets = ['example_data']`

If you have other more customized needs, just modify this file for free.


