# data/

This directory stores raw and intermediate data for the 16S amplicon sequencing analysis.

## Expected Structure

```
data/
├── raw/                    # Raw FASTQ files (paired-end, gzipped)
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── ...
├── trimmed/                # Primer-trimmed FASTQ files (output of cutadapt)
├── figaro/                 # rFigaro output (trimParameters.json, EE plots)
├── dada2/                  # DADA2 intermediate objects
│   ├── seqtab.rds
│   ├── seqtab_nochim.rds
│   ├── taxa.rds
│   └── track_reads.csv
├── phyloseq/               # Assembled phyloseq objects
│   └── ps.rds
└── metadata/               # Sample metadata
    └── metadata.csv        # Required: sample IDs, group, collection info
```

## Metadata format

The `metadata/metadata.csv` file should contain at minimum:

| Column       | Description                                      |
|-------------|--------------------------------------------------|
| SampleID    | Unique sample identifier (matches FASTQ filename)|
| Group       | Experimental group / condition                   |
| SampleType  | "Sample", "NegControl", or "PosControl"          |

Additional columns (e.g., host variables, batch, collection date) are welcome.

## Notes

- Raw FASTQ files are **not** tracked by git (add to `.gitignore`).
- Large RDS objects are **not** tracked by git.
- Only `metadata.csv` and small summary tables should be committed.
