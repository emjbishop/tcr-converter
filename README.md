# tcr-converter

Convert TCRrep data from 10X, Adaptive, or AIRR format to IMGT. Most of this code was borrowed from [tcrdist3](https://github.com/kmayerb/tcrdist3)
  
**Adaptive/AIRR**: Can only take single-chain data and won't output all columns by default. Specify extra columns with `--extra`

## Usage

```
tcr_converter.py [-h HELP] -i INPUT -c CHAINS -o OUTPUT [-f INPUT FORMAT] [-e EXTRA COLUMNS] [-s SPECIES]
```

**Example:**
```
tcr_converter.py -i filtered_contig_annotations.csv -c 'alpha-beta' -o out.tsv -e ['ptid', 'cohort']
```

## Requirements

* pandas
* numpy
