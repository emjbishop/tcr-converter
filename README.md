# tcr-converter

Convert TCRrep data from 10X, Adaptive, or AIRR format to IMGT.

**Adaptive/AIRR**: Can only take single-chain data and won't output all columns by default. Specify extra columns with `-e`/`--extra`

Most of this code was borrowed from [tcrdist3](https://github.com/kmayerb/tcrdist3)

## Usage

```tcr_converter.py [-h HELP] -i INPUT -c CHAINS -o OUTPUT [-t INPUT FORMAT] [-e EXTRA COLUMNS] [-s SPECIES]```

## Requirements

* pandas
* numpy
