# tcr-converter

## **NOTE: This tool is under development and is not guaranteed to work at this time.**  
  
Convert TCRrep data from Adaptive v4 or AIRR format to IMGT. Most of this code was 
borrowed from [tcrdist3](https://github.com/kmayerb/tcrdist3)
  
Can only take single-chain data and won't output all columns by 
default. Specify extra columns with `--extra`

## Usage  
```
tcr_converter.py [-h HELP] -i INPUT [-c CHAIN] [-o OUTPUT] [-e EXTRA COLUMNS] [-s SPECIES]

optional arguments:
  -h, --help       show this help message and exit
  -i , --input     Input TCR filepath
  -c , --chain     Input TCR chain. Options are: "alpha", "beta", "gamma", "delta". Default is "beta"
  -o , --output    Output TCR filepath. Extension should be .csv
  -e , --extra     List of extra input Adaptive/AIRR columns to keep in addition to the CDR3/V/J ones
  -s , --species   Options are: "human", "mouse". Default is "human"
```

### Example  
  
Beta chain data with user-added 'ptid' and 'cohort' columns. Output to be saved 
in current working directory.
```
tcr_converter.py -i adaptive_tcrs.tsv -c 'beta' -e ['ptid', 'cohort']
```

## Requirements

* pandas
* numpy
* tcrdist3
