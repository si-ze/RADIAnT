Configuration files for RADIAnT. 

Pass the config file corresponding to your sequencing reads to the snakemake pipeline the following way (e.g. for RADICL-seq): 

```
snakemake -s /path/to/RADIAnT/workflow/RADIAnT.smk --configfile=/path/to/RADIAnT/config/config_RADICL_mESCs.yaml --cores 32
```

Make sure to **replace** ```/path/to/RADIAnT/``` with the **absolute path of the local directory** you have cloned the repository to. **Both** in 
* the **command** issued to run snakemake and
* in **all paths referenced in the config file**.
