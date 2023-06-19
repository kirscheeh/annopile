#!/usr/bin/env python

import sys, glob 
import pandas as pd 

RESULT_PATH=sys.argv[1]

samples = pd.read_excel(sys.argv[2])

samples["library"] = "Agilent_"+samples["Library ID"]
samples["sample"] = samples["Tumor ID"]

samples = samples[["sample", "library"]]

with open("data/SOI.csv", "w") as SOI:

    for index, (sample, lib) in samples.iterrows():
        pa = glob.glob(f"{RESULT_PATH}/{sample}-DNA*/{lib}/bam_v2.0/{sample}*MERGED.bam")[0]
        SOI.write(f"{pa}\n")