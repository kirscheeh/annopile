#!/usr/bin/env python

import sys 
import pandas as pd

DETAILS = pd.read_excel(sys.argv[1], sheet_name=["bed_file", "coverage", "snv", "indel", "details"])
ANNOVAR=  pd.read_csv(sys.argv[2])

ANNOVAR.drop(['Otherinfo2','Otherinfo3', 'Otherinfo4', 'Otherinfo5', 'Otherinfo6', 'Otherinfo7','Otherinfo8', 'Otherinfo9', 'Otherinfo10', 'Otherinfo11', 'Otherinfo12','Otherinfo13', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG',"cosmic68", "CLNALLELEID"], inplace=True, axis=1)

# filter synonymous
ANNOVAR = ANNOVAR[ANNOVAR["ExonicFunc.refGene"] != "synonymous SNV"]
# filter intronic
ANNOVAR = ANNOVAR[ANNOVAR["Func.refGene"] != "intronic"]

writer = pd.ExcelWriter(sys.argv[1], engine="xlsxwriter")

# add annovar sheet
for name, df in DETAILS.items():
    df.to_excel(writer, index=False, sheet_name=name)
    if name == "snv":
        ANNOVAR.to_excel(writer, index=False, sheet_name="snv_annovar")

writer.close()
