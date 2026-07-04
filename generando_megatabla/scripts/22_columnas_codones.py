import pandas as pd
megatabla = pd.read_csv("megatablas/megatabla_7.tsv", sep="\t")
frecs = pd.read_csv("frec_relativa_codones.tsv", sep="\t")
nueva = megatabla.merge(right=frecs, left_on="associated_transcript", right_on="transcript_id", how="left")
nueva.to_csv("megatablas/megatabla_8.tsv", sep="\t",index=False)
