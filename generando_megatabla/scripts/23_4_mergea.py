import polars as pl
import glob
# Última ejecución le costó 34 mins
print("Lazy scan of our master table")
master_df = pl.scan_csv("megatablas/megatabla_8.tsv", separator="\t")
print("Lazy scan of kmer tables...")

parquet_files = glob.glob("pivotadas/*.parquet")

# for each parquet file
for file_path in parquet_files:
    kmer_lazy = pl.scan_parquet(file_path).with_columns(
        pl.col("^kmer_.*$").cast(pl.Float32)
    )
    master_df = master_df.join(kmer_lazy,left_on="associated_transcript", right_on="SeqID", how="left")

print("Execution graph ready. Starting with streaming...")

# 5. Guardar el resultado de forma eficiente
output_file = "megatablas/megatabla_9/*"
print(f"Saving result in  {output_file}...")
master_df.sink_parquet(output_file, compression="zstd", maintain_order = False)

print("End of script")
