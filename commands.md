```
Nedeed: # Graphviz 
brew install graphviz 

########

# Emits Graphviz DOT text
snakemake --dag > current_pipeline.dot

# Render
dot -Tsvg current_pipeline.dot > current_pipeline.svg
```


---

Epigenomics
Trimming after first QC:

```
cd osd47-spaceflight-multiomics/epigenomics
snakemake -p --cores 4 --runtime-source-cache-path .snakemake-cache/runtime-source-cache-20260615b
```

only trim:

```snakemake -p --cores 4 results/trimmed/.trimmed.ok --runtime-source-cache-path .snakemake-cache/runtime-source-cache-20260615b```