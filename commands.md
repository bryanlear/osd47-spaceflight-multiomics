```
Nedeed: # Graphviz 
brew install graphviz 

########

# Emits Graphviz DOT text
snakemake --dag > current_pipeline.dot

# Render
dot -Tsvg current_pipeline.dot > current_pipeline.svg
```