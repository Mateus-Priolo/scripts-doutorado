# Dependências do ambiente conda `methscan`

## Python

Pacotes externos usados pelos scripts Python do pacote:

- `pandas`
- `methscan`

Imports da biblioteca padrão usados: `argparse`, `csv`, `gzip`, `os`, `re`, `sys`, `pathlib`, `typing`.

## R

Pacotes exigidos pelos scripts downstream em R:

- `data.table`
- `dplyr`
- `tidyr`
- `tibble`
- `ggplot2`
- `irlba`
- `uwot`
- `igraph`

## Sistema / CLI

Obrigatórios no cenário atual:

- `bash`
- `sbatch`
- `conda`
- `python`
- `Rscript`
- `find`
- `sort`

Opcionais:

- `job-nanny`

Somente para a trilha SCNA futura:

- `samtools`
- `bedtools`

## Observação

Eu não consigo verificar daqui se o seu ambiente `methscan` já contém todos esses pacotes. O pacote atualizado apenas passa a ativar **um único ambiente**; a presença real das bibliotecas precisa ser confirmada no servidor.
