
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 

# Install the development version of catGenes from GitHub, 
# together with its required dependencies 
BiocManager::install("DBOSlab/catGenes", dependencies = TRUE)

library(catGenes)

### Buscando sequencias do genebank do gênero Anthurium e do gene específico
### Salva num arquivo .fasta

# mineTaxa usa qualquer termo que poderia ser buscado no NCBI
# pode usar número de acesso tb, com "OR" na busca booleana
# pode fazer um loop pra ir baixando cada sequencia dos números de acessos de um vetor

# Exemplo 1
result <- mineTaxa(
  term = "Anthurium[Organism] AND rbcL[Gene]",
  filename = "Anthurium_rbcL.fasta",
  clean.taxa = TRUE,
  add.voucher = TRUE,          # Don't add voucher info
  original.query = TRUE,        # Keep original unprocessed file
  plastome.apart = TRUE,       # Don't separate plastomes
  rm.duplicated = TRUE,        # Keep all sequences
  retmax = 1000,
  verbose = TRUE,
  save = TRUE,
  dir = "RESULTS_mineTaxa"
)

# Exemplo 2
result.matK <- mineTaxa(
  term = "Anthurium[Organism] AND matK[Gene]",
  filename = "Anthurium_matK.fasta",
  clean.taxa = TRUE,
  add.voucher = TRUE,          # Don't add voucher info
  original.query = TRUE,        # Keep original unprocessed file
  plastome.apart = TRUE,       # Don't separate plastomes
  rm.duplicated = TRUE,        # Keep all sequences
  retmax = 1000,
  verbose = TRUE,
  save = TRUE,
  dir = "RESULTS_mineTaxa"
)

# Exemplo 3 - outgroup
mineTaxa(
  term = "Pothos scandens[Organism] AND rbcL[Gene]",
  filename = "Pothos_rbcL.fasta",
  clean.taxa = TRUE,
  add.voucher = TRUE,          # Don't add voucher info
  original.query = TRUE,        # Keep original unprocessed file
  plastome.apart = TRUE,       # Don't separate plastomes
  rm.duplicated = TRUE,        # Keep all sequences
  retmax = 1000,
  verbose = TRUE,
  save = TRUE,
  dir = "RESULTS_mineTaxa"
)

mineTaxa(
  term = "Pothos scandens[Organism] AND matK[Gene]",
  filename = "Pothos_matK.fasta",
  clean.taxa = TRUE,
  add.voucher = TRUE,          # Don't add voucher info
  original.query = TRUE,        # Keep original unprocessed file
  plastome.apart = TRUE,       # Don't separate plastomes
  rm.duplicated = TRUE,        # Keep all sequences
  retmax = 1000,
  verbose = TRUE,
  save = TRUE,
  dir = "RESULTS_mineTaxa"
)


# Tabela criada pelo Domingos com o grupo que ele estuda
#Planilha com dados já curados, aí usa mineSeq
data(GenBank_accessions)
mineSeq(inputdf = GenBank_accessions,
        gb.colnames = c("ETS", "ITS", "matK", "petBpetD", "trnTF", "Xdh"),
        as.character = FALSE,
        verbose = TRUE,
        save = TRUE,
        filename = "GenBanK_seqs",
        dir = "RESULTS_mineSeq")

## Combinar os arquivos fasta
combineFASTA(
   input.files = c("RESULTS_mineTaxa/Anthurium/Anthurium_matK.fasta", "RESULTS_mineTaxa/Pothos/Pothos_matK.fasta"),
   output.file = "combined.matK.fasta"
 )

combineFASTA(
  input.files = c("RESULTS_mineTaxa/Anthurium/Anthurium_rbcL.fasta", "RESULTS_mineTaxa/Pothos/Pothos_rbcL.fasta"),
  output.file = "combined.rbcL.fasta"
)

## Alinhar sequencias
alignSeqs(filepath = "RESULTS_combineFASTA/03Sep2026/",
          method = "ClustalW",
          gapOpening = "default",
          format = "NEXUS",
          verbose = TRUE,
          dir = "RESULTS_alignSeqs")

## Exemplo catmultGenes - concatenate

data(Luetzelburgia)
catdf <- catmultGenes(Luetzelburgia,
                      maxspp = TRUE,
                      shortaxlabel = TRUE,
                      missdata = TRUE,
                      verbose = TRUE)

outgrouptaxa <- c("Vataireopsis_araroba", "Vataireopsis_speciosa")
catdf <- catmultGenes(Luetzelburgia,
                      maxspp = FALSE,
                      shortaxlabel = TRUE,
                      missdata = FALSE,
                      outgroup = outgrouptaxa,
                      verbose = TRUE)

