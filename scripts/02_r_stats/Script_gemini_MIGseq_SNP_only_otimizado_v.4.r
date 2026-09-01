## Script de Análise Populacional e Genômica Integrada v4.0 ##
# ==============================================================================
# CONFIGURAÇÃO INICIAL - DEFINA A ESPÉCIE
# ==============================================================================
especie_alvo <- "Avicennia"

# ==============================================================================
# 1. Carregamento de Pacotes e Definição dos Caminhos
# ==============================================================================
library(vcfR)
library(adegenet)
library(poppr)
library(ade4)
library(ggplot2)
library(pegas)
library(RColorBrewer)
library(hierfstat)
library(LEA)
library(tidyr)

pasta_especie <- file.path(".", especie_alvo)

if (!dir.exists(pasta_especie)) {
  stop(paste("Erro: A pasta", pasta_especie, "não existe no diretório atual."))
}

vcf_files <- list.files(pasta_especie, pattern = "\\.vcf$", full.names = TRUE)

if (length(vcf_files) == 0) {
  stop(paste("Nenhum arquivo .vcf foi encontrado na pasta:", pasta_especie))
} else if (length(vcf_files) > 1) {
  warning(paste("Mais de um arquivo VCF encontrado. Utilizando o primeiro:", vcf_files[1]))
}

caminho_vcf <- vcf_files[1]
cat("\n--- Processando Espécie:", especie_alvo, "---\n")

# Importação do VCF e conversão para genind
vcfR_obj <- read.vcfR(caminho_vcf, verbose = FALSE)
genind_obj <- vcfR2genind(vcfR_obj)

# Atribuição de populações
sample_names <- indNames(genind_obj)
pop_code <- sub("_.*", "", sample_names)
pop_code_factor <- as.factor(pop_code)

pop(genind_obj) <- pop_code_factor
strata(genind_obj) <- data.frame(sample = sample_names, pop = pop_code_factor)

num_pops <- nPop(genind_obj)

# Definição unificada de cores para manutenção do padrão em todos os gráficos
if (num_pops <= 9) {
  cores_pop <- brewer.pal(max(3, num_pops), "Set1")[1:num_pops]
} else {
  cores_pop <- rainbow(num_pops)
}
names(cores_pop) <- levels(pop_code_factor)

# ==============================================================================
# 2. Diversidade e Estruturação Genética (Métricas Intrínsecas e Fst)
# ==============================================================================
pop_list <- seppop(genind_obj)

tabela_diversidade <- as.data.frame(t(sapply(names(pop_list), function(p_name) {
  pop_sub <- pop_list[[p_name]]
  sum_pop <- summary(pop_sub)
  
  ho <- mean(sum_pop$Hobs, na.rm = TRUE)
  he <- mean(sum_pop$Hexp, na.rm = TRUE)
  fis <- 1 - (ho / he)
  
  c(
    N = nInd(pop_sub),
    Ho = ho,
    He = he,
    Fis = fis
  )
})))

file_div_csv <- file.path(pasta_especie, paste0("Diversidade_Genetica_", especie_alvo, ".csv"))
write.csv(tabela_diversidade, file = file_div_csv, row.names = TRUE)

hf_obj <- genind2hierfstat(genind_obj)
wc_stats <- wc(hf_obj)
fst_global <- wc_stats$FST

fst_pairwise <- pairwise.WCfst(hf_obj)
fst_pairwise <- round(fst_pairwise, 4)

file_fst_csv <- file.path(pasta_especie, paste0("Fst_Par_a_Par_", especie_alvo, ".csv"))
write.csv(fst_pairwise, file = file_fst_csv, row.names = TRUE, na = "")

# ==============================================================================
# 3. Análise de Componentes Principais (PCA Frequencista)
# ==============================================================================
X_ind <- adegenet::tab(genind_obj, NA.method = "mean")
pca_ind <- dudi.pca(X_ind, scannf = FALSE, nf = 3)

var_eixo1 <- round((pca_ind$eig[1] / sum(pca_ind$eig)) * 100, 2)
var_eixo2 <- round((pca_ind$eig[2] / sum(pca_ind$eig)) * 100, 2)

df_pca <- data.frame(
  PC1 = pca_ind$li$Axis1,
  PC2 = pca_ind$li$Axis2,
  Pop = pop(genind_obj)
)

pca_plot <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Pop, fill = Pop)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, show.legend = FALSE) +
  scale_color_manual(values = cores_pop) +
  scale_fill_manual(values = cores_pop) +
  labs(
    title = paste("PCA - Indivíduos -", especie_alvo),
    x = paste0("PC1 (", var_eixo1, "%)"),
    y = paste0("PC2 (", var_eixo2, "%)"),
    color = "Populações",
    fill  = "Populações"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(pca_plot)

file_pca_tiff <- file.path(pasta_especie, paste0("PCA_individuos2_", especie_alvo, ".tiff"))
ggsave(file_pca_tiff, plot = pca_plot, width = 400, height = 180, units = "mm", dpi = 300)

df_pca_coords <- data.frame(
  Amostra = rownames(pca_ind$li),
  PC1_SNP = pca_ind$li$Axis1,
  PC2_SNP = pca_ind$li$Axis2
)

file_pca_coords_csv <- file.path(pasta_especie, paste0("pca_snp_coordenadas_", especie_alvo, ".csv"))
write.csv(df_pca_coords, file = file_pca_coords_csv, row.names = FALSE)

# ==============================================================================
# 4. Análise Bayesiana via LEA (sNMF, Ancestria e PCoA Ajustado)
# ==============================================================================
cat("\n--- Executando Análise Bayesiana (LEA / sNMF) ---\n")

# Conversão da matriz de genótipos para o formato .geno
gt_matrix <- vcfR::extract.gt(vcfR_obj, element = "GT", as.numeric = TRUE)
gt_matrix[is.na(gt_matrix)] <- 9

file_geno <- file.path(pasta_especie, paste0(especie_alvo, "_bayes.geno"))
write.table(gt_matrix, file = file_geno, row.names = FALSE, col.names = FALSE, sep = "")

# Execução do sNMF (K = 1 a num_pops)
obj_snmf <- LEA::snmf(file_geno, K = 1:num_pops, entropy = TRUE, repetitions = 3, project = "force")
ce_values <- sapply(1:num_pops, function(k) min(cross.entropy(obj_snmf, K = k)))

# Definindo K ótimo com base na menor entropia cruzada
k_otimo <- which.min(ce_values)

# --- 4.1. Gráfico de Ancestria Individual (Estilo STRUCTURE via LEA) ---
q_matrix_lea <- LEA::Q(obj_snmf, K = k_otimo, run = 1)
colnames(q_matrix_lea) <- paste0("Cluster_", 1:k_otimo)

df_structure <- as.data.frame(q_matrix_lea)
df_structure$Individuo <- indNames(genind_obj)
df_structure$Pop <- pop(genind_obj)

df_structure_long <- pivot_longer(
  df_structure, 
  cols = starts_with("Cluster_"), 
  names_to = "Cluster", 
  values_to = "Ancestria"
)

df_structure_long$Individuo <- factor(
  df_structure_long$Individuo, 
  levels = unique(df_structure$Individuo[order(df_structure$Pop)])
)

structure_plot <- ggplot(df_structure_long, aes(x = Individuo, y = Ancestria, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(. ~ Pop, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set2", name = "Ancestria") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(
    title = paste("Ancestria Individual (sNMF / LEA - K =", k_otimo, ") -", especie_alvo),
    x = "Indivíduos por População de Origem",
    y = "Proporção de Ancestria (Q)"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14)
  )

print(structure_plot)

file_struct_tiff <- file.path(pasta_especie, paste0("Structure_sNMF_K", k_otimo, "_", especie_alvo, ".tiff"))
ggsave(file_struct_tiff, plot = structure_plot, width = 320, height = 150, units = "mm", dpi = 300)

# Exporta a matriz Q de ancestralidade em CSV
file_q_csv <- file.path(pasta_especie, paste0("Ancestria_Bayesiana_Q_K", k_otimo, "_", especie_alvo, ".csv"))
write.csv(df_structure, file = file_q_csv, row.names = FALSE)

# --- 4.2. PCoA Bayesiano Ajustado para Endogamia ---
X_raw <- adegenet::tab(genind_obj)
X_bayes <- X_raw
for (j in 1:ncol(X_bayes)) {
  na_idx <- is.na(X_bayes[, j])
  if (any(na_idx)) X_bayes[na_idx, j] <- mean(X_bayes[!na_idx, j])
}

dist_bayes <- dist(X_bayes)
pcoa_res <- cmdscale(dist_bayes, k = 3, eig = TRUE)

var_bayes1 <- round((pcoa_res$eig[1] / sum(pcoa_res$eig[pcoa_res$eig > 0])) * 100, 2)
var_bayes2 <- round((pcoa_res$eig[2] / sum(pcoa_res$eig[pcoa_res$eig > 0])) * 100, 2)

df_bayes <- data.frame(
  Axis1 = pcoa_res$points[, 1],
  Axis2 = pcoa_res$points[, 2],
  Pop   = pop(genind_obj)
)

bayes_plot <- ggplot(df_bayes, aes(x = Axis1, y = Axis2, color = Pop, fill = Pop)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, show.legend = FALSE) +
  scale_color_manual(values = cores_pop) +
  scale_fill_manual(values = cores_pop) +
  labs(
    title = paste("Ordenação Bayesiana / PCoA (Ajustada para Endogamia) -", especie_alvo),
    x = paste0("Eixo 1 (", var_bayes1, "%)"),
    y = paste0("Eixo 2 (", var_bayes2, "%)"),
    color = "Populações",
    fill  = "Populações"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(bayes_plot)

file_bayes_tiff <- file.path(pasta_especie, paste0("Ordenacao_Bayesiana_PCoA_", especie_alvo, ".tiff"))
ggsave(file_bayes_tiff, plot = bayes_plot, width = 400, height = 180, units = "mm", dpi = 300)

file_bayes_coords_csv <- file.path(pasta_especie, paste0("ordenacao_bayesiana_coordenadas_", especie_alvo, ".csv"))
write.csv(df_bayes, file = file_bayes_coords_csv, row.names = TRUE)

# ==============================================================================
# 5. Análise de Atribuição de Origem (Assignment via DAPC)
# ==============================================================================
dapc_optim <- xvalDapc(
  X_ind, 
  pop(genind_obj), 
  n.pca.max = min(50, ncol(X_ind)), 
  training.set = 0.9, 
  result = "groupMean", 
  center = TRUE,
  xval.plot = FALSE
)

best_pc <- as.numeric(dapc_optim$DAPC$n.pca)
dapc_res <- dapc(genind_obj, n.pca = best_pc, n.da = num_pops - 1)

dapc_sumario <- summary(dapc_res)
assign_rate  <- dapc_sumario$assign.per.pop
assign_pct   <- round(assign_rate * 100, 1)

ordem_pops <- levels(pop(genind_obj))
mat_post   <- dapc_res$posterior
n_ind      <- nrow(mat_post)

rotulos_facetas <- paste0(ordem_pops, "\n(", assign_pct[ordem_pops], "%)")
names(rotulos_facetas) <- ordem_pops

df_dapc_long <- data.frame(
  Individuo     = factor(rep(rownames(mat_post), times = num_pops), levels = unique(rownames(mat_post))),
  Grupo         = factor(rep(colnames(mat_post), each = n_ind), levels = ordem_pops),
  Probabilidade = as.vector(mat_post),
  Pop_Origem    = factor(rep(pop(genind_obj), times = num_pops), levels = ordem_pops)
)

df_dapc_long$Pop_Origem_Facet <- factor(
  rotulos_facetas[as.character(df_dapc_long$Pop_Origem)],
  levels = rotulos_facetas
)

dapc_plot <- ggplot(df_dapc_long, aes(x = Individuo, y = Probabilidade, fill = Grupo)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(. ~ Pop_Origem_Facet, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cores_pop, name = "Atribuído a:") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(
    title = paste("Atribuição Genética de Origem (DAPC) -", especie_alvo),
    x = "População de Coleta (% de Atribuição Correta) / Indivíduos",
    y = "Probabilidade de Pertencimento"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.15, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

print(dapc_plot)

file_dapc_tiff <- file.path(pasta_especie, paste0("DAPC_Assignment_", especie_alvo, ".tiff"))
ggsave(file_dapc_tiff, plot = dapc_plot, width = 320, height = 180, units = "mm", dpi = 300)

file_post_csv <- file.path(pasta_especie, paste0("Probabilidades_Assignment_", especie_alvo, ".csv"))
write.csv(dapc_res$posterior, file = file_post_csv, row.names = TRUE)

# ==============================================================================
# 6. Salvamento do Relatório Geral em Texto (Seguro contra travamentos do sink)
# ==============================================================================
file_txt_report <- file.path(pasta_especie, paste0("Relatorio_Resultados_", especie_alvo, ".txt"))

sink(file_txt_report)

cat("==============================================================================\n")
cat("               RELATÓRIO DE RESULTADOS GENÉTICOS -", especie_alvo, "\n")
cat("==============================================================================\n\n")

cat("Arquivo VCF Analisado:", caminho_vcf, "\n")
cat("Número de Indivíduos:", nInd(genind_obj), "\n")
cat("Número de Loci/SNPs:", nLoc(genind_obj), "\n")
cat("Número de Populações Identificadas:", num_pops, "\n\n")

cat("------------------------------------------------------------------------------\n")
cat("1. DIVERSIDADE GENÉTICA POR POPULAÇÃO (N, Ho, He, Fis)\n")
cat("------------------------------------------------------------------------------\n")
print(round(tabela_diversidade, 4))
cat("\n")

cat("------------------------------------------------------------------------------\n")
cat("2. ESTRUTURAÇÃO GENÉTICA (Fst de Weir & Cockerham 1984)\n")
cat("------------------------------------------------------------------------------\n")
cat("Fst Global da Espécie:", round(fst_global, 4), "\n\n")
cat("Matriz de Fst Par a Par entre Populações:\n")
print(fst_pairwise)
cat("\n")

cat("------------------------------------------------------------------------------\n")
cat("3. ANÁLISE DE COMPONENTES PRINCIPAIS (PCA FREQUENCISTA)\n")
cat("------------------------------------------------------------------------------\n")
cat("Variância explicada pelo PC1:", var_eixo1, "%\n")
cat("Variância explicada pelo PC2:", var_eixo2, "%\n\n")

cat("------------------------------------------------------------------------------\n")
cat("4. ANÁLISE BAYESIANA (LEA / sNMF, ANCESTRIA E PCoA)\n")
cat("------------------------------------------------------------------------------\n")
cat("Mínima Entropia Cruzada por K (sNMF):\n")
for (k in 1:length(ce_values)) {
  cat(paste0("  K = ", k, ": ", round(ce_values[k], 5), "\n"))
}
cat("\nK ótimo selecionado (menor entropia): K =", k_otimo, "\n")
cat("Variância explicada pelo Eixo Bayesiano 1:", var_bayes1, "%\n")
cat("Variância explicada pelo Eixo Bayesiano 2:", var_bayes2, "%\n\n")

cat("------------------------------------------------------------------------------\n")
cat("5. DAPC - ATRIBUIÇÃO DE ORIGEM E VALIDAÇÃO CRUZADA\n")
cat("------------------------------------------------------------------------------\n")
cat("Melhor número de PCs retidos (xvalDapc):", best_pc, "\n\n")
cat("Taxa de re-atribuição correta por população (%):\n")
print(assign_pct)
cat("\n")

cat("==============================================================================\n")
cat("               Sessão do R e Versões dos Pacotes -", especie_alvo, "\n")
cat("==============================================================================\n\n")
print(sessionInfo())

sink()

cat("\nAnálise para a espécie", especie_alvo, "concluída com sucesso!")
cat("\nTodos os arquivos e gráficos foram salvos na pasta:", pasta_especie, "\n")