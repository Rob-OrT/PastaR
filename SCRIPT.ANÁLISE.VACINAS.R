# TCC – ROBSON LUIZ DE OLIVEIRA
# ANÁLISE DE COBERTURA VACINAL INFANTIL NO BRASIL (2020–2024)
# BLOCO 1: ANÁLISE REGIONAL TEMPORAL
# BLOCO 2: CLUSTERIZAÇÃO ESTADUAL (2024)
# Linguagem: R (versão 4.3.3)
# =====================================================================

# Carrego todos os pacotes que utilizei ao longo da minha análise:
library(tidyverse)
library(reshape2)
library(cluster)
library(factoextra)
library(ggplot2)
library(dendextend)

# ====================================================
# BLOCO 1 – ANÁLISE REGIONAL TEMPORAL (2020–2024)
# ====================================================

#REGIAO
dados_regionais <- tribble(
  ~Regiao, ~Ano, ~BCG, ~HBV, ~DTP_HB_Hib, ~IPV, ~PCV10, ~HRV, ~MenC,
  "Centro-Oeste", 2020, 80.50, 71.50, 80.42, 80.47, 86.57, 81.72, 83.71,
  "Centro-Oeste", 2021, 78.80, 72.42, 74.50, 74.22, 79.17, 75.74, 76.15,
  "Centro-Oeste", 2022, 90.44, 85.22, 80.69, 80.50, 87.59, 81.19, 83.84,
  "Centro-Oeste", 2023, 91.66, 88.06, 85.69, 86.39, 90.50, 87.36, 90.24,
  "Centro-Oeste", 2024, 93.91, 92.01, 89.63, 91.39, 95.92, 91.70, 89.28,
  "Nordeste", 2020, 74.83, 70.34, 70.11, 73.11, 79.72, 74.72, 76.09,
  "Nordeste", 2021, 75.13, 72.12, 69.47, 68.53, 72.67, 69.16, 69.41,
  "Nordeste", 2022, 97.54, 89.34, 78.86, 78.50, 82.35, 75.90, 78.86,
  "Nordeste", 2023, 88.95, 84.17, 89.41, 90.26, 92.08, 88.74, 90.45,
  "Nordeste", 2024, 94.36, 93.20, 89.45, 89.12, 91.73, 88.13, 87.00,
  "Norte", 2020, 81.19, 74.18, 64.19, 65.69, 76.19, 68.54, 71.53,
  "Norte", 2021, 80.66, 75.73, 62.35, 62.29, 69.44, 63.94, 66.06,
  "Norte", 2022, 96.63, 87.52, 71.59, 71.24, 78.93, 70.09, 74.64,
  "Norte", 2023, 87.59, 83.55, 79.13, 81.40, 86.43, 80.10, 83.75,
  "Norte", 2024, 93.34, 89.92, 83.84, 83.44, 83.35, 81.33, 81.67,
  "Sudeste", 2020, 73.41, 57.91, 83.37, 78.28, 81.29, 78.66, 79.16,
  "Sudeste", 2021, 71.18, 60.48, 71.86, 71.53, 73.83, 71.76, 71.77,
  "Sudeste", 2022, 83.31, 76.20, 74.79, 75.14, 77.85, 75.27, 75.98,
  "Sudeste", 2023, 83.32, 77.47, 86.82, 87.61, 88.78, 87.28, 89.33,
  "Sudeste", 2024, 91.90, 89.04, 90.52, 90.47, 91.94, 89.44, 89.21,
  "Sul", 2020, 87.45, 69.02, 87.99, 86.50, 90.79, 87.52, 89.16,
  "Sul", 2021, 78.42, 64.66, 80.89, 79.98, 83.94, 81.34, 81.54,
  "Sul", 2022, 88.33, 82.31, 83.31, 83.10, 88.42, 84.23, 85.61,
  "Sul", 2023, 91.09, 86.60, 92.06, 92.31, 94.45, 92.75, 91.75,
  "Sul", 2024, 94.94, 94.30, 92.74, 92.81, 94.38, 91.66, 89.69
)

# Transformar em formato longo para visualização
dados_longos <- melt(dados_regionais, id.vars = c("Regiao", "Ano"),
                     variable.name = "Vacina", value.name = "Cobertura")

# Gráfico de linhas: evolução temporal por vacina
ggplot(dados_longos, aes(x = Ano, y = Cobertura, color = Regiao)) +
  geom_line(size = 1.2) +
  facet_wrap(~ Vacina) +
  theme_minimal() +
  labs(title = "Evolução da Cobertura Vacinal por Região e Imunizante (2020–2024)",
       x = "Ano", y = "Cobertura (%)", color = "Região")

# Padronização por Z-score para análise de clustering
dados_z <- dados_regionais %>%
  select(-Regiao, -Ano) %>%
  scale()

# Criar objeto com informações de região e ano
info <- dados_regionais %>%
  select(Regiao, Ano)

# Clustering hierárquico com método de Ward
dist_matrix <- dist(dados_z)
hc <- hclust(dist_matrix, method = "ward.D2")

# Dendrograma com rótulos
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 3)
plot(dend, main = "Dendrograma Hierárquico – Séries Regionais")

# Extração dos clusters (k = 3)
grupos <- cutree(hc, k = 3)
info$Cluster <- factor(grupos)

# Média por cluster (padronizada)
medias_cluster1 <- aggregate(dados_z, by = list(info$Cluster), FUN = mean)
print(medias_cluster1)

# ANOVA por vacina
for (vac in colnames(dados_z)) {
  print(anova(aov(dados_z[, vac] ~ info$Cluster)))
}


# BLOCO 2 – CLUSTERIZAÇÃO DOS ESTADOS EM 2024
# ====================================================

# estados
dados_estaduais <- tribble(
  ~Estado, ~BCG, ~HBV, ~DTP_HB_Hib, ~IPV, ~PCV10, ~HRV, ~MenC,
  "AC", 95.22, 93.21, 87.04, 85.78, 93.06, 83.24, 87.36,
  "AL", 104.80, 105.42, 90.98, 91.05, 93.35, 90.68, 90.59,
  "AM", 108.83, 105.64, 87.55, 86.41, 91.53, 81.94, 82.76,
  "AP", 105.94, 105.26, 71.41, 72.09, 80.39, 71.42, 73.35,
  "BA", 81.22, 79.73, 88.55, 88.05, 90.93, 87.07, 82.11,
  "CE", 113.55, 112.97, 94.92, 94.26, 95.80, 93.05, 92.37,
  "DF", 119.93, 90.12, 92.00, 90.71, 97.22, 92.68, 93.90,
  "ES", 89.03, 87.38, 95.49, 94.96, 98.53, 94.96, 77.30,
  "GO", 77.30, 73.93, 83.09, 83.33, 88.68, 84.51, 82.36,
  "MA", 76.59, 73.66, 85.45, 84.89, 83.47, 82.80, 83.47,
  "MG", 97.38, 95.30, 94.53, 94.11, 95.04, 93.46, 92.76,
  "MS", 104.10, 102.32, 99.11, 98.37, 101.43, 97.17, 100.00,
  "MT", 98.22, 95.56, 91.85, 90.30, 94.45, 90.46, 90.79,
  "PA", 83.93, 79.01, 80.22, 80.14, 85.01, 79.32, 78.56,
  "PB", 100.06, 99.15, 86.75, 86.73, 89.15, 86.46, 86.10,
  "PE", 102.13, 101.52, 88.72, 89.01, 92.34, 88.50, 87.96,
  "PI", 97.02, 96.98, 95.01, 94.12, 95.75, 92.11, 93.82,
  "PR", 99.97, 103.23, 94.65, 94.17, 95.33, 93.04, 90.27,
  "RJ", 90.52, 85.27, 80.37, 80.12, 84.44, 80.87, 80.12,
  "RN", 98.65, 96.75, 89.08, 88.91, 93.03, 90.04, 90.73,
  "RO", 92.06, 90.46, 94.86, 94.39, 97.01, 93.03, 95.37,
  "RR", 78.60, 75.58, 72.67, 72.91, 81.27, 68.71, 63.63,
  "RS", 100.62, 97.08, 92.52, 90.20, 95.01, 91.71, 90.20,
  "SC", 81.49, 79.31, 90.04, 91.56, 92.20, 89.37, 88.07,
  "SE", 109.72, 108.20, 90.83, 90.62, 92.05, 89.38, 90.21,
  "SP", 90.73, 88.30, 92.00, 91.97, 92.73, 90.21, 90.37,
  "TO", 106.61, 106.02, 92.00, 91.55, 93.73, 89.24, 92.12
)

# Padronização dos dados estaduais
vacinas_uf <- scale(dados_estaduais[, -1])
rownames(vacinas_uf) <- dados_estaduais$Estado

# Método do cotovelo
fviz_nbclust(vacinas_uf, kmeans, method = "wss") +
  labs(title = "Método do Cotovelo – Estados (2024)")

# K-means com k = 3
set.seed(42)
km <- kmeans(vacinas_uf, centers = 3, nstart = 25)

# Anexar cluster aos dados
dados_estaduais$Cluster <- factor(km$cluster)

# Boxplot padronizado por cluster
dados_longos_uf <- cbind(Estado = dados_estaduais$Estado,
                         Cluster = dados_estaduais$Cluster,
                         as.data.frame(vacinas_uf)) %>%
  pivot_longer(cols = BCG:MenC, names_to = "Vacina", values_to = "Z")

ggplot(dados_longos_uf, aes(x = Vacina, y = Z, fill = Cluster)) +
  geom_boxplot() +
  labs(title = "Boxplot das Coberturas Vacinais Padronizadas por Cluster (2024)",
       y = "Z-score", x = "Vacina") +
  theme_minimal()

# ANOVA por imunizante
for (vac in colnames(vacinas_uf)) {
  print(anova(aov(vacinas_uf[, vac] ~ dados_estaduais$Cluster)))
}


# ---- Análise Descritiva Regional (2020–2024) ----
cat("===== ANÁLISE DESCRITIVA REGIONAL =====\n")
resumo_regional <- dados_longos %>%
  group_by(Regiao, Vacina) %>%
  summarise(
    Média = round(mean(Cobertura), 2),
    Mediana = round(median(Cobertura), 2),
    Mínimo = min(Cobertura),
    Máximo = max(Cobertura),
    Desvio_Padrão = round(sd(Cobertura), 2),
    Amplitude = round(max(Cobertura) - min(Cobertura), 2),
    .groups = "drop"
  )



