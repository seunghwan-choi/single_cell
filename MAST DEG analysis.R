rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(Matrix)
library(MAST)
library(SingleCellExperiment)

# Date : 24.11.26

# Load Data (Scanpy Python Data)
counts = readMM('./Downloads/nac_scrna/data_normalizedX_neuronal_normalized_Lv0Celltype.mtx')
cell_metadata <- read.csv("./Downloads/nac_scrna/data_obs_neuronal_normalized_Lv0Celltype.csv", row.names = 1)
gene_metadata <- read.csv("./Downloads/nac_scrna/data_var_neuronal_normalized_Lv0Celltype.csv", row.names = 1)

# Counts colnames and rownames setting
rownames(counts) = rownames(cell_metadata)
colnames(counts) = rownames(gene_metadata)
counts_t=t(counts)

# Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = counts_t)
seurat_obj <- AddMetaData(seurat_obj, metadata = cell_metadata)
seurat_obj[['RNA']]@misc$gene_metadata = gene_metadata

# Convert SeuratObject to SingleCellExperiment object
sce_obj = as.SingleCellExperiment(seurat_obj)

# metadata check
colData(sce_obj)$Condition%>%head()
colData(sce_obj)$siteid%>%head()
colData(sce_obj)$sample%>%head()

# Remove vHipp_BLA_Both and VTA_PFC_Both Cells => sce_obj_filtered 생성
cells_to_keep <- which(!(seurat_obj@meta.data$siteid %in% c("vHipp_BLA_Both", "VTA_PFC_Both")))
seurat_obj_filtered <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[cells_to_keep])
sce_obj_filtered = as.SingleCellExperiment(seurat_obj_filtered)

# Convert SingleCellExperiment data into SingleCellAssay object
sca_obj_filtered = SceToSingleCellAssay(sce_obj_filtered)

# Check siteid data
colData(sca_obj_filtered)$siteid

# DEG analysis (Site vs Others)
# 비교할 siteid 목록
siteid_list <- c("vHipp_BLA_EYFP", "vHipp_BLA_mCherry", "VTA_PFC_EYFP", "VTA_PFC_mCherry")

# 결과를 저장할 리스트 초기화
deg_results_list <- list()

# 각 siteid에 대해 DEG 분석 수행
for (site in siteid_list) {
  
  # siteid_group 변수 생성: 현재 site와 나머지 비교
  colData(sca_obj_filtered)$siteid_group <- ifelse(colData(sca_obj_filtered)$siteid == site, site, "Other")
  
  # zlm 모델 생성: siteid_group을 주요 변수로 설정하고 sample을 covariate로 추가
  zlmCond <- zlm(~ siteid_group + sample, sca_obj_filtered)
  
  # DEG 결과 추출: siteid_group에서 현재 site와 나머지 그룹의 차이 검정
  summaryCond <- summary(zlmCond, doLRT = paste0("siteid_group", site))
  summaryCondTable <- summaryCond$datatable
  
  # 유의미한 DEG 필터링
  deg_results <- summaryCondTable[contrast == paste0("siteid_group", site) & `Pr(>Chisq)` < 0.05, ]
  
  # 결과를 리스트에 저장
  deg_results_list[[site]] <- deg_results
}


# 병렬 처리 패키지 로드
library(BiocParallel)

# 병렬 처리 환경 설정
BPPARAM <- MulticoreParam(workers = 7)  # 사용하려는 코어 수 지정 (예: 4)

# 각 siteid에 대해 DEG 분석 수행
for (site in siteid_list) {
  
  # siteid_group 변수 생성: 현재 site와 나머지 비교
  colData(sca_obj_filtered)$siteid_group <- ifelse(colData(sca_obj_filtered)$siteid == site, site, "Other")
  
  # zlm 모델 생성: siteid_group을 주요 변수로 설정하고 sample을 covariate로 추가
  zlmCond <- zlm(~ siteid_group + sample, sca_obj_filtered, method = "bayesglm", BPPARAM = BPPARAM)
  
  # DEG 결과 추출: siteid_group에서 현재 site와 나머지 그룹의 차이 검정
  summaryCond <- summary(zlmCond, doLRT = paste0("siteid_group", site))
  summaryCondTable <- summaryCond$datatable
  
  # FDR 보정 (전체 p-value를 기준으로 FDR 계산)
  summaryCondTable$FDR <- p.adjust(summaryCondTable$`Pr(>Chisq)`, method = "BH")
  
  # 유의미한 DEG 필터링 (FDR < 0.05)
  deg_results <- summaryCondTable[contrast == paste0("siteid_group", site) & FDR < 0.05, ]
  
  # 결과를 리스트에 저장
  deg_results_list[[site]] <- deg_results
}

# 각 사이트별 DEG 결과 확인
deg_results_list$vHipp_BLA_EYFP # Site1
deg_results_list$vHipp_BLA_mCherry # Site2
deg_results_list$VTA_PFC_EYFP # Site3
deg_results_list$VTA_PFC_mCherry # Site4

## save DEG results (!!)
saveRDS(deg_results_list, file = './Downloads/nac_scrna/deg_results_list_fdr_v2.Rdata', version=2) # version2 setting: able to read also in R4.3.3 

k= readRDS('./Downloads/nac_scrna/deg_results_list_fdr_v2.Rdata')

# DEG results for each site
site1 = deg_results_list$vHipp_BLA_EYFP # Site1 
site2 = deg_results_list$vHipp_BLA_mCherry # Site2
site3 = deg_results_list$VTA_PFC_EYFP # Site3
site4 = deg_results_list$VTA_PFC_mCherry # Site4

deg_reults_list
 