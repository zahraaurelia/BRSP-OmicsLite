#Modul: Analisis Ekspresi Gen Osteoartritis
#Dataset: GSE55457 (OA vs Normal)
#Platform: Microarray (Affymetrix Human Genome U133A Array)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 

#PART A. PENGANTAR KONSEP 

#Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen 
#antara dua kondisi biologis (misalnya OA vs normal) 
#Pada modul ini kita menggunakan pendekatan statistik limma (Linear Models
#for Microarray Data), yang merupakan standar emas untuk data microarray. 

#PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE) 

#Apa itu package? 
#Package adalah kumpulan fungsi siap pakai di R
#Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor 

#1. Install BiocManager (manajer paket Bioconductor) 
#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi”

if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

# 2. Install paket Bioconductor (GEOquery & limma) 
#GEOquery: mengambil data dari database GEO 
#limma: analisis statistik ekspresi gen 

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

#Install annotation package sesuai platform
#GPL96 = Affymetrix Human Genome U133 Plus 2.0
BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)

#3. Install paket CRAN untuk visualisasi dan manipulasi data 
#phetmap: heatmap ekspresi gen 
#ggplot2: grafik (volcano plot)
#dplyr: manipulasi tabel data 

install.packages(c("pheatmap", "ggplot2", "dplyr"))

#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

#Install preprocessCore
BiocManager::install("preprocessCore")
library(preprocessCore)

#4. Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan 
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
library(umap)
library(preprocessCore)

#PART C. PENGAMBILAN DATA DARI GEO 


#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE55457", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

#PART D. PRE-PROCESSING DATA EKSPRESI 

#pData(): metadata sampel
#source_name_ch1 berisi informasi kondisi biologis sampel
#perlu pre-processing karena hanya dibutuhkan data pasien kontrol normal dan OA
group_info <- pData(gset)[["source_name_ch1"]]
groups <- sub(".*\\((.*)\\).*", "\\1", group_info)
unique(groups)
keep_samples <- grepl("osteoarthrisits", groups, ignore.case = TRUE) |
  grepl("normal controls", groups, ignore.case = TRUE)
gset <- gset[, keep_samples]
groups <- groups[keep_samples]

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
ex <- exprs(gset)
exprs(gset) <- ex

#Mengapa perlu log2 transformasi?
#Data microarray mentah memiliki rentang nilai sangat besar.
#Log2 digunakan untuk:
#1. Menstabilkan varians
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change

#quantile(): menghitung nilai kuantil (persentil)
#as.numeric(): mengubah hasil quantile (yang berupa named vector)
#menjadi vektor numerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

#LogTransform adalah variabel logika (TRUE / FALSE)
#Operator logika:
#>  : lebih besar dari
#|| : OR (atau)
#&& : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#normalisasi quantile agar distribusi seragam
ex <- normalizeBetweenArrays(ex, method = "quantile")


#PART E. DEFINISI KELOMPOK SAMPEL 

#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(groups)

#factor():
#Mengubah data kategorik menjadi faktor
#Faktor sangat penting untuk analisis statistik di R
gset$group <- factor(groups)

#levels(): melihat kategori unik dalam faktor
nama_grup <- levels(gset$group)
print(nama_grup)

#PART F. DESIGN MATRIX (KERANGKA STATISTIK) 

#model.matrix():
#Membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + gset$group)

#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)

#Menentukan perbandingan biologis
grup_normal <- nama_grup[1]
grup_OA <- nama_grup[2]



contrast_formula <- paste(grup_OA, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)

#lmFit():
#Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)

#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)

#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01  -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)



#PART H. ANOTASI NAMA GEN 

#Penting:
#Pada data microarray Affymetrix, unit analisis awal adalah PROBE,
#bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan
#database resmi Bioconductor.

#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI 

#Boxplot digunakan untuk:
#- Mengecek distribusi nilai ekspresi antar sampel
#- Melihat apakah ada batch effect
#- Mengevaluasi apakah normalisasi/log-transform sudah wajar

# Reset bersih dari gset
ex <- exprs(gset)

# Cek dulu, jangan log2 lagi kalau sudah log2
print(range(ex, na.rm = TRUE))
# Kalau range 2-16 = JANGAN log2 lagi
# Kalau range > 100  = baru log2
# Data butuh di log2 kembali
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#Penentuan warna
group_colors <- c("cornflowerblue", "darkorchid2")

#Mapping warna ke tiap sampel
group_numeric <- as.numeric(gset$group)
sample_colors <- group_colors[group_numeric]

boxplot(ex, col = sample_colors, las = 2, outline = FALSE,
        main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
        ylab = "Expression Value (log2)", cex.axis = 0.5)
legend("topright", legend = levels(gset$group),
       fill = group_colors, cex = 0.8)


#PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 

#Density plot menunjukkan sebaran global nilai ekspresi gen
#Digunakan untuk:
#- Mengecek efek log-transform
#- Membandingkan distribusi antar grup

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

print(qx)
print(LogTransform)
print(range(ex, na.rm = TRUE))

# Simpan filtered ke variabel baru, jangan timpa ex
ex_filtered <- ex[rowMeans(ex, na.rm = TRUE) > 3, ]

# Gunakan ex_filtered untuk density plot
expr_long <- data.frame(
  Expression = as.vector(ex_filtered),
  Group = rep(gset$group, each = nrow(ex_filtered))
)

# Buat density plot
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (Filtered)",
    x = "Expression Value (log2)",
    y = "Density"
  )

#PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)

#UMAP digunakan untuk:
#- Mereduksi ribuan gen menjadi 2 dimensi
#- Melihat pemisahan sampel secara global
#- Alternatif PCA (lebih sensitif ke struktur lokal)

#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#Jalankan UMAP
umap_result <- umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )


#PART J.1 VISUALISASI VOLCANO PLOT 

#Volcano plot menggabungkan:
#- Log fold change (efek biologis)
#- Signifikansi statistik

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG OA")


#PART J.2 VISUALISASI HEATMAP 

#Heatmap digunakan untuk melihat pola ekspresi gen
#antar sampel berdasarkan gen-gen paling signifikan

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

# Reset ex dari gset asli
ex <- exprs(gset)

# Log2 transform ulang jika perlu
ex[ex <= 0] <- NA
ex <- log2(ex)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)


#PART K. MENYIMPAN HASIL 

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE55457_DEG.csv")



#PART L. Analisis KEGG dan Enrichment

#Install ClusterProfiler
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "AnnotationDbi"),
                     ask = FALSE, update = FALSE)

install.packages(c("ggrepel"))

library(AnnotationDbi)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

sig_deg <- topTableResults %>%
  filter(!is.na(SYMBOL), SYMBOL != "",
         adj.P.Val < 0.05, abs(logFC) > 1)

cat("\nGen untuk enrichment:", nrow(sig_deg), "\n")

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys      = sig_deg$SYMBOL,
                     column    = "ENTREZID",
                     keytype   = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)
cat("Entrez IDs:", length(entrez_ids), "\n")

# GO Biological Process
go_bp <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

# GO BP - tampil langsung di RStudio
barplot(go_bp, showCategory = 20,
        title = "GO Biological Process – Osteoartritis (GSE55457)")

dotplot(go_bp, showCategory = 20,
        title = "GO Biological Process – Dot Plot")

# GO Molecular Function
go_mf <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

# 2. LANGKAH PEMASTIAN (Cek apakah ada hasil yang signifikan)
# Gunakan fungsi nrow() untuk melihat jumlah baris hasil
jumlah_mf <- nrow(as.data.frame(go_mf))
print(paste("Jumlah kategori MF yang signifikan:", jumlah_mf))

# 3. Visualisasi Bersyarat
if (jumlah_mf > 0) {
  # Jika ada data, buat plot
  barplot(go_mf, showCategory = 20,
          title = "GO Molecular Function – Osteoartritis (GSE55457)")
} else {
  # Jika kosong, berikan pesan agar Anda tidak bingung mencari plot yang tidak ada
  message("Peringatan: Tidak ada istilah GO MF yang signifikan pada p < 0.05.")

# KEGG Pathway
kegg_res <- enrichKEGG(gene          = entrez_ids,
                       organism      = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05)
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
barplot(kegg_res, showCategory = 20,
        title = "KEGG Pathway – Osteoartritis (GSE55457)")

dotplot(kegg_res, showCategory = 20,
        title = "KEGG Pathway – Dot Plot")

cnetplot(kegg_res, showCategory = 10)

# Simpan tabel hasil enrichment
write.csv(as.data.frame(go_bp),    "GO_BP_results.csv",  row.names = FALSE)
write.csv(as.data.frame(go_mf),    "GO_MF_results.csv",  row.names = FALSE)
write.csv(as.data.frame(kegg_res), "KEGG_results.csv",   row.names = FALSE)


message("✅ Analisis selesai! Semua file telah disimpan.")


