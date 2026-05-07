
#1.Descarga de genoma 
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

#2.Copie las secuencias...
cp difexpp_human /home/bio.pt/data/marianarozor/clases/transcriptomica/expreson_genoma_ref_anotation
#NO HICE FILTROSSSSS.... SERÁ QUE TOCABA O YA ESTABAN FILTRADOS??



#3. Codigo para todo el index_mapeo:
nano index_mapeo.sh
#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -o salida_mapeo_humano_22.out
#SBATCH -e error_mapeo_humano_22.err
#SBATCH --mail-user=mariana.rozor@urosario.edu.co
#SBATCH --mail-type=ALL

#Ruta de donde estan los scribs de HISAT2
APP=/datacnmat01/ciencias/appsbio/conda/envs/appsb/bin

# Creacion de archivos de sites
$APP/hisat2_extract_splice_sites.py chr22_with_ERCC92.gtf > splicesites.tsv

# Creacion de archivos de exones
$APP/hisat2_extract_exons.py chr22_with_ERCC92.gtf > exons.tsv

# Creación del index
$APP/hisat2-build -p 8 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_genome_tran

# Nombres de muestras
ls *fastq.gz | sed 's/\.read.*//' > nombres_archivos.txt

# Alineamientos
for line in $(cat nombres_archivos.txt)
do
$APP/hisat2 -p 8 -x chr22_genome_tran \
-1 $line.read1.fastq.gz \
-2 $line.read2.fastq.gz \
-S $line.sam \
2>> resumen_aln.txt
done

# SAM a BAM
for line in $(cat nombres_archivos.txt)
do
/datacnmat01/ciencias/appsbio/conda/envs/appsb/bin/samtools sort -@ 8 -o $line.bam $line.sam
rm $line.sam
done

rm -f tmp_*


#Para Analizar si todo salio bien:
ls
salloc
module load samtools
for i in *.bam
do
  echo "===== $i ====="
  samtools flagstat $i
done
less less resumen_aln.txt
#Anlaiszar loa nyerior sobre todo loq eusalio despeus del for o ademas ver si es lo msivo que esta en el archivo resumen_aln.txt

#4. Análisis de expresión (matrices y R studio)
nano conteos.sh
#!/bin/bash
#SBATCH -p dev
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -o salida_conteo_matrix.out
#SBATCH -e error_conteo_matrix.err
#SBATCH --mail-user=mariana.rozor@urosario.edu.co
#SBATCH --mail-type=ALL

APP=/datacnmat01/ciencias/appsbio/conda/envs/appsb/bin

for line in $(cat nombres_archivos.txt)
do
$APP/htseq-count \
  --format=bam \
  --order=pos \
  --mode=intersection-strict \
  --stranded=reverse \
  --minaqual=1 \
  --type=exon \
  --idattr=gene_id \
  $line.bam \
  chr22_with_ERCC92.gtf > ${line}_gene.tsv
done

join UHR_Rep1_gene.tsv UHR_Rep2_gene.tsv | join - UHR_Rep3_gene.tsv | join - HBR_Rep1_gene.tsv | join - HBR_Rep2_gene.tsv | join - HBR_Rep3_gene.tsv > gene_read_counts_table_all.tsv
echo "GeneID UHR_Rep1 UHR_Rep2 UHR_Rep3 HBR_Rep1 HBR_Rep2 HBR_Rep3" > header.txt
cat header.txt gene_read_counts_table_all.tsv > head_gene_read_counts_table_all.tsv
head -n -5 head_gene_read_counts_table_all.tsv > final_gene_read_counts_table_all.tsv



# 8. CONCLUSIÓN SIMPLE
Error: falta índice .bai
Impacto: bajo (no rompe el análisis)
Solución: correr samtools index
Resultado: tu pipeline está bien
#ENTENDERRR ESTE ERROR PUES UGUAL FUNCIONO PERO EN  MI PROYECTO NO DEBE PASAR ESTO



module load R 3.6.3
R
library(edgeR)

library(ggrepel)

library(gplots)

library(RColorBrewer)

library(dplyr)

library(viridis)

x <- read.table("final_gene_read_counts_table_all.tsv", header=TRUE)

row.names(x) <- x$GeneID

x1 <- subset(x, select=-GeneID)

x1[] <- lapply(x1, as.numeric)

x <- as.matrix(x1)

group <- factor(c("UHR","UHR","UHR","HBR","HBR","HBR"))

y <- DGEList(counts=x, group=group)

levels(y$samples$group)

design <- model.matrix(~0+group, data=y$samples)

colnames(design) <- levels(y$samples$group)

y <- calcNormFactors(y)

y <- estimateDisp(y, design)

fit <- glmFit(y, design)

my.contrasts <- makeContrasts(HBR_UHR=HBR-UHR, levels=design)

qlf <- glmLRT(fit, contrast=my.contrasts[,"HBR_UHR"])

summary(decideTests(qlf, p.value=0.05))

topTags(qlf)

points <- 23

colors <- c("red","blue")

pdf("plotMDS_2.pdf")

plotMDS(y, col=colors[group], pch=15)

legend("top", legend=levels(group), pch=15, col=colors)

dev.off()

resLRTfilt <- topTags(qlf, n=nrow(y$counts))

volcanoData <- cbind(resLRTfilt$table$logFC, -log10(resLRTfilt$table$FDR))

colnames(volcanoData) <- c("logFC","negLogPval")

DEGs <- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC) > 1

point.col <- ifelse(DEGs, "red", "black")

pdf("volcano.pdf")

plot(volcanoData, pch=16, col=point.col, cex=0.5)

dev.off()

logcounts <- cpm(y, log=TRUE)

var_genes <- apply(logcounts, 1, var)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]

highly_variable_lcpm <- logcounts[select_var,]

dim(highly_variable_lcpm)

mypalette <- scale_fill_viridis(discrete=TRUE)

pdf("heatmap.pdf")

heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples", key=TRUE)

dev.off()


scp -i bio.pt.pem -P 53841 bio.pt@loginpub-hpc.urosario.edu.co:/home/bio.pt/data/marianarozor/clases/transcriptomica/expreson_genoma_ref_anotation/difexpp_human/*.pdf /mnt/c/Users/Mariana\ Rozo/Downloads/


debo analisasr :
resumen_aln.txt o los resultados del for después del alineamiento usando los bam
y los 3 graficoss.odf que descarge
