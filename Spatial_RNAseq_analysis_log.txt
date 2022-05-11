Slided ID: V10N16-076
###### 1. Generating FASTQs with spaceranger mkfastq (done by Alberto)
##===================================================
spaceranger mkfastq --id=peanut_mkfastq \
                    --filter-dual-index\
                    --run=rawData \
                    --csv=spaceranger-tiny-bcl-samplesheet-1.0.0.csv\
                    --localcores=30\
                    --localmem=200\
                    --force-single-index
# normal demultiplex failed, as many reads are undetermined, and sample C1 is lost. It is because of the mismatch on i5 index
# use --force-single-index to solve this issue, which only uses i7 index, and delete i5 index in the sample file
# result shows that not many extra reads are identified, and C1 is still missing

spaceranger mkfastq --id=peanut_mkfastq3 \
                    --filter-dual-index\
                    --run=rawData \
                    --csv=spaceranger-tiny-bcl-samplesheet-1.0.0.csv\
                    --localcores=40\
                    --localmem=200\
                    # --force-single-index \
                    --barcode-mismatches=2



# used bcl2fastq to demultiplex, using i7+i5 index
bcl2fastq  \
  --minimum-trimmed-read-length=7 \
  --mask-short-adapter-reads=7 \
  --create-fastq-for-index-reads \
  -r 40 -w 40 \
  -R rawData \
  --output-dir=peanut_bcl2fastq_3_2lanes \
  --barcode-mismatches=2\
  --adapter-stringency=0.7\
  --sample-sheet=./spaceranger-tiny-bcl-samplesheet-1.0.0.csv \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  # --ignore-missing-bcls \
  # --interop-dir=${INTEROP_DIR} \


bcl2fastq  \
  --minimum-trimmed-read-length=7 \
  --mask-short-adapter-reads=7 \
  --create-fastq-for-index-reads \
  -r 40 -w 40 \
  -R rawData \
  --output-dir=all_pooled_libraries \
  --barcode-mismatches=1\
  --adapter-stringency=0.7\
  --sample-sheet=./SampleSheet-Lanes1_2.csv \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \

#using the truncated 8bp barcode
bcl2fastq  \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --create-fastq-for-index-reads \
          -r 20 -w 20 \
          -R bcl_secondDownload \
          --output-dir=peanut_bcl2fastq_8bp \
          --sample-sheet=./spaceranger-samplesheet-8bp.csv \
          --use-bases-mask=Y28n*,I8n*,I8n*,Y* \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \


fastqc -t 30 -o QC *.fastq.gz

fix the unpaired reads
ml bbmap
repair.sh in=D1_S4_L001_R1_001.fastq.gz in2=D1_S4_L001_R2_001.fastq.gz out=D1_S4_L001_R1_001_clean.fastq.gz out2=D1_S4_L001_R2_001_clean.fastq.gz

java -jar validatefastq-assembly-0.1.1.jar -i peanut_bcl2fastq_2/NS2612-JWang-S4-HFMGGDSX3-Lane1-2/D1_S4_L001_R1_001.fastq.gz -j peanut_bcl2fastq_2/NS2612-JWang-S4-HFMGGDSX3-Lane1-2/D1_S4_L001_R2_001.fastq.gz 
###### 2. aignment and count using spaceranger count
##===================================================
##prepare transcriptome

# Filter transcriptome with mkgtf
# spaceranger mkgtf GCF_003086295.2_arahy.Tifrunner.gnm1.KYV3_genomic.gtf GCF_003086295.2_arahy.Tifrunner.gnm1.KYV3_genomic.filtered.gtf --attribute=gene_biotype:protein_coding \
#                    --attribute=gene_biotype:lincRNA \
#                    --attribute=gene_biotype:protein_coding \
#                    --attribute=gene_biotype:lincRNA \
#                    --attribute=gene_biotype:antisense \
#                    --attribute=gene_biotype:IG_LV_gene \
#                    --attribute=gene_biotype:IG_V_gene \
#                    --attribute=gene_biotype:IG_V_pseudogene \
#                    --attribute=gene_biotype:IG_D_gene \
#                    --attribute=gene_biotype:IG_J_gene \
#                    --attribute=gene_biotype:IG_J_pseudogene \
#                    --attribute=gene_biotype:IG_C_gene \
#                    --attribute=gene_biotype:IG_C_pseudogene \
#                    --attribute=gene_biotype:TR_V_gene \
#                    --attribute=gene_biotype:TR_V_pseudogene \
#                    --attribute=gene_biotype:TR_D_gene \
#                    --attribute=gene_biotype:TR_J_gene \
#                    --attribute=gene_biotype:TR_J_pseudogene \
#                    --attribute=gene_biotype:TR_C_gene

# the chromosome name of GCF_003086295 gtf file is different, if use this newly download version, need to change all chr names.


# Index with spaceranger mkref
grep -v "Pltd\|spliceosomal RNA\|small nucleolar RNA\|tRNAscan" NCBI_arahy.Tifrunner.gnm1.KYV3.gtf > NCBI_chr_gene_only.gtf
spaceranger mkref --nthreads=30 --memgb=100 --genome=NCBI_spaceranger --fasta=arahy.Tifrunner.gnm1.KYV3.genome_main.fa --genes=NCBI_chr_gene_only.gtf

#count
spaceranger count --id=ahy_spt_a1 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/Fastq/outs/fastq_path  \
               --sample=A1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/A1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=A1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-A1.json \
               --localcores=20 \
               --localmem=150
spaceranger count --id=ahy_spt_a1 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/Fastq/outs/fastq_path  \
               --sample=A1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/A1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=A1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-A1.json \
               --localcores=20 \
               --localmem=150

spaceranger count --id=ahy_spt_d1_all \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/Fastq/outs/fastq_path  \
               --sample=D1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/D1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=D1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-D1_all.json \
               --localcores=30 \
               --localmem=100

spaceranger count --id=ahy_spt_d1 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/Fastq/outs/fastq_path  \
               --sample=D1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/D1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=D1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-D1.json \
               --localcores=30 \
               --localmem=100

#align the bcl2fq extracted data
spaceranger count --id=ahy_spt_d1_bcl2fq1 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/peanut_bcl2fastq_2 \
               --sample=D1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/D1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=D1 \
               --r1-length=28 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-D1.json \
               --localcores=30 \
               --localmem=100

#error:
FASTQ header mismatch detected at line 3946268 of input files "/blue/wang/luoziliang/nodulation/SpatialSeq/peanut_bcl2fastq_2/NS2612-JWang-S4-HFMGGDSX3-Lane1-2/D1_S4_L001_R1_001.fastq.gz" and "/blue/wang/luoziliang/nodulation/SpatialSeq/peanut_bcl2fastq_2/NS2612-JWang-S4-HFMGGDSX3-Lane1-2/D1_S4_L001_I1_001.fastq.gz"
R1: @A00916:214:HFMGGDSX3:1:2443:18484:1031 1:N:0:NTCTAGCGAG+GATGAAGAGT
I1: @A00916:214:HFMGGDSX3:1:2415:12201:1078 1:N:0:CTCTAGCGAG+GATGAAGAGT
# manually change the N to C


#align the mkfastq 2mm extracted data
spaceranger count --id=ahy_spt_d1_mkfastq3 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/peanut_mkfastq3 \
               --sample=D1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/D1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=D1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-D1.json \
               --localcores=40 \
               --localmem=100
               # --r1-length=28 \

spaceranger count --id=ahy_spt_d1_whole_mkfastq3 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/peanut_mkfastq3 \
               --sample=D1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/D1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=D1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-D1_whole.json \
               --localcores=40 \
               --localmem=100
               # --r1-length=28 \
               # --r2-length=90\


spaceranger count --id=ahy_spt_a1_mkfastq3 \
               --transcriptome=/blue/wang/luoziliang/nodulation/genome/NCBI_spaceranger \
               --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/peanut_mkfastq3 \
               --sample=D1 \
               --image=/blue/wang/luoziliang/nodulation/SpatialSeq/A1_stiched_BF.tif \
               --slide=V10N16-076 \
               --area=A1 \
               --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-A1.json \
               --localcores=40 \
               --localmem=100


# ##test tifrunner version2 
  # # make the reference 
  # cd /blue/wang/luoziliang/nodulation/genome/tif_v2/

  # gffread arahy.Tifrunner.gnm2.ann1.4K0L.gene_models_main.gff3 -T -o arahy.Tifrunner.gnm2.ann1.4K0L.gene_models_main.gtf

  # spaceranger mkref --nthreads=30 --memgb=100 --genome=Tifv2_spaceranger --fasta=arahy.Tifrunner.gnm2.J5K5.genome_main.fna --genes=arahy.Tifrunner.gnm2.ann1.4K0L.gene_models_main.gtf

  # # run the counting
  # cd /blue/wang/luoziliang/nodulation/SpatialSeq

  # spaceranger count --id=ahy_spt_a1_test2 \
  #                --transcriptome=/blue/wang/luoziliang/nodulation/genome/tif_v2/Tifv2_spaceranger \
  #                --fastqs=/blue/wang/luoziliang/nodulation/SpatialSeq/Fastq/outs/fastq_path  \
  #                --sample=A1 \
  #                --image=/blue/wang/luoziliang/nodulation/SpatialSeq/A1_stiched_BF.tif \
  #                --slide=V10N16-076 \
  #                --area=A1 \
  #                --loupe-alignment=/blue/wang/luoziliang/nodulation/SpatialSeq/V10N16-076-A1.json \
  #                --localcores=30 \
  #                --localmem=100

  # # The alignment rate to peanut genome is low ~50%, will directly mapp the reads to peanut genome using other aligner
  # STAR --runThreadN 30 --readFilesCommand zcat --genomeDir /blue/wang/luoziliang/nodulation/genome/star_index/ --outFileNamePrefix A1 --readFilesIn ./A1_S1_L001_R2_001.fastq.gz --sjdbGTFfile /blue/wang/luoziliang/nodulation/genome/arahy.Tifrunner.gnm1.ann1.CCJH.gene_models_main.gff3 --twopassMode Basic --outSAMstrandField intronMotif --outSAMtype BAM Unsorted
  # # result shows that the alignment rate is still low at ~50%
  # # will trim the reads before alignment.
  # cutadapt -a file:/blue/wang/luoziliang/nodulation/genome/adaptor.1.19.fasta -j 30 -m 50 -e 0.2 -o trimmed_reads/A1_S1_L001_R2_clean.fastq.gz A1_S1_L001_R2_001.fastq.gz 
  # # === Summary ===

  # # Total reads processed:                 865,763
  # # Reads with adapters:                   715,447 (82.6%)
  # # Reads that were too short:                 543 (0.1%)
  # # Reads written (passing filters):       865,220 (99.9%)

  # # Total basepairs processed:    78,784,433 bp
  # # Total written (filtered):     75,004,453 bp (95.2%)

  # #align trimmed reads
  # STAR --runThreadN 30 --readFilesCommand zcat --genomeDir /blue/wang/luoziliang/nodulation/genome/star_index/ --outFileNamePrefix A1_trimmed --readFilesIn trimmed_reads/A1_S1_L001_R2_clean.fastq.gz --sjdbGTFfile /blue/wang/luoziliang/nodulation/genome/arahy.Tifrunner.gnm1.ann1.CCJH.gene_models_main.gff3 --twopassMode Basic --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

  # # merge reads from two lanes
  # cat D1_S4_L001_R2_001.fastq.gz D1_S4_L002_R2_001.fastq.gz > D1_merged_R2.fastq.gz
  # # align to genome
  # STAR --runThreadN 40 --readFilesCommand zcat --outFilterScoreMinOverLread 0.6 --genomeDir /blue/wang/luoziliang/nodulation/genome/star_index/ --outFileNamePrefix D1 --readFilesIn D1_merged_R2.fastq.gz --sjdbGTFfile /blue/wang/luoziliang/nodulation/genome/arahy.Tifrunner.gnm1.ann1.CCJH.gene_models_main.gff3 --twopassMode Basic --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

############check undertermined nthreads##################

ml seqtk
seqtk seq -a Undetermined_S0_L001_I1_001.fastq.gz > Undetermined_S0_L001_I1_001.fasta
formatdb -i Undetermined_S0_L001_I1_001.fasta -p F
module load  ncbi_blast
blastall -a 10 -p blastn -d Undetermined_S0_L001_I1_001.fasta -i index.fa -o index2read.blst -e 1e-6 -m 9



######3. Analyze the spanceranger count data using Seurat################
##=======================================================================

ml R/4.1
R
#######R environment
#####load libraries
library(Seurat)
# library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
#####load the data
d1 <- Load10X_Spatial("/blue/wang/luoziliang/nodulation/SpatialSeq/ahy_spt_d1_mkfastq1/outs",
       filename = "filtered_feature_bc_matrix.h5",
       assay = "Spatial",
       slice = "D1",
       filter.matrix = TRUE,
       to.upper = FALSE,
       image = NULL)
# using sctransform (Hafemeister and Satija, Genome Biology 2019), which which builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance. 
d1 <- SCTransform(d1, assay = "Spatial", verbose = FALSE)

######plot the molecular(RNA) counts 
pdf("molecule_count.pdf")
plot1 <- VlnPlot(d1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(d1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.off()

###### Gene expression visualization
pdf('gene_expression.pdf')
SpatialFeaturePlot(d1, features = c("LOC112703831", "LOC112775219"), alpha = c(0.1, 1),pt.size.factor = 1.6)
dev.off()
# features: Name of the feature to visualize. Provide either group.by OR features, not both. (just put gene name here)
# group.by: Name of meta.data column to group the data by
# pt.size.factor- This will scale the size of the spots. Default is 1.6
# alpha - minimum and maximum transparency. Default is c(1, 1).
# Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression

###### Dimensionality reduction, clustering, and visualization

d1 <- RunPCA(d1, assay = "SCT", verbose = FALSE)
d1 <- FindNeighbors(d1, reduction = "pca", dims = 1:30)
d1 <- FindClusters(d1, verbose = FALSE)
d1 <- RunUMAP(d1, reduction = "pca", dims = 1:30)

# visualize the results of the clustering either in UMAP 
pdf('umap.pdf')
p1 <- DimPlot(d1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(d1, label = TRUE, label.size = 3)
# label parameter places a colored box at the median of each cluste
p1 + p2
dev.off()

pdf('umap_cell_types.pdf')
SpatialDimPlot(d1, cells.highlight = CellsByIdentities(object = d1, idents = c(0,1,2)), 
	facet.highlight = TRUE,
	alpha = c(0.5, 1))
dev.off()
#use cells.highlight parameter to demarcate particular cells of interest, idents are the cluster #

## Interactive plotting
SpatialDimPlot(d1, interactive = TRUE)


# Identification of Spatially Variable Features
# To identify molecular features that correlate with spatial location. 
# 1.perform differential expression based on pre-annotated anatomical regions.(clusters)
de_markers <- FindMarkers(d1, ident.1 = 0, ident.2 = 1)
#take top 3 marker to differentiate the tissue types
pdf('top_marker_gene_predefined_expression.pdf')
SpatialFeaturePlot(object = d1, 
				features = rownames(de_markers)[1:6],
				alpha = c(0.1, 1), 
				# ncol = 3,
				)
dev.off()
# 2. alternative approach, implemented in FindSpatiallyVariables(), is to search for features exhibiting spatial patterning in the absence of pre-annotation
d1 <- FindSpatiallyVariableFeatures(d1, assay = "SCT", features = VariableFeatures(d1)[1:1000],
    selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(d1, selection.method = "markvariogram"), 6)
pdf('top_marker_gene_spatial_variable_expression.pdf')
SpatialFeaturePlot(d1, features = top.features, ncol = 3, alpha = c(0.1, 1))
dev.off()


###Subset out anatomical regions
cortex <- subset(d1, idents = c(0))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2


