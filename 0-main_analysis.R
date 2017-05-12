library(rtracklayer) # 导入导出各类格式
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
setwd("../data")

#//////////////////////////////////////////
# 利用良定的区域
#/////////////////////////////////////////
refine.region = keepStandardChromosomes(import.bed("refine_sort_hg19.bed"))
NCBI.style =  mapSeqlevels(seqlevels(refine.region), "NCBI")
refine.region = renameSeqlevels(refine.region, NCBI.style)
raw.data = import("R_raw.bed", format = "BED") 

cds.final = refine.region[refine.region$name == "CDS"]
intergenic.final = refine.region[refine.region$name == "Intergenic"]
promoter.final = refine.region[refine.region$name == "Promoter"]
enhancer.final = refine.region[refine.region$name == "Enhancer"]
utr3.final = refine.region[refine.region$name == "UTR3"]
utr5.final = refine.region[refine.region$name == "UTR5"]
intron.final = refine.region[refine.region$name == "Intron"]

#/////////////////
# 获得全局突变率
#////////////////
SAMPLE_SIZE = 2954
FINALS = grep(".final", ls(), value = T)
for(i in FINALS){
  assign(gsub("final", "rate", i), 
         sum(countOverlaps(raw.data, reduce(get(i)))/
               (SAMPLE_SIZE * sum(width(reduce(get(i)))))))
}
pdf("rate.pdf", width = 20, height = 10)
barplot(sort(unlist(mget(grep("rate", ls(), value = T)))))
dev.off()
#/////////////////
# 读取簇信信息
#/////////////////
sample.global.prop = read.table("sample_prop.tsv")
names(sample.global.prop) = c("sample.id", "global.prop")
merged.snv = read_tsv("for_R_analysis.bed", col_types = "ciiiic", 
                      col_names = c("chr", "start", "end", "mutation.number", 
                                    "mutation.sample", "sample.id")) # many metadata
sample.id = merged.snv$sample.id
prop = vector(mode = "double", length = length(sample.id))
for (i in seq_along(sample.id)){
  sps = unlist(strsplit(sample.id[i], split = ","))
  prop[i] = mean(sample.global.prop[sample.global.prop$sample.id 
                                    %in% sps, "global.prop"])
}
merged.snv$prop = prop

merged.snv$failure.number = merged.snv$end - merged.snv$start - 
  merged.snv$mutation.number

merged.snv$sig = vector(mode = "numeric", length = dim(merged.snv)[1])
merged.snv$sig = pnbinom(merged.snv$failure.number, merged.snv$mutation.number, 
                         merged.snv$prop) # vectorization

merged.snv$adjust.sig = p.adjust(merged.snv$sig, method = "BH")
save(merged.snv, file = "merged.snv.new.RData") # save the file 

all.cluster = 
  makeGRangesFromDataFrame(merged.snv[, c("chr", "start", "end", 
                                          "adjust.sig", "mutation.sample", "sample.id")], 
                           keep.extra.columns = T, ignore.strand = T, 
                           seqnames.field = "chr", starts.in.df.are.0based = T) 
all.cluster = all.cluster[width(all.cluster) >= 50,]
#///////////////////////////
# 确定不同簇所属的元件类型
#//////////////////////////
load("regional-mutation-rate.RData") # 各簇突变信息
cluster.enhancer = subsetByOverlaps(all.cluster, enhancer.final)
next.cluster = all.cluster[!overlapsAny(all.cluster, enhancer.final)]

cluster.promoter = subsetByOverlaps(next.cluster, promoter.final)
next.cluster = next.cluster[!overlapsAny(next.cluster, promoter.final)]

cluster.utr5 = subsetByOverlaps(next.cluster, utr5.final)
next.cluster = next.cluster[!overlapsAny(next.cluster, utr5.final)]

cluster.utr3 = subsetByOverlaps(next.cluster, utr3.final)
next.cluster = next.cluster[!overlapsAny(next.cluster, utr3.final)]

cluster.intron = subsetByOverlaps(next.cluster, intron.final)
next.cluster = next.cluster[!overlapsAny(next.cluster, intron.final)]

cluster.intergenic = subsetByOverlaps(next.cluster, intergenic.final)
##//////////////////////////////////////////////////////////
# 全局突变率假设检
#//////////////////////////////////////////////////////////
GlobalSigTest = function(specific.cluster, cluster.rate){
  specific.cluster$global.rate = 1 - (1-cluster.rate)^(width(specific.cluster))
  specific.cluster$global.p = pbinom(q = specific.cluster$mutation.sample, 
                                     size = SAMPLE_SIZE, lower.tail = F,
                                     prob = specific.cluster$global.rate)
  mcols(specific.cluster) = mcols(specific.cluster)[c(1, 2, 3, 5)]
  return(specific.cluster)
}
##//////////////////////////////////////////////////////////
# 局部突变率假设检
#//////////////////////////////////////////////////////////
normalized.mut.rate = mut.rate/SAMPLE_SIZE
ReginalSigTest = function(cluster.GR, mutation.rate){
  # original cluster and normalized mutation rate
  cluster.GR$reginal.rate = 1 - (1-mutation.rate)^(width(all.cluster))
  cluster.GR$reginal.p = pbinom(q = cluster.GR$mutation.sample, 
                                size = SAMPLE_SIZE, lower.tail = F,
                                prob = cluster.GR$reginal.rate)
  mcols(cluster.GR) = mcols(cluster.GR)[c(1, 2, 3, 5)]
  return(cluster.GR)
}
#//////////////////////////////////////////////////////////
# 进行检验和结果合并
# "reginal.test.all.cluster"
# 1. only preserve "adjust.sig", "mutation.sample", "reginal.p".
# 2. combine all the global test and new column "element".
# 3. combine global test and reginal test.
# 4. do BH adjustment
# 5. convert to GRange: "g.r" is the final GRanges!
# 6. get sigs p<0.001
#//////////////////////////////////////////////////////////
reginal.test.all = ReginalSigTest(all.cluster, normalized.mut.rate)

global.test.enhancer = GlobalSigTest(cluster.enhancer, enhancer.rate)
global.test.promoter = GlobalSigTest(cluster.promoter, promoter.rate)
global.test.intergenic = GlobalSigTest(cluster.intergenic, intergenic.rate)
global.test.intron = GlobalSigTest(cluster.intron, intron.rate)
global.test.utr5 = GlobalSigTest(cluster.utr5, utr5.rate)
global.test.utr3 = GlobalSigTest(cluster.utr3, utr3.rate)

global.test.enhancer$element = "enhancer"
global.test.promoter$element = "promoter"
global.test.intergenic$element = "intergenic"
global.test.intron$element = "intron"
global.test.utr3$element = "utr3"
global.test.utr5$element = "utr5"

global.test.enhancer = as.data.frame(global.test.enhancer)
global.test.intergenic = as.data.frame(global.test.intergenic)
global.test.promoter = as.data.frame(global.test.promoter)
global.test.intron = as.data.frame(global.test.intron)
global.test.utr3 = as.data.frame(global.test.utr3)
global.test.utr5 = as.data.frame(global.test.utr5)

global.test.all = Reduce(rbind, mget(ls(pattern = "global.test.*")))
reginal.test.all = as.data.frame(reginal.test.all)

global.reginal = inner_join(global.test.all, reginal.test.all)
global.reginal$adjust.global = p.adjust(global.reginal$global.p, method = "BH")
global.reginal$adjust.reginal = p.adjust(global.reginal$reginal.p, method = "BH")
global.reginal = global.reginal[c("seqnames", "start", "end", "strand", 
                                  "adjust.sig", "mutation.sample", "element", 
                                  "adjust.global", "adjust.reginal", "sample.id")]
names(global.reginal)[7] = c("type")
g.r = 
  makeGRangesFromDataFrame(global.reginal, keep.extra.columns = T,
                           starts.in.df.are.0based = T) # from bed

g.r.sig = g.r[g.r$adjust.global < 0.001 & g.r$adjust.reginal < 0.001 & 
                g.r$adjust.reginal != 0, ] # we have del 0, wrong or right?
g.r.sig = as.data.frame(g.r.sig)
#////////////////
# 进行可视化
#////////////////
reverselog_trans = function(base = exp(1)){
  trans = function(x) -log(x, base)
  inv = function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), domain = c(1e-100, Inf))
}
ggplot(g.r.sig,
       aes(x = adjust.global, y = adjust.reginal, shape = type)) + 
  geom_point(size = 2) + 
  scale_x_continuous(trans = reverselog_trans(10), 
                     labels = trans_format("log10", math_format(10^.x)), name = "Regional recurrence, global model (adjusted P value)") + 
  scale_y_continuous(trans = reverselog_trans(10), 
                     labels = trans_format("log10", math_format(10^.x)), name = "Regional recurrence, reginal model (adjusted P value)")
