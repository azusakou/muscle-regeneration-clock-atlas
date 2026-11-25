suppressPackageStartupMessages({library(data.table);library(Matrix)})
set.seed(20261125)
args<-commandArgs(TRUE);root<-if(length(args))args[1]else"."
raw<-Sys.getenv("MUSCLE_LIFESPAN_COUNTS")
stopifnot(dir.exists(raw))
out<-file.path(root,"results");dir.create(out,recursive=TRUE,showWarnings=FALSE)
md<-as.data.table(readRDS(file.path(root,"input/core_cell_metadata.rds")))
md<-md[dataset_id=="GSE232106" & passes_qc_and_singlet==TRUE & !grepl("^Unknown",unified_cell_type)]
md[,unit:=sample_id]
md[age_group=="geriatric",unit:=paste(age_group,time_days,ifelse(replicate%in%c("A","B"),"AB","CD"),sep="::")]
sets<-list(Schwann_identity=c("Sox10","S100b","Plp1","Mpz","Ngfr","Gfap"),
 Schwann_repair=c("Jun","Sox2","Ngfr","Shh"),
 NMJ_postsynaptic=c("Musk","Lrp4","Rapsn","Dok7","Chrne","Ache","Colq"),
 inflammatory_exercise=c("Nos2","Ccl2","Ptgs2"),
 inflammatory_extended=c("Tnf","Il1b","Il6","Ccr2","Ly6c2"),
 restorative=c("Mertk","Axl","Marco","Cd163","Il10","Trem2","Lpl","Apoe"),
 FAP_matrix=c("Col1a1","Col3a1","Fn1","Postn","Ctgf"),
 FAP_receptors=c("Tnfrsf1a","Tnfrsf1b","Tgfbr1","Tgfbr2"),
 macrophage_ligands=c("Tnf","Tgfb1","Vegfa"))
genes<-unique(unlist(sets))
fwrite(rbindlist(lapply(names(sets),function(nm)data.table(module=nm,gene=sets[[nm]]))),file.path(out,"biological_gene_sets.tsv"),sep="\t")
rows<-list()
for(sid in unique(md$sample_id)) {
 checkpoint<-file.path(out,paste0("genes_",sid,".tsv"))
 if(file.exists(checkpoint)){rows[[sid]]<-fread(checkpoint);next}
 message("Scoring ",sid)
 if(grepl("Ger",sid)) {
  prefix<-file.path(raw,paste0(sid,"_soupx"));mf<-paste0(prefix,"_matrix.mtx.gz")
  gf<-paste0(prefix,"_genes.tsv.gz");bf<-paste0(prefix,"_barcodes.tsv.gz")
 }else{
  p<-file.path(raw,"young_soupx",sid);mf<-file.path(p,"matrix.mtx")
  gf<-file.path(p,"genes.tsv");bf<-file.path(p,"barcodes.tsv")
 }
 stopifnot(all(file.exists(c(mf,gf,bf))))
 feat<-fread(gf,header=FALSE);bc<-fread(bf,header=FALSE)[[1]]
 conn<-if(grepl("gz$",mf))gzfile(mf,"rt")else file(mf,"rt")
 mat<-as(readMM(conn),"dgCMatrix");close(conn)
 rownames(mat)<-make.unique(feat[[2]]);colnames(mat)<-paste(sid,bc,sep="__")
 dd<-md[sample_id==sid];ii<-match(dd$source_cell_id,colnames(mat));stopifnot(!anyNA(ii))
 totals<-colSums(mat)[ii];mat<-mat[intersect(genes,rownames(mat)),ii,drop=FALSE]
 expression<-as.matrix(mat);expression<-log1p(sweep(expression,2,pmax(totals,1),"/")*10000)
 a<-rbindlist(lapply(unique(dd$unified_cell_type),function(ct){
  jj<-which(dd$unified_cell_type==ct)
  data.table(sample_id=sid,unit=dd$unit[1],age_class=dd$age_group[1],time_days=dd$time_days[1],
   cell_type=ct,n_cells=length(jj),gene=rownames(mat),mean_logCP10K=rowMeans(expression[,jj,drop=FALSE]),
   fraction_detected=rowMeans(mat[,jj,drop=FALSE]>0))
 }))
 fwrite(a,checkpoint,sep="\t");rows[[sid]]<-a
 rm(mat,expression);gc()
}
g<-rbindlist(rows);fwrite(g,file.path(out,"lineage_gene_expression.tsv.gz"),sep="\t")
units<-g[,.(mean_logCP10K=weighted.mean(mean_logCP10K,n_cells),
 fraction_detected=weighted.mean(fraction_detected,n_cells),n_cells=sum(n_cells)),
 by=.(unit,age_class,time_days,cell_type,gene)]
fwrite(units,file.path(out,"biological_unit_gene_expression.tsv.gz"),sep="\t")
receptors<-units[cell_type=="FAP" & n_cells>=20 & gene%in%sets$FAP_receptors,
 .(n_groups=.N,mean_expression=mean(mean_logCP10K),mean_fraction_detected=mean(fraction_detected)),by=gene]
fwrite(receptors,file.path(out,"FAP_receptor_expression.tsv"),sep="\t")
modules<-rbindlist(lapply(names(sets),function(nm)units[gene%in%sets[[nm]],
 .(score=mean(mean_logCP10K),n_genes=.N,n_cells=unique(n_cells)[1]),
 by=.(unit,age_class,time_days,cell_type)][,module:=nm]))
fwrite(modules,file.path(out,"biological_module_scores.tsv"),sep="\t")

# Published exercise markers are evaluated continuously without assigning
# mutually exclusive M1/M2 labels.
mac<-modules[cell_type=="Macrophage" & n_cells>=20]
fwrite(mac,file.path(out,"macrophage_state_scores.tsv"),sep="\t")
pair<-merge(units[cell_type=="Macrophage" & n_cells>=20 & gene%in%c("Tnf","Tgfb1","Vegfa")],
 modules[cell_type=="FAP" & n_cells>=20 & module=="FAP_matrix",.(unit,matrix_score=score)],by="unit")
assoc<-pair[,{
 fit1<-lm(mean_logCP10K~factor(time_days)+age_class)
 fit2<-lm(matrix_score~factor(time_days)+age_class)
 r<-cor(residuals(fit1),residuals(fit2),method="spearman")
 # Exact permutations preserve injury-time/age strata. Permuting within a
 # stratum leaves the fitted additive covariate projection unchanged.
 permute_indices<-function(v) {
  if(length(v)==1L)return(matrix(v,nrow=1L))
  do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],permute_indices(v[-i]))))
 }
 strata<-split(seq_len(.N),interaction(time_days,age_class,drop=TRUE))
 pp<-lapply(strata,permute_indices)
 gg<-expand.grid(lapply(pp,function(a)seq_len(nrow(a))))
 r1<-rank(residuals(fit1));r2<-rank(residuals(fit2))
 null<-vapply(seq_len(nrow(gg)),function(i){
  idx<-seq_len(.N)
  for(j in seq_along(strata))idx[strata[[j]]]<-pp[[j]][gg[i,j],]
  cor(r1,r2[idx])
 },0.0)
 exact_p<-mean(abs(null)>=abs(r)-1e-12)
 .(n_units=.N,raw_rho=cor(mean_logCP10K,matrix_score,method="spearman"),
 adjusted_rho=r,adjusted_p=exact_p,n_permutations=length(null))
},by=gene]
assoc[,BH_q:=p.adjust(adjusted_p,"BH")]
fwrite(assoc,file.path(out,"macrophage_FAP_associations.tsv"),sep="\t")
fwrite(pair,file.path(out,"macrophage_FAP_pairs.tsv"),sep="\t")
print(assoc)
print(modules[module%in%c("Schwann_identity","NMJ_postsynaptic"),
 .(n_units=uniqueN(unit),mean_score=mean(score)),by=.(module,cell_type)])
