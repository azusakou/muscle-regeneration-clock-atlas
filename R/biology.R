lifespan_validation <- function() {
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

}

# Independent validation cohorts are not used to select the primary clock.
independent_validation <- function(root) {
 suppressPackageStartupMessages({library(data.table);library(Matrix);library(SingleCellExperiment);library(DropletUtils);library(scDblFinder);library(BiocParallel)})
 set.seed(20261125)
 out<-file.path(root,'results');dir.create(out,recursive=TRUE,showWarnings=FALSE)
 raw<-Sys.getenv('MUSCLE_AGE_VALIDATION',file.path(root,'data/GSE180225'))
 ids<-c(paste0('GSM545628',1:8),paste0('GSM75964',58:62))
 
 # Author analysis code Seurat_1.Rmd identifies A_day4PI/A_day7PI as young
 # and OU_R1/OU_R2 as aged baseline, correcting conflicting GEO titles.
 # Labels are fixed from provenance before any age comparison.
 design<-data.table(sample_id=ids,age_class=c('young','young','young','young','old','old','old','old','young','young','old','old','old'),time_days=c(0,0,4,7,0,0,4,7,4,7,4,4,7))
 markers<-list(MuSC_myogenic=c('Pax7','Calcr','Sdc4','Vcam1','Myod1','Myog','Des'),myonucleus=c('Ttn','Myh1','Myh2','Myh4','Ckm','Acta1'),FAP=c('Pdgfra','Col1a1','Col3a1','Dcn','Lum'),macrophage=c('Ptprc','Lyz2','Adgre1','Cd68','C1qa','C1qb'),neutrophil=c('S100a8','S100a9','Retnlg','Camp','Ly6g'),T_NK=c('Cd3d','Cd3e','Trbc1','Nkg7','Klrd1'),endothelial=c('Pecam1','Kdr','Emcn','Cdh5','Klf2'),pericyte=c('Rgs5','Cspg4','Pdgfrb','Notch3','Des'),Schwann=c('Sox10','S100b','Plp1','Mpz','Ngfr','Gfap'))
 programs<-list(Activation=c('Fos','Jun','Junb','Egr1','Myod1'),Proliferation=c('Mki67','Top2a','Ccnb1','Pcna','Ube2c'),Commitment=c('Myod1','Myf5','Des','Cdkn1c'),Fusion=c('Myog','Mymk','Mymx','Tnnt1','Ckm'),Remodeling=c('Col1a1','Col3a1','Fn1','Mmp14','Tgfb1'),Resolution=c('Mrc1','C1qa','C1qb','Il10','Apoe'))
 all<-list();qcs<-list()
 for(sid in ids) {
  checkpoint<-file.path(out,paste0('age_alias_validation_',sid,'.rds'))
  if(file.exists(checkpoint)){z<-readRDS(checkpoint);all[[sid]]<-z$scores;qcs[[sid]]<-z$qc;next}
  f<-list.files(raw,pattern=paste0('^',sid,'_.*\\.h5$'),full.names=TRUE);stopifnot(length(f)==1)
  message('Independent age: ',sid)
  pointer<-rhdf5::h5read(f,"matrix/indptr");shape<-rhdf5::h5read(f,"matrix/shape")
  candidate<-which(diff(pointer)>=200)
  m<-new("dgCMatrix",i=as.integer(rhdf5::h5read(f,"matrix/indices")),p=as.integer(pointer),x=as.numeric(rhdf5::h5read(f,"matrix/data")),Dim=as.integer(shape))
  m<-m[,candidate,drop=FALSE];symbols<-as.character(rhdf5::h5read(f,"matrix/features/name"));if(!'Mymk'%in%symbols)symbols[symbols=='Tmem8c']<-'Mymk';if(!'Mymx'%in%symbols)symbols[symbols=='Gm7325']<-'Mymx';rownames(m)<-make.unique(symbols)
  s<-SingleCellExperiment(list(counts=m));umi<-Matrix::colSums(m);det<-Matrix::colSums(m>0)
  mt<-if(any(grepl('^mt-',symbols,ignore.case=TRUE)))Matrix::colSums(m[grepl('^mt-',symbols,ignore.case=TRUE),,drop=FALSE])/pmax(umi,1)*100 else rep(0,length(umi))
  keep<-umi>=500&umi<=50000&det>=200&det<=6000&mt<=20
  n_raw<-shape[2];s<-s[,keep];pre_doublet<-ncol(s)
  message(sid,': ',pre_doublet,' nuclei before doublet removal');s<-scDblFinder(s,BPPARAM=SerialParam(),verbose=FALSE)
  s<-s[,s$scDblFinder.class=='singlet'];m<-counts(s)
  lc<-m;lc@x<-log1p(lc@x/rep(Matrix::colSums(m),diff(lc@p))*10000)
  score<-function(gs){hit<-intersect(gs,rownames(lc));if(!length(hit))return(rep(0,ncol(lc)));Matrix::colMeans(lc[hit,,drop=FALSE])}
  ms<-sapply(markers,score);top<-max.col(ms,ties.method='first');ord<-t(apply(ms,1,sort,decreasing=TRUE))
  ct<-names(markers)[top];ct[ord[,1]-ord[,2]<.08|ord[,1]<.02]<-'Unknown_or_ambiguous'
  d<-as.data.table(sapply(programs,score));d[,cell_type:=ct]
  sm<-d[cell_type!='Unknown_or_ambiguous',c(list(n_cells=.N),lapply(.SD,mean)),by=cell_type,.SDcols=names(programs)]
  sm[,sample_id:=sid]
  qc<-design[sample_id==sid];qc[,`:=`(raw_barcodes=n_raw,core_qc=pre_doublet,singlets=ncol(s),annotated=sum(ct!='Unknown_or_ambiguous'),file=basename(f))]
  saveRDS(list(scores=sm,qc=qc),checkpoint);all[[sid]]<-sm;qcs[[sid]]<-qc
  rm(s,m,lc,ms,d);gc()
 }
 d<-merge(rbindlist(all),design,by='sample_id');qc<-rbindlist(qcs);qc[,c('age_class','time_days'):=NULL];qc<-merge(qc,design,by='sample_id')
 fwrite(d,file.path(out,'independent_age_lineage_scores.tsv'),sep='\t');fwrite(qc,file.path(out,'independent_age_QC.tsv'),sep='\t')
 evaluate_independent_age(d,root,names(programs))
}
evaluate_independent_age <- function(d,root,programs=c('Activation','Proliferation','Commitment','Fusion','Remodeling','Resolution')) {
 suppressPackageStartupMessages(library(data.table));out<-file.path(root,'results');d<-copy(d[n_cells>=20]);design<-unique(d[,.(sample_id,age_class,time_days)])
 soft<-function(z){z<-sweep(z,1,apply(z,1,max),'-')/.8;p<-exp(z);as.numeric((p/rowSums(p))%*%c(.12,.28,.45,.64,.80,.94))}
 fold<-function(train,test){
  mu<-colMeans(train[,..programs]);ss<-sapply(train[,..programs],sd);ss[!is.finite(ss)|ss<1e-8]<-1
  project<-function(a){z<-sweep(sweep(as.matrix(a[,..programs]),2,mu,'-'),2,ss,'/');a<-copy(a);a[,`:=`(phase=soft(z),a=log1p(time_days)/log(8))];a}
  tr<-project(train);te<-project(test);te[,expected:=NA_real_]
  for(ct in unique(te$cell_type)){
   a<-tr[cell_type==ct];ii<-which(te$cell_type==ct)
   if(nrow(a)<3||uniqueN(a$time_days)<2)next
   te$expected[ii]<-as.numeric(predict(lm(phase~a,data=a),te[ii]))
  };te[,deviation:=phase-expected];te
 }
 score<-function(labels){
  dd<-copy(d);dd[,age_class:=unname(labels[sample_id])];young<-dd[age_class=='young']
  old<-fold(young,dd[age_class=='old']);yy<-rbindlist(lapply(unique(young$sample_id),function(id)fold(young[sample_id!=id],young[sample_id==id])))
  dev<-rbindlist(list(yy,old));sig<-max(.05,sd(yy$deviation,na.rm=TRUE));if(!is.finite(sig))sig<-.05
  res<-rbindlist(lapply(c('all_lineages','common_four'),function(panel){
   a<-dev[is.finite(deviation)];if(panel=='common_four')a<-a[cell_type%in%c('FAP','macrophage','endothelial','myonucleus')]
   z<-a[,.(RSI=if(.N>=4)exp(-sum((deviation-mean(deviation))^2)/.N/sig^2)else NA_real_,n_lineages=.N,mean_phase=mean(phase)),by=.(sample_id,age_class,time_days)]
   z[,panel:=panel];z
  }),fill=TRUE)
  list(samples=res,lineages=dev)
 }
 labels<-setNames(design$age_class,design$sample_id);obs<-score(labels)
 effect<-function(s,pn,tt){a<-s[panel==pn & time_days==tt & is.finite(RSI)];if(!all(c('young','old')%in%a$age_class))return(NA_real_);mean(a[age_class=='old',RSI])-mean(a[age_class=='young',RSI])}
 tests<-list();times<-sort(unique(design$time_days))
 for(tt in times){
  id<-design[time_days==tt,sample_id];ny<-design[time_days==tt&age_class=='young',.N];pp<-combn(id,ny,simplify=FALSE)
  null<-sapply(pp,function(ys){lab<-labels;lab[id]<-ifelse(id%in%ys,'young','old');z<-score(lab)$samples;c(all_lineages=effect(z,'all_lineages',tt),common_four=effect(z,'common_four',tt))})
  for(pn in c('all_lineages','common_four')){
   e<-effect(obs$samples,pn,tt);vals<-null[pn,];a<-obs$samples[panel==pn&time_days==tt&is.finite(RSI)]
   tests[[paste(pn,tt)]]<-data.table(panel=pn,time_days=tt,n_young=sum(a$age_class=='young'),n_old=sum(a$age_class=='old'),effect=e,exact_refit_p=if(is.finite(e)&&all(is.finite(vals)))mean(abs(vals)>=abs(e)-1e-12)else NA_real_,n_allocations=length(vals),n_invalid=sum(!is.finite(vals)))
  }
 }
 tests<-rbindlist(tests);tests[,holm_p:=p.adjust(exact_refit_p,'holm',n=3),by=panel]
 # A separate age-label-independent myonuclear phase projection uses all 13
 # libraries, without a cross-lineage eligibility requirement.
 mu<-colMeans(d[,..programs]);ss<-sapply(d[,..programs],sd);ss[!is.finite(ss)|ss<1e-8]<-1
 z<-sweep(sweep(as.matrix(d[,..programs]),2,mu,'-'),2,ss,'/');d[,pooled_phase:=soft(z)]
 myo<-d[cell_type=='myonucleus'];myo[,age_class:=factor(age_class,c('young','old'))]
 coef_age<-function(ag){X<-model.matrix(~factor(myo$time_days)+factor(ag,levels=c('young','old')));tail(lm.fit(X,myo$pooled_phase)$coefficients,1)}
 actual<-unname(coef_age(myo$age_class));strata<-split(seq_len(nrow(myo)),myo$time_days);alloc<-lapply(strata,function(ii)combn(ii,sum(myo$age_class[ii]=='young'),simplify=FALSE));grid<-expand.grid(lapply(alloc,seq_along))
 null<-vapply(seq_len(nrow(grid)),function(k){ag<-rep('old',nrow(myo));for(j in seq_along(alloc))ag[alloc[[j]][[grid[[j]][k]]]]<-'young';unname(coef_age(ag))},0.0)
 phase_test<-data.table(cell_type='myonucleus',effect=actual,exact_p=mean(abs(null)>=abs(actual)-1e-12),n_libraries=nrow(myo),n_allocations=length(null))
 fwrite(obs$lineages,file.path(out,'independent_age_phase_deviations.tsv'),sep='\t');fwrite(obs$samples,file.path(out,'independent_age_RSI.tsv'),sep='\t');fwrite(tests,file.path(out,'independent_age_tests.tsv'),sep='\t');fwrite(myo,file.path(out,'independent_age_myonuclear_phase.tsv'),sep='\t');fwrite(phase_test,file.path(out,'independent_age_myonuclear_test.tsv'),sep='\t');print(tests);print(phase_test)
}

denervation_validation <- function(root) {
 suppressPackageStartupMessages({library(SingleCellExperiment);library(Matrix);library(scuttle);library(scDblFinder);library(BiocParallel);library(data.table);library(edgeR)})
 set.seed(20261125);out<-file.path(root,'results')
 checkpoint<-Sys.getenv('MUSCLE_DENERVATION_QC')
 if(nzchar(checkpoint))s<-readRDS(checkpoint) else {
  suppressPackageStartupMessages(library(Seurat))
  raw<-Sys.getenv('MUSCLE_DENERVATION_RDS');stopifnot(file.exists(raw))
  con<-gzcon(gzfile(raw,open="rb"));obj<-readRDS(con);close(con)
  s<-as.SingleCellExperiment(obj);rm(obj);gc()
  s$sample_id<-as.character(s$stim)
  q<-perCellQCMetrics(s,subsets=list(Mito=grepl("^mt-",rownames(s))))
  ls<-log10(q$sum+1);lg<-log10(q$detected+1)
  keep<-ls>=median(ls)-3*mad(ls)&ls<=median(ls)+5*mad(ls)&lg>=median(lg)-3*mad(lg)&lg<=median(lg)+5*mad(lg)&q$subsets_Mito_percent<=median(q$subsets_Mito_percent)+3*mad(q$subsets_Mito_percent)
  s<-s[,keep];set.seed(20260715)
  s<-scDblFinder(s,samples="sample_id",BPPARAM=SerialParam(),verbose=FALSE);s$passes_QC<-s$scDblFinder.class=="singlet"
 }
 s<-s[,as.logical(s$passes_QC)];m<-counts(s);ct<-as.character(s$ident)
 body<-ct%in%c('TypeI','TypeIIA','TypeIIx','TypeIIB');sid<-as.character(s$sample_id)
 sample_ids<-sort(unique(sid));group<-ifelse(grepl('^N',sample_ids),'denervated','control')
 mat<-sapply(sample_ids,function(id)Matrix::rowSums(m[,sid==id&body,drop=FALSE]));rownames(mat)<-rownames(m)
 meta<-data.table(sample_id=sample_ids,condition=group,n_nuclei=vapply(sample_ids,function(id)sum(sid==id&body),0L))
 y<-DGEList(mat);design<-model.matrix(~factor(group,levels=c('control','denervated')))
 y<-y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE)
 fit<-glmQLFit(y,design,robust=TRUE);test<-glmQLFTest(fit,coef=2);tab<-as.data.table(topTags(test,n=Inf)$table,keep.rownames='gene')
 sets<-list(denervation_response=c('Chrna1','Chrnb1','Chrnd','Chrng','Musk','Ncam1','Myog','Hdac4','Runx1'),postsynaptic_NMJ=c('Musk','Lrp4','Rapsn','Dok7','Chrne','Ache','Colq'))
 logcpm<-cpm(y,log=TRUE,prior.count=2);zs<-t(scale(t(logcpm)));zs[!is.finite(zs)]<-0
 idx<-lapply(sets,function(g)which(rownames(y)%in%g));cam<-as.data.table(camera(y,index=idx,design=design,contrast=2,inter.gene.cor=.01,sort=FALSE),keep.rownames='module');cam[,BH_q:=p.adjust(PValue,'BH')]
 scores<-rbindlist(lapply(names(sets),function(nm){i<-idx[[nm]];cbind(meta,data.table(module=nm,score=colMeans(zs[i,,drop=FALSE]),mean_logCPM=colMeans(logcpm[i,,drop=FALSE])))}))
 fwrite(data.table(gene=rownames(logcpm),logcpm)[gene%in%unlist(sets)],file.path(out,'denervation_sample_gene_expression.tsv'),sep='\t');
 fwrite(tab[gene%in%unlist(sets)],file.path(out,'denervation_gene_tests.tsv'),sep='\t');fwrite(cam,file.path(out,'denervation_pathway_tests.tsv'),sep='\t');fwrite(scores,file.path(out,'denervation_sample_scores.tsv'),sep='\t')
 fwrite(data.table(gene=rownames(mat),mat),file.path(out,'denervation_body_pseudobulk.tsv.gz'),sep='\t')
 fwrite(meta,file.path(out,'denervation_sample_design.tsv'),sep='\t')

 # A fixed four-lineage panel separates the synchrony contrast from recovery
 # weights and lineage-composition changes in the denervation experiment.
 programs<-list(Activation=c("Fos","Jun","Junb","Egr1","Myod1"),Proliferation=c("Mki67","Top2a","Ccnb1","Pcna","Ube2c"),Commitment=c("Myod1","Myf5","Des","Cdkn1c"),Fusion=c("Myog","Mymk","Mymx","Tnnt1","Ckm"),Remodeling=c("Col1a1","Col3a1","Fn1","Mmp14","Tgfb1"),Resolution=c("Mrc1","C1qa","C1qb","Il10","Apoe"))
 broad<-ifelse(body,"myonucleus",ifelse(ct=="FAPs","FAP",ifelse(ct=="Macrophage","macrophage",ifelse(ct=="Endothelial","endothelial","excluded"))))
 lc<-m;lc@x<-log1p(lc@x/rep(Matrix::colSums(m),diff(lc@p))*10000)
 sc<-as.data.table(sapply(programs,function(g)Matrix::colMeans(lc[intersect(g,rownames(lc)),,drop=FALSE])))
 sc[,`:=`(sample_id=sid,cell_type=broad)]
 phase_input<-sc[cell_type!="excluded",c(list(n_cells=.N),lapply(.SD,mean)),by=.(sample_id,cell_type),.SDcols=names(programs)]
 phase_input<-merge(phase_input,meta[,.(sample_id,condition)],by="sample_id")
 fwrite(phase_input,file.path(out,"denervation_lineage_programs.tsv"),sep="\t")
 denervation_RSI(phase_input,root)
 print(cam);print(tab[gene%in%unlist(sets),.(gene,logFC,FDR)])
}

denervation_RSI <- function(d,root) {
 suppressPackageStartupMessages(library(data.table));d<-copy(d[n_cells>=20]);pcols<-c("Activation","Proliferation","Commitment","Fusion","Remodeling","Resolution")
 eligible<-d[,.N,by=sample_id][N==4,sample_id];d<-d[sample_id%in%eligible]
 stopifnot(uniqueN(d$sample_id)==6)
 project<-function(train,test){mu<-colMeans(train[,..pcols]);ss<-sapply(train[,..pcols],sd);ss[!is.finite(ss)|ss<1e-8]<-1
  z<-sweep(sweep(as.matrix(test[,..pcols]),2,mu,"-"),2,ss,"/");z<-sweep(z,1,apply(z,1,max),"-")/.8;pr<-exp(z)
  a<-copy(test);a[,phase:=as.numeric((pr/rowSums(pr))%*%c(.12,.28,.45,.64,.80,.94))];a}
 calc<-function(ref_ids){ref<-d[sample_id%in%ref_ids];other<-d[!sample_id%in%ref_ids]
  rr<-project(ref,ref);mu<-rr[,.(expected=mean(phase)),by=cell_type]
  oo<-merge(project(ref,other),mu,by="cell_type");oo[,deviation:=phase-expected]
  yy<-rbindlist(lapply(ref_ids,function(id){tr<-ref[sample_id!=id];m<-project(tr,tr)[,.(expected=mean(phase)),by=cell_type];a<-merge(project(tr,ref[sample_id==id]),m,by="cell_type");a[,deviation:=phase-expected];a}))
  sigma<-max(.05,sd(yy$deviation));dev<-rbindlist(list(yy,oo))
  a<-dev[,.(RSI=exp(-var(deviation)*3/4/sigma^2)),by=sample_id];a[,reference:=sample_id%in%ref_ids]
  list(samples=a,effect=mean(a[reference==FALSE,RSI])-mean(a[reference==TRUE,RSI]))}
 ids<-sort(unique(d$sample_id));ref<-unique(d[condition=="control",sample_id]);observed<-calc(ref);ii<-combn(ids,3);null<-apply(ii,2,function(k)calc(k)$effect)
 fwrite(merge(observed$samples,unique(d[,.(sample_id,condition)]),by="sample_id"),file.path(root,"results/denervation_RSI.tsv"),sep="\t")
 fwrite(data.table(effect=observed$effect,exact_refit_p=mean(abs(null)>=abs(observed$effect)-1e-12),n_allocations=length(null),n_control=3,n_denervated=3,n_lineages=4),file.path(root,"results/denervation_RSI_test.tsv"),sep="\t")
}

depletion_gene_validation <- function(root) {
 suppressPackageStartupMessages({library(data.table);library(Matrix)})
 raw<-Sys.getenv('MUSCLE_DEPLETION_COUNTS');mf<-Sys.getenv('MUSCLE_DEPLETION_METADATA','results/GSE291782_failure/cell_metadata.tsv.gz')
 md<-fread(mf);out<-file.path(root,'results')
 sets<-list(FAP_matrix=c('Col1a1','Col3a1','Fn1','Postn','Ctgf'),TNF_response=c('Nfkbia','Tnfaip3','Relb','Birc3','Icam1'),TGFbeta_response=c('Serpine1','Ctgf','Smad7','Tagln'),macrophage_ligands=c('Tnf','Tgfb1','Vegfa'))
 if(!'cell_id'%in%names(md))stop('Cell identifiers missing')
 rows<-list()
 for(sid in unique(md$sample_id)) {
  p<-file.path(raw,sid);feat<-fread(paste0(p,'_features.tsv.gz'),header=FALSE);bc<-fread(paste0(p,'_barcodes.tsv.gz'),header=FALSE)[[1]]
  m<-as(readMM(gzfile(paste0(p,'_matrix.mtx.gz'))),'dgCMatrix');rownames(m)<-make.unique(feat[[2]]);colnames(m)<-paste(sid,bc,sep='::')
  dd<-md[sample_id==sid & passes_qc_and_singlet==TRUE & coarse_cell_type%in%c('FAP','macrophage')]
  ii<-match(dd$cell_id,colnames(m));stopifnot(!anyNA(ii));tot<-colSums(m)[ii];sm<-m[intersect(unique(unlist(sets)),rownames(m)),ii,drop=FALSE]
  lc<-log1p(sweep(as.matrix(sm),2,pmax(tot,1),'/')*10000)
  rows[[sid]]<-rbindlist(lapply(unique(dd$coarse_cell_type),function(ct){jj<-which(dd$coarse_cell_type==ct);data.table(sample_id=sid,condition=dd$condition[1],cell_type=ct,n_cells=length(jj),gene=rownames(sm),mean_logCP10K=rowMeans(lc[,jj,drop=FALSE]),fraction_detected=rowMeans(sm[,jj,drop=FALSE]>0))}))
 }
 g<-rbindlist(rows);fwrite(g,file.path(out,'depletion_gene_expression.tsv'),sep='\t')
 z<-rbindlist(lapply(names(sets),function(nm)g[gene%in%sets[[nm]],.(score=mean(mean_logCP10K),n_cells=unique(n_cells)[1]),by=.(sample_id,condition,cell_type)][,module:=nm]))
 z<-z[n_cells>=20 & (cell_type=='FAP'&module!='macrophage_ligands')]
 tests<-z[,{v<-score;ids<-which(condition=='clodronate');ii<-combn(seq_along(v),length(ids));obs<-mean(v[ids])-mean(v[-ids]);null<-apply(ii,2,function(k)mean(v[k])-mean(v[-k]));list(effect=obs,exact_p=mean(abs(null)>=abs(obs)-1e-12),n_sham=sum(condition=='sham'),n_clodronate=sum(condition=='clodronate'))},by=.(cell_type,module)]
 tests[,BH_q:=p.adjust(exact_p,'BH')]
 fwrite(z,file.path(out,'depletion_programme_scores.tsv'),sep='\t');fwrite(tests,file.path(out,'depletion_programme_tests.tsv'),sep='\t');print(tests)
}

harmonized_age_validation <- function(root) {
 suppressPackageStartupMessages({library(data.table);library(Matrix)})
 out<-file.path(root,'results');rawroot<-Sys.getenv('MUSCLE_DATA_ROOT','data');lifespan<-Sys.getenv('MUSCLE_LIFESPAN_COUNTS');historic<-Sys.getenv('MUSCLE_NOTEXIN_COUNTS')
 md<-as.data.table(readRDS(file.path(root,'input/core_cell_metadata.rds')))[dataset_id%in%c('GSE232106','GSE143437')&passes_qc_and_singlet==TRUE]
 map<-c(Myogenic_progenitor='MuSC_myogenic',Myonucleus='myonucleus',FAP='FAP',Macrophage='macrophage',Neutrophil='neutrophil',Endothelial='endothelial',Pericyte='pericyte',T_cell='T_NK',B_cell='other_immune',Schwann='Schwann')
 md[,cell_type:=unname(map[unified_cell_type])];md<-md[!is.na(cell_type)]
 programs<-list(Activation=c('Fos','Jun','Junb','Egr1','Myod1'),Proliferation=c('Mki67','Top2a','Ccnb1','Pcna','Ube2c'),Commitment=c('Myod1','Myf5','Des','Cdkn1c'),Fusion=c('Myog','Mymk','Mymx','Tnnt1','Ckm'),Remodeling=c('Col1a1','Col3a1','Fn1','Mmp14','Tgfb1'),Resolution=c('Mrc1','C1qa','C1qb','Il10','Apoe'))

 # Use the same available genes in all four NTX accessions.
 probe<-list.files(lifespan,pattern="_genes.tsv.gz$",full.names=TRUE)[1]
 canonical<-function(g){if(!'Mymk'%in%g)g[g=='Tmem8c']<-'Mymk';if(!'Mymx'%in%g)g[g=='Gm7325']<-'Mymx';g};g1<-canonical(fread(probe,header=FALSE)[[2]])
 gf<-list.files(file.path(rawroot,"extracted/GSE159500"),pattern="_soupx_genes.tsv.gz$",full.names=TRUE)[1];g2<-fread(gf,header=FALSE)[[2]]
 zf<-list.files(file.path(rawroot,"extracted/GSE162172"),pattern="_soupx.zip$",full.names=TRUE)[1];member<-unzip(zf,list=TRUE)$Name;member<-member[grepl("/genes.tsv$",member)]
 g3<-read.delim(unz(zf,member),header=FALSE)[[1]]
 g4<-fread(file.path(historic,"GSE143437_DeMicheli_MuSCatlas_rawdata.txt.gz"),select=1)[[1]]
 common<-Reduce(intersect,lapply(list(g1,g2,g3,g4),canonical))
 coverage<-rbindlist(lapply(names(programs),function(nm)data.table(programme=nm,gene=programs[[nm]],shared=programs[[nm]]%in%common)))
 fwrite(coverage,file.path(out,"harmonized_programme_coverage.tsv"),sep="\t")
 programs<-lapply(programs,function(g)intersect(g,common));stopifnot(all(lengths(programs)>=3))
 genes<-unique(unlist(programs));res<-list()
 summarize<-function(m,dd,ids,totals=NULL){
  rownames(m)<-make.unique(canonical(rownames(m)))
  ii<-match(dd$source_cell_id,ids);stopifnot(!anyNA(ii));tot<-if(is.null(totals))Matrix::colSums(m)[ii]else totals[ii]
  hit<-intersect(genes,rownames(m));lc<-log1p(sweep(as.matrix(m[hit,ii,drop=FALSE]),2,pmax(tot,1),'/')*10000)
  ss<-as.data.table(sapply(programs,function(g)colMeans(lc[intersect(g,rownames(lc)),,drop=FALSE])))
  ss[,cell_type:=dd$cell_type]
  z<-ss[,c(list(n_cells=.N),lapply(.SD,mean)),by=cell_type,.SDcols=names(programs)]
  z[,`:=`(sample_id=dd$sample_id[1],dataset_id=dd$dataset_id[1],time_days=dd$time_days[1],age_class=dd$age_group[1])];z
 }
 for(sid in unique(md[dataset_id=='GSE232106',sample_id])) {
  ck<-file.path(out,paste0('harmonized_alias_',sid,'.rds'));if(file.exists(ck)){res[[sid]]<-readRDS(ck);next}
  message('Common programmes: ',sid)
  dd<-md[sample_id==sid];if(grepl('Ger',sid)){
   p<-file.path(lifespan,paste0(sid,'_soupx'));mf<-paste0(p,'_matrix.mtx.gz');gf<-paste0(p,'_genes.tsv.gz');bf<-paste0(p,'_barcodes.tsv.gz')
  }else {p<-file.path(lifespan,'young_soupx',sid);mf<-file.path(p,'matrix.mtx');gf<-file.path(p,'genes.tsv');bf<-file.path(p,'barcodes.tsv')}
  con<-if(grepl('gz$',mf))gzfile(mf,open='rt')else file(mf,open='rt');m<-as(readMM(con),'CsparseMatrix');close(con);rownames(m)<-make.unique(fread(gf,header=FALSE)[[2]]);bc<-fread(bf,header=FALSE)[[1]]
  z<-summarize(m,dd,paste(sid,bc,sep='__'));saveRDS(z,ck);res[[sid]]<-z;rm(m);gc()
 }
 # The historical wide matrix is read once, retaining only the common panel.
 ck<-file.path(out,'harmonized_alias_GSE143437.rds')
 if(file.exists(ck))res[['GSE143437']]<-readRDS(ck) else {
  con<-gzfile(file.path(historic,'GSE143437_DeMicheli_MuSCatlas_rawdata.txt.gz'),'rt');bc<-strsplit(readLines(con,n=1),'\t',fixed=TRUE)[[1]][-1];mat<-list();totals<-numeric(length(bc))
  repeat {lines<-readLines(con,n=250);if(!length(lines))break
   for(line in lines){parts<-strsplit(line,'\t',fixed=TRUE)[[1]];g<-parts[1];if(g=='Tmem8c')g<-'Mymk';if(g=='Gm7325')g<-'Mymx';v<-as.numeric(parts[-1]);totals<-totals+v;if(g%in%genes)mat[[g]]<-v}
  };close(con);m<-do.call(rbind,mat)
  z<-rbindlist(lapply(unique(md[dataset_id=='GSE143437',sample_id]),function(sid)summarize(m,md[sample_id==sid],bc,totals)))
  saveRDS(z,ck);res[['GSE143437']]<-z
 }
 # Retain the original cell QC/author mapping in the two additional NTX cohorts.
 for(acc in c('GSE159500','GSE162172')) {
  meta<-Sys.getenv(paste0('MUSCLE_',acc,'_METADATA'));stopifnot(file.exists(meta));ddall<-fread(meta)[passes_qc_and_singlet==TRUE]
  ddall[,`:=`(source_cell_id=cell_id,cell_type=coarse_cell_type,dataset_id=acc,age_group=ifelse(acc=='GSE162172','old','young'))]
  for(sid in unique(ddall$sample_id)){
   ck<-file.path(out,paste0('harmonized_alias_',acc,'_',sid,'.rds'));if(file.exists(ck)){res[[paste(acc,sid)]]<-readRDS(ck);next}
   dd<-ddall[sample_id==sid];p<-file.path(rawroot,'extracted',acc)
   if(acc=='GSE162172'){
    f<-list.files(p,pattern=paste0('_',sid,'_soupx.zip$'),full.names=TRUE);stopifnot(length(f)==1);members<-unzip(f,list=TRUE)$Name
    read_member<-function(nm)unzip(f,files=members[grepl(paste0('/',nm,'$'),members)],exdir=file.path(root,'data/age_harmonization'),junkpaths=TRUE)
    mf<-read_member('matrix.mtx');gf<-read_member('genes.tsv');bf<-read_member('barcodes.tsv');m<-as(readMM(mf),'CsparseMatrix');gs<-fread(gf,header=FALSE)[[1]];bc<-fread(bf,header=FALSE)[[1]];rownames(m)<-make.unique(gs);ids<-paste(sid,bc,sep='::');unlink(c(mf,gf,bf))
   }else{
    # Source metadata sample names are embedded in the author bundle prefixes.
    files<-list.files(p,pattern='_soupx_matrix.mtx.gz$',full.names=TRUE);found<-character()
    for(f in files){meta_f<-sub('_soupx_matrix.mtx.gz$','_meta.csv.gz',f);if(file.exists(meta_f)&&fread(meta_f,nrows=1)$sample[1]==sid)found<-f}
    stopifnot(length(found)==1);prefix<-sub('_matrix.mtx.gz$','',found);m<-as(readMM(gzfile(found)),'CsparseMatrix');rownames(m)<-make.unique(fread(paste0(prefix,'_genes.tsv.gz'),header=FALSE)[[2]]);ids<-fread(paste0(prefix,'_barcodes.tsv.gz'),header=FALSE)[[1]]
   }
   message('Common programmes: ',acc,' ',sid);z<-summarize(m,dd,ids);saveRDS(z,ck);res[[paste(acc,sid)]]<-z;rm(m);gc()
  }
 }
 z<-rbindlist(res);z[dataset_id=='GSE143437',age_class:='young'];fwrite(z,file.path(out,'harmonized_age_lineage_programs.tsv'),sep='\t');harmonized_age_RSI(z,root)
}

harmonized_age_RSI <- function(d,root) {
 suppressPackageStartupMessages({library(data.table);library(mgcv)})
 set.seed(20261125);out<-file.path(root,'results');d<-copy(d[n_cells>=20]);d[dataset_id=='GSE143437',age_class:='young'];stopifnot(all(d$age_class%in%c('young','old','geriatric')))
 pcols<-c('Activation','Proliferation','Commitment','Fusion','Remodeling','Resolution')
 mat<-scale(as.matrix(d[,..pcols]));mat[!is.finite(mat)]<-0;mat<-sweep(mat,1,apply(mat,1,max),'-')/.8;pp<-exp(mat)
 d[,`:=`(phase=as.numeric((pp/rowSums(pp))%*%c(.12,.28,.45,.64,.80,.94)),a=pmin(log1p(time_days)/log(8),1),w=1-exp(-n_cells/50))]
 d[,unit:=paste(dataset_id,sample_id,sep='::')];d[age_class=='geriatric',unit:=paste('geriatric',time_days,ifelse(grepl('_[AB]$',sample_id),'AB','CD'),sep='::')]
 calculate<-function(data,reference='pooled') {
  z<-copy(data);ref<-if(reference=='young')z[age_class=='young']else copy(z);z[,expected:=NA_real_]
  pred<-function(tr,te){if(nrow(tr)<6||uniqueN(tr$a)<4)return(rep(NA_real_,nrow(te)));k<-min(5,uniqueN(tr$a)-1);fit<-gam(phase~s(a,k=k),data=tr,weights=w,method='REML');as.numeric(predict(fit,te))}
  for(ct in unique(z$cell_type)){
   tr<-ref[cell_type==ct];ii<-which(z$cell_type==ct);z$expected[ii]<-pred(tr,z[ii])
   for(u in unique(tr$unit)){jj<-ii[z$unit[ii]==u];z$expected[jj]<-pred(tr[unit!=u],z[jj])}
  }
  z[,deviation:=phase-expected];sig<-max(.05,sd(z[unit%in%ref$unit,deviation],na.rm=TRUE))
  s<-z[is.finite(deviation),.(RSI=if(.N>=4)exp(-sum(w*(deviation-weighted.mean(deviation,w))^2)/sum(w)/sig^2)else NA_real_,n_lineages=.N),by=.(sample_id,unit,age_class,time_days)]
  units<-s[is.finite(RSI),.(RSI=weighted.mean(RSI,n_lineages),n_captures=.N,n_lineages=min(n_lineages)),by=.(unit,age_class,time_days)]
  list(samples=s,units=units,lineages=z)
 }
 tests<-list();all<-list()
 for(ref in c('young','pooled')) {
  r<-calculate(d,ref);all[[ref]]<-r$units
  fwrite(r$lineages,file.path(out,paste0('harmonized_age_',ref,'_lineages.tsv')),sep='\t')
 }
 # The pooled-age reference is independent of age labels. All groups are
 # cross-fitted symmetrically; label permutations do not require model refits.
 a<-copy(all$pooled);a[,age_class:=factor(age_class,c('young','old','geriatric'))]
 # Every included injury time must contain all three age groups.
 keep<-a[,.(n=uniqueN(age_class)),by=time_days][n==3,time_days];a<-a[time_days%in%keep]
 X0<-model.matrix(~factor(time_days),a);y<-a$RSI
 calc<-function(ages){X<-model.matrix(~factor(a$time_days)+factor(ages,levels=c('young','old','geriatric')));fit<-lm.fit(X,y);rss<-sum(fit$residuals^2);F<-((sum(lm.fit(X0,y)$residuals^2)-rss)/2)/(rss/(length(y)-ncol(X)));c(omnibus=F,old=unname(fit$coefficients[ncol(X)-1]),geriatric=unname(fit$coefficients[ncol(X)]))}
 obs<-calc(as.character(a$age_class));strata<-split(seq_len(nrow(a)),a$time_days)
 null<-replicate(10000,{aa<-as.character(a$age_class);for(ii in strata)aa[ii]<-sample(aa[ii]);calc(aa)})
 p<-c((1+sum(null[1,]>=obs[1]-1e-12))/10001,(1+rowSums(abs(null[2:3,,drop=FALSE])>=abs(obs[2:3])-1e-12))/10001)
 test<-data.table(contrast=c('omnibus_age','old_minus_young','geriatric_minus_young'),estimate=obs,permutation_p=p,n_units=nrow(a),n_times=length(keep));test[,holm_p:=p.adjust(permutation_p,'holm')]
 fwrite(data.table(statistic=rep(c("omnibus_age","old_minus_young","geriatric_minus_young"),each=ncol(null)),value=as.vector(t(null))),file.path(out,"harmonized_age_permutations.tsv.gz"),sep="\t")
 fwrite(rbindlist(all,idcol='reference'),file.path(out,'harmonized_age_RSI_units.tsv'),sep='\t');fwrite(test,file.path(out,'harmonized_age_tests.tsv'),sep='\t');print(test)
}

harmonized_lineage_age <- function(root) {
 suppressPackageStartupMessages(library(data.table));set.seed(20261125)
 out<-file.path(root,'results');d<-fread(file.path(out,'harmonized_age_pooled_lineages.tsv'))[is.finite(deviation)]
 d<-d[,.(deviation=weighted.mean(deviation,w),n_cells=sum(n_cells)),by=.(unit,age_class,time_days,cell_type)]
 d[,departure:=abs(deviation)];res<-list();coverage<-list()
 for(ct in unique(d$cell_type)) {
  a<-copy(d[cell_type==ct]);counts<-a[,.N,by=.(time_days,age_class)]
  times<-counts[N>=2,.(n_groups=uniqueN(age_class)),by=time_days][n_groups==3,time_days];a<-a[time_days%in%times]
  coverage[[ct]]<-data.table(cell_type=ct,n_times=length(times),n_units=nrow(a))
  if(length(times)<2||nrow(a)<12)next
  y<-a$departure
  fit_coeff<-function(ag){X<-model.matrix(~factor(a$time_days)+factor(ag,levels=c('young','old','geriatric')));tail(lm.fit(X,y)$coefficients,2)}
  obs<-fit_coeff(a$age_class);strata<-split(seq_len(nrow(a)),a$time_days)
  null<-replicate(10000,{ag<-a$age_class;for(ii in strata)ag[ii]<-sample(ag[ii]);fit_coeff(ag)})
  p<-(1+rowSums(abs(null)>=abs(obs)-1e-12))/10001
  res[[ct]]<-data.table(cell_type=ct,contrast=c('old_minus_young','geriatric_minus_young'),effect=unname(obs),permutation_p=p,n_units=nrow(a),n_times=length(times))
 }
 z<-rbindlist(res);z[,BH_q:=p.adjust(permutation_p,'BH')]
 fwrite(z,file.path(out,'harmonized_lineage_age_tests.tsv'),sep='\t');fwrite(rbindlist(coverage),file.path(out,'harmonized_lineage_age_coverage.tsv'),sep='\t');print(z[order(BH_q)])
}

args<-commandArgs(TRUE)
root<-if(length(args))args[1]else"."
mode<-if(length(args)>1)args[2]else"lifespan"
switch(mode,lifespan=lifespan_validation(),independent=independent_validation(root),denervation=denervation_validation(root),depletion=depletion_gene_validation(root),harmonized=harmonized_age_validation(root),lineage_age=harmonized_lineage_age(root),harmonized_scores=harmonized_age_RSI(data.table::fread(file.path(root,"results/harmonized_age_lineage_programs.tsv")),root),stop("Unknown biological validation: ",mode))
