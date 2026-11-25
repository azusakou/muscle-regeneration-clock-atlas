suppressPackageStartupMessages({library(data.table); library(mgcv)})
set.seed(20261125)
args <- commandArgs(TRUE)
root <- if (length(args)) args[1] else "."
input <- file.path(root, "input")
output <- file.path(root, "results")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
write_table <- function(x, name) fwrite(x, file.path(output, paste0(name,".tsv")), sep="\t")
wm <- function(x,w) {o<-order(x);x<-x[o];w<-w[o];x[which(cumsum(w)>=sum(w)/2)[1]]}
wv <- function(x,w) sum(w*(x-weighted.mean(x,w))^2)/sum(w)
scale_safe <- function(x) {s<-sd(x);if (!is.finite(s)||s<1e-8) return(x*0);(x-mean(x))/s}
soft <- function(z, centers=c(.12,.28,.45,.64,.80,.94), temperature=.8) {
 z<-sweep(z,1,apply(z,1,max),"-")/temperature;p<-exp(z);as.numeric((p/rowSums(p))%*%centers)
}
x <- fread(file.path(input,"lineage_scores.tsv.gz"))
scores <- paste0("z_state_score_",c("Activation","Proliferation","Commitment","Fusion","Remodeling","Resolution"))
stopifnot(all(scores %in% names(x)))
x[, unit:=paste(dataset_id,biological_group,sep="::")]

# Anchors are fitted separately for each lineage; cross-fitting removes the
# evaluated young biological unit from its own reference trajectory.
calculate <- function(d, anchor="injury", weights="recovery", floor=.05, min_lineages=4L,
                      reference="all_young", crossfit=TRUE, supplied_phase=NULL) {
 d<-copy(d)
 d[,phase:=if(is.null(supplied_phase)) phase_harmonized else supplied_phase]
 d[,w:=switch(weights,recovery=1-exp(-n_cells/50),equal=rep(1,.N),sqrt=sqrt(n_cells))]
 d[,consensus:=wm(phase,w),by=sample_uid]
 d[,a:=switch(anchor,injury=target_tau,consensus=consensus,lolo=vapply(seq_len(.N),function(i)wm(phase[-i],w[-i]),0.0)),by=sample_uid]
 d[,expected:=NA_real_]
 ref<-d[age_class=="young" & injury_model_std %in% c("NTX","CTX","BaCl2")]
 if(reference=="NTX") ref<-ref[injury_model_std=="NTX"]
 for(ct in unique(d$cell_type)) {
  rr<-ref[cell_type==ct];if(nrow(rr)<6||uniqueN(rr$a)<4) next
  predict_ref<-function(tr,te) {
   if(nrow(tr)<6||uniqueN(tr$a)<4) return(rep(NA_real_,nrow(te)))
   k<-min(5,uniqueN(tr$a)-1)
   fit<-gam(phase~s(a,k=k),data=tr,weights=w,method="REML")
   as.numeric(predict(fit,newdata=te))
  }
  ii<-which(d$cell_type==ct);d$expected[ii]<-predict_ref(rr,d[ii])
  if(crossfit) for(id in unique(rr$unit)) {
   jj<-ii[d$unit[ii]==id];d$expected[jj]<-predict_ref(rr[unit!=id],d[jj])
  }
 }
 d[,deviation:=phase-expected]
 refdev<-d[age_class=="young" & injury_model_std %in% c("NTX","CTX","BaCl2"),deviation]
 sigma<-max(floor,sd(refdev,na.rm=TRUE))
 r<-d[is.finite(deviation),.(dispersion=wv(deviation,w),n_lineages=.N,n_cells=sum(n_cells)),
        by=.(dataset_id,sample_uid,unit,time_days,age_class,injury_model_std)]
 r[,RSI:=ifelse(n_lineages>=min_lineages,exp(-dispersion/sigma^2),NA_real_)]
 r[,reference_sigma:=sigma]
 list(samples=r,lineages=d)
}
base <- calculate(x,crossfit=FALSE)
write_table(base$samples,"injury_anchored_RSI")
old<-fread(file.path(input,"age_RSI_original.tsv"))
check<-merge(old,base$samples,by="sample_uid")
stopifnot(max(abs(check$RSI.x-check$RSI.y),na.rm=TRUE)<1e-9)

variants<-list(injury_crossfit=list(),consensus_crossfit=list(anchor="consensus"),
 lolo_crossfit=list(anchor="lolo"),equal_weights=list(weights="equal"),
 sqrt_weights=list(weights="sqrt"),NTX_reference=list(reference="NTX"),
 floor_003=list(floor=.03),floor_008=list(floor=.08),minimum_3=list(min_lineages=3L),
 minimum_5=list(min_lineages=5L))
for(t in c(.4,1.2,1.6)) variants[[paste0("temperature_",t)]]<-list(supplied_phase=soft(as.matrix(x[,..scores]),temperature=t))
variants$swap_remodel_resolution<-list(supplied_phase=soft(as.matrix(x[,..scores]),centers=c(.12,.28,.45,.64,.94,.80)))
variants$uniform_centers<-list(supplied_phase=soft(as.matrix(x[,..scores]),centers=seq(.1,.9,length.out=6)))
for(j in seq_along(scores)) variants[[paste0("omit_",sub("z_state_score_","",scores[j]))]]<-
 list(supplied_phase=soft(as.matrix(x[,scores[-j],with=FALSE]),centers=c(.12,.28,.45,.64,.80,.94)[-j]))
pc<-prcomp(x[,..scores],center=TRUE,scale.=FALSE)$x[,1]
if(cor(pc,x$target_tau)<0)pc<- -pc
variants$PC1<-list(supplied_phase=(pc-min(pc))/diff(range(pc)))
all<-list(original=base$samples);summ<-list();lineage<-list()
for(nm in names(variants)) {
 z<-do.call(calculate,c(list(d=x),variants[[nm]]));all[[nm]]<-z$samples
 a<-merge(base$samples[,.(sample_uid,baseline=RSI)],z$samples[,.(sample_uid,RSI)],by="sample_uid")
 summ[[nm]]<-data.table(variant=nm,n_complete=sum(complete.cases(a)),
  spearman=cor(a$baseline,a$RSI,method="spearman",use="complete.obs"),
  median_absolute_difference=median(abs(a$baseline-a$RSI),na.rm=TRUE))
 if(nm %in% c("injury_crossfit","consensus_crossfit","lolo_crossfit")) lineage[[nm]]<-z$lineages[,.(sample_uid,unit,time_days,age_class,injury_model_std,cell_type,n_cells,phase,a,expected,deviation)]
}
write_table(rbindlist(summ),"RSI_sensitivity")
write_table(rbindlist(all,idcol="variant"),"RSI_sensitivity_samples")
write_table(rbindlist(lineage,idcol="variant"),"phase_anchor_sensitivity")

aggregate_units <- function(d) d[is.finite(RSI),.(RSI=weighted.mean(RSI,n_lineages),
 dispersion=weighted.mean(dispersion,n_lineages),n_captures=.N,n_lineages=min(n_lineages)),
 by=.(unit,age_class,time_days)]
age <- aggregate_units(all$injury_crossfit[injury_model_std=="NTX"])
write_table(age,"age_RSI_biological_units")
exact<-function(a,b) {
 v<-c(a,b);obs<-mean(b)-mean(a);n<-length(a)
 perms<-combn(seq_along(v),n)
 e<-apply(perms,2,function(ii)mean(v[-ii])-mean(v[ii]))
 c(effect=obs,p=mean(abs(e)>=abs(obs)-1e-12),allocations=ncol(perms))
}
contrasts<-list()
for(tt in sort(unique(age$time_days))) {
 a<-age[time_days==tt & age_class=="young",RSI];b<-age[time_days==tt & age_class=="geriatric",RSI]
 if(!length(a)||!length(b))next
 e<-exact(a,b)
 boots<-replicate(10000,mean(sample(b,length(b),TRUE))-mean(sample(a,length(a),TRUE)))
 ci<-if(min(length(a),length(b))>=2)quantile(boots,c(.025,.975)) else c(NA,NA)
 contrasts[[as.character(tt)]]<-data.table(time_days=tt,n_young=length(a),n_geriatric=length(b),
  effect=e[1],lower=ci[1],upper=ci[2],exact_p=e[2],n_allocations=e[3],
  young_sd=sd(a),geriatric_sd=sd(b))
}
contrasts<-rbindlist(contrasts);contrasts[,holm_p:=p.adjust(exact_p,"holm")]
write_table(contrasts,"age_RSI_inference")
vage<-rbindlist(lapply(names(all),function(nm) {
 a<-aggregate_units(all[[nm]][injury_model_std=="NTX"])
 b<-a[,.(mean_RSI=mean(RSI),n=.N),by=.(time_days,age_class)]
 w<-dcast(b,time_days~age_class,value.var="mean_RSI")
 if(!all(c("young","geriatric") %in% names(w)))return(NULL)
 w[,.(variant=nm,time_days,geriatric_minus_young=geriatric-young)]
}));write_table(vage,"age_contrasts_sensitivity")

# Sample-level associations describe technical dependence; they are not causal
# corrections and include only eligible young/old/geriatric NTX samples.
ntx<-base$samples[injury_model_std=="NTX" & is.finite(RSI)]
tech<-rbindlist(lapply(c("n_cells","n_lineages"),function(v) {
 ct<-cor.test(ntx$RSI,ntx[[v]],method="spearman",exact=FALSE)
 data.table(variable=v,n=nrow(ntx),rho=unname(ct$estimate),p=ct$p.value)
}));write_table(tech,"RSI_recovery_associations")

# Age-associated lineage departures at the two within-accession matched times.
d<-lineage$injury_crossfit[injury_model_std=="NTX" & time_days %in% c(1,3.5)]
d<-d[,.(deviation=weighted.mean(deviation,1-exp(-n_cells/50))),by=.(unit,age_class,time_days,cell_type)]
lr<-list()
for(ct in unique(d$cell_type))for(tt in c(1,3.5)) {
 a<-d[cell_type==ct & time_days==tt & age_class=="young" & is.finite(deviation),deviation]
 b<-d[cell_type==ct & time_days==tt & age_class=="geriatric" & is.finite(deviation),deviation]
 if(!length(a)||!length(b))next
 e<-exact(a,b);lr[[paste(ct,tt)]]<-data.table(cell_type=ct,time_days=tt,n_young=length(a),n_geriatric=length(b),difference=e[1],exact_p=e[2])
}
lr<-rbindlist(lr);lr[,BH_q:=p.adjust(exact_p,"BH")];write_table(lr,"lineage_age_comparisons")

# Uncertainty at the independent-study level for the original grouped benchmark.
m<-fread(file.path(input,"validation_metrics.tsv"))[validation=="leave_study_out"]
w<-dcast(m,held_out~model,value.var="MAE_tau")
diffs<-w$coupled_regeneration_clock-w$independent_celltype_clocks
sg<-as.matrix(expand.grid(rep(list(c(-1,1)),length(diffs))))
write_table(data.table(comparison="coupled_minus_independent",mean_study_difference=mean(diffs),
 lower=mean(diffs)-qt(.975,4)*sd(diffs)/sqrt(5),upper=mean(diffs)+qt(.975,4)*sd(diffs)/sqrt(5),
 exact_signflip_p=mean(abs(sg%*%diffs/length(diffs))>=abs(mean(diffs))-1e-12),n_studies=5),"model_comparison_uncertainty")
print(contrasts);print(rbindlist(summ));print(lr[order(-abs(difference))]);print(tech)

# Hold the evaluated lineage panel fixed, excluding composition-dependent
# weights. This evaluates recovery sensitivity rather than tissue abundance.
panel<-c("FAP","macrophage","endothelial","T_NK")
fx<-x[cell_type%in%panel];ids<-fx[,.N,by=sample_uid][N==4,sample_uid]
fixed<-calculate(fx[sample_uid%in%ids],weights="equal")
write_table(fixed$samples,"fixed_panel_RSI")

# Equal cell recovery within sample-lineage strata. Each replicate starts
# from existing per-cell scores, recomputes medians and study scaling, and
# refits references. The evaluated strata are identical in each replicate.
if(file.exists(file.path(input,"core_cell_metadata.rds"))) {
 md<-as.data.table(readRDS(file.path(input,"core_cell_metadata.rds")))
 map<-c(Myogenic_progenitor="MuSC_myogenic",Myonucleus="myonucleus",FAP="FAP",
 Macrophage="macrophage",Neutrophil="neutrophil",Endothelial="endothelial",
 Pericyte="pericyte",T_cell="T_NK",B_cell="other_immune",Schwann="other_major")
 md<-md[passes_qc_and_singlet==TRUE];md[,cell_type:=unname(map[unified_cell_type])]
 md<-md[!is.na(cell_type)];md[,sample_uid:=paste(dataset_id,sample_id,sep="::")]
 rawcols<-sub("^z_","",scores)
 groups<-md[,.(indices=list(.I),n=.N),by=.(sample_uid,cell_type,dataset_id)][n>=50]
 groups[,key:=paste(sample_uid,cell_type,sep="|")]
 keys<-paste(x$sample_uid,x$cell_type,sep="|")
 keep<-keys%in%groups$key
 xx<-x[keep]
 rebuild<-function(rarefy=FALSE) {
  a<-rbindlist(lapply(seq_len(nrow(groups)),function(i){
   ii<-groups$indices[[i]];if(rarefy)ii<-sample(ii,50)
   z<-as.data.table(lapply(md[ii,..rawcols],median,na.rm=TRUE))
   z[,`:=`(sample_uid=groups$sample_uid[i],cell_type=groups$cell_type[i],dataset_id=groups$dataset_id[i])];z
  }))
  for(j in seq_along(scores)) a[,(scores[j]):=scale_safe(get(rawcols[j])),by=dataset_id]
  a[,phase:=soft(as.matrix(.SD)),.SDcols=scores]
  a<-a[match(paste(xx$sample_uid,xx$cell_type,sep="|"),paste(sample_uid,cell_type,sep="|"))]
  calculate(xx,weights="equal",supplied_phase=a$phase)$samples
 }
 full<-rebuild();depth<-rbindlist(lapply(seq_len(100),function(i){
  a<-rebuild(TRUE);a[,iteration:=i];a
 }))
 depth<-merge(depth,full[,.(sample_uid,full_RSI=RSI)],by="sample_uid")
 write_table(depth,"equal_depth_RSI_samples")
 write_table(depth[is.finite(RSI)&is.finite(full_RSI),.(n_samples=.N,
   spearman=cor(RSI,full_RSI,method="spearman"),median_absolute_difference=median(abs(RSI-full_RSI))),by=iteration],"equal_depth_RSI_summary")
}

if(file.exists(file.path(input,"depletion_phase.tsv"))) {
 d<-fread(file.path(input,"depletion_phase.tsv"));d<-d[n_cells>=20 & is.finite(deviation)]
 sig<-max(.05,sqrt(mean(d[condition=="sham",deviation]^2)))
 # Reproduce the original scale from all sham strata, including low recovery.
 orig<-fread(file.path(input,"depletion_phase.tsv"))
 sig<-max(.05,sqrt(mean(orig[condition=="sham",deviation]^2,na.rm=TRUE)))
 panels<-list(common_four=c("FAP","macrophage","endothelial","myonucleus"),
  without_macrophage=c("FAP","endothelial","myonucleus"))
 res<-rbindlist(lapply(names(panels),function(p) {
  dd<-d[coarse_cell_type%in%panels[[p]],.(RSI=exp(-wv(deviation,rep(1,.N))/sig^2),n=.N),by=.(sample_id,condition)]
  dd<-dd[n==length(panels[[p]])];e<-exact(dd[condition=="sham",RSI],dd[condition=="clodronate",RSI])
  data.table(panel=p,n_animals=nrow(dd),n_lineages=length(panels[[p]]),effect=e[1],exact_p=e[2])
 }))
 write_table(res,"depletion_fixed_panel")
}
