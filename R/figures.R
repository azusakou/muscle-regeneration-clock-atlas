suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(TRUE);root<-if(length(args))args[1]else"."
res<-file.path(root,"results");fig<-Sys.getenv("MUSCLE_FIGURE_DIR",file.path(root,"figures"));dir.create(fig,recursive=TRUE,showWarnings=FALSE)
theme_set(theme_classic(base_size=8,base_family="Helvetica")+
 theme(plot.title=element_text(size=9,face="bold"),plot.tag=element_text(face="bold",size=11),
  strip.background=element_blank(),strip.text=element_text(face="bold",size=8),
  legend.title=element_blank(),axis.text=element_text(colour="black",size=7)))
savefig<-function(p,name,w=7.2,h=5)ggsave(file.path(fig,paste0(name,".pdf")),p,width=w,height=h,device=pdf,useDingbats=FALSE)
ages<-c(young="#2878B5",old="#8599B1",geriatric="#D1495B")
a<-fread(file.path(res,"age_RSI_biological_units.tsv"));a[,age_class:=factor(age_class,names(ages))]
means<-a[,.(RSI=mean(RSI),n=.N),by=.(age_class,time_days)]
p1<-ggplot(a,aes(time_days,RSI,colour=age_class))+geom_point(size=1.5,position=position_jitter(width=.10,height=0,seed=11),alpha=.7)+
 geom_line(data=means,aes(group=age_class),linewidth=.5)+geom_point(data=means,size=2.3)+
 scale_colour_manual(values=ages)+scale_x_continuous(breaks=c(0,1,2,3.5,5,7))+ylim(0,1)+
 labs(x="Days after injury",y="Cross-fitted RSI",title="Biological groups",tag="a")+theme(legend.position="top")
c<-fread(file.path(res,"age_RSI_inference.tsv"));c[,day:=factor(time_days,rev(time_days))]
p2<-ggplot(c,aes(effect,day))+geom_vline(xintercept=0,colour="grey65",linewidth=.35)+
 geom_errorbar(aes(xmin=lower,xmax=upper),width=.16,orientation="y",na.rm=TRUE,linewidth=.5)+geom_point(size=2,colour=ages[3])+
 geom_text(aes(label=paste0("n=",n_young,"/",n_geriatric)),x=.64,hjust=1,size=2.2)+
 coord_cartesian(xlim=c(-.45,.68),clip="off")+labs(x="Geriatric - young RSI (95% interval)",y="Days after injury",title="Observed contrasts",tag="b")
l<-fread(file.path(res,"phase_anchor_sensitivity.tsv"))[variant=="injury_crossfit" & injury_model_std=="NTX"]
l<-l[,.(deviation=weighted.mean(deviation,1-exp(-n_cells/50),na.rm=TRUE)),by=.(unit,age_class,time_days,cell_type)]
l<-l[,.(deviation=mean(deviation,na.rm=TRUE)),by=.(age_class,time_days,cell_type)]
l[,cell_type:=factor(cell_type,rev(c("MuSC_myogenic","myonucleus","macrophage","FAP","endothelial","pericyte","T_NK","neutrophil","other_immune","other_major")))]
p3<-ggplot(l,aes(factor(time_days),cell_type,fill=deviation))+geom_tile(colour="white",linewidth=.25)+
 facet_wrap(~factor(age_class,levels=names(ages)),nrow=1)+
 scale_fill_gradient2(low="#2878B5",mid="white",high="#D1495B",midpoint=0,limits=c(-.3,.3),oob=scales::squish,na.value="white")+
 labs(x="Days after injury",y=NULL,fill="Deviation",tag="c",title="Lineage deviation from cross-fitted reference")+
 scale_y_discrete(labels=function(x)gsub("_"," ",x))+theme(legend.position="right",axis.text.y=element_text(size=6.5))
savefig((p1+p2)/p3+plot_layout(heights=c(1,1.1)),"Figure3_age_RSI",h=5.5)

# Uncertainty is defined across five held-out studies, not across cells.
m<-fread(file.path(root,"input/validation_metrics.tsv"))[validation=="leave_study_out"]
w<-dcast(m,held_out~model,value.var="MAE_tau")
w[,difference:=coupled_regeneration_clock-independent_celltype_clocks]
u<-fread(file.path(res,"model_comparison_uncertainty.tsv"))
q1<-ggplot(w,aes(difference,reorder(held_out,difference)))+geom_vline(xintercept=0,colour="grey65")+
 geom_segment(aes(x=0,xend=difference,yend=held_out),linewidth=.6,colour="#8599B1")+geom_point(size=2)+
 labs(x="Coupled - independent MAE",y=NULL,title="Paired study differences",tag="a")
q2<-ggplot(u,aes(mean_study_difference,comparison))+geom_vline(xintercept=0,colour="grey65")+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.1)+geom_point(size=2.5,colour="#2878B5")+
 labs(x="Study-mean difference (95% t interval)",y=NULL,title="Five-study uncertainty",tag="b")+scale_y_discrete(labels="Coupled -\nindependent")
savefig(q1+q2+plot_layout(widths=c(1.1,1)),"SupplementaryFigure1",h=2.5)

sens<-fread(file.path(res,"RSI_sensitivity.tsv"))
labels<-c(injury_crossfit="Injury anchor, cross-fitted",consensus_crossfit="Consensus anchor",lolo_crossfit="Leave one lineage out",
 equal_weights="Equal weights",sqrt_weights="Square-root count weights",NTX_reference="Notexin-only reference",
 floor_003="Dispersion floor 0.03",floor_008="Dispersion floor 0.08",minimum_3="Minimum 3 lineages",minimum_5="Minimum 5 lineages",
 temperature_0.4="Temperature 0.4",temperature_1.2="Temperature 1.2",temperature_1.6="Temperature 1.6",
 swap_remodel_resolution="Swap remodelling and resolution",uniform_centers="Equally spaced centres",PC1="First principal component")
sens[,label:=ifelse(variant%in%names(labels),labels[variant],paste("Omit",sub("omit_","",variant))) ]
r1<-ggplot(sens,aes(spearman,reorder(label,spearman)))+geom_segment(aes(x=0,xend=spearman,yend=label),colour="#8599B1")+
 geom_point(size=1.5,colour="#2878B5")+xlim(0,1)+labs(x="Rank correlation with original injury anchor",y=NULL,tag="a")+theme(axis.text.y=element_text(size=6))
ac<-fread(file.path(res,"age_contrasts_sensitivity.tsv"))[variant%in%c("injury_crossfit","consensus_crossfit","lolo_crossfit")]
ac[,variant:=factor(variant,c("injury_crossfit","consensus_crossfit","lolo_crossfit"),c("Injury anchor","Consensus","Leave one lineage out"))]
r2<-ggplot(ac,aes(time_days,geriatric_minus_young,colour=variant))+geom_hline(yintercept=0,colour="grey70")+
 geom_line(linewidth=.5)+geom_point(size=1.4)+scale_colour_manual(values=c("#2878B5","#8599B1","#D1495B"))+
 labs(x="Days after injury",y="Geriatric - young RSI",tag="b")+theme(legend.position="top",legend.text=element_text(size=6))
dep<-fread(file.path(res,"equal_depth_RSI_summary.tsv"))
r3<-ggplot(dep,aes(x="50 cells per lineage",y=spearman))+geom_boxplot(width=.35,outlier.shape=NA)+
 geom_point(position=position_jitter(width=.12,height=0,seed=1),size=.55,alpha=.55,colour="#2878B5")+
 coord_cartesian(ylim=c(0,1))+labs(x=NULL,y="Rank correlation with full depth",tag="c")
savefig(r1/(r2+r3)+plot_layout(heights=c(1.8,1)),"SupplementaryFigure5",h=7)

b<-fread(file.path(res,"biological_module_scores.tsv"))
mac<-b[cell_type=="Macrophage" & n_cells>=20 & module%in%c("inflammatory_exercise","restorative")]
mac[,module:=factor(module,c("inflammatory_exercise","restorative"),c("Inflammatory (exercise markers)","Restorative"))]
mm<-mac[,.(score=mean(score)),by=.(time_days,age_class,module)]
s1<-ggplot(mac,aes(time_days,score,colour=age_class))+geom_point(size=1,alpha=.7,position=position_jitter(width=.08,height=0,seed=12))+
 geom_line(data=mm,aes(group=age_class),linewidth=.5)+facet_wrap(~module,nrow=1)+
 scale_colour_manual(values=ages)+labs(x="Days after injury",y="Mean log-normalized expression",tag="a")+theme(legend.position="top")
ctrl<-b[module%in%c("Schwann_identity","NMJ_postsynaptic") & n_cells>=20]
ctrl[,module:=factor(module,c("NMJ_postsynaptic","Schwann_identity"),c("Postsynaptic NMJ","Schwann identity"))]
ctrl[,cell_type:=factor(cell_type,rev(c("Schwann","Myogenic_progenitor","Myonucleus","FAP","Macrophage","Endothelial","Pericyte","T_cell","B_cell","Neutrophil")))]
cm<-ctrl[,.(score=mean(score)),by=.(module,cell_type)]
s2<-ggplot(ctrl,aes(score,cell_type))+geom_point(size=.6,alpha=.35,position=position_jitter(width=0,height=.13,seed=2),colour="#8599B1")+
 geom_point(data=cm,size=1.8,colour="#2878B5")+facet_wrap(~module,scales="free_x",nrow=1)+
 labs(x="Module expression (biological groups)",y=NULL,tag="b")+scale_y_discrete(labels=function(x)gsub("_"," ",x))
assoc<-fread(file.path(res,"macrophage_FAP_associations.tsv"))
aa<-melt(assoc,id.vars="gene",measure.vars=c("raw_rho","adjusted_rho"),variable.name="estimate",value.name="rho")
s3<-ggplot(aa,aes(rho,gene,colour=estimate))+geom_vline(xintercept=0,colour="grey65")+geom_point(position=position_dodge(width=.35),size=2)+
 scale_colour_manual(values=c("#8599B1","#2878B5"),labels=c("Unadjusted","Time/age adjusted"))+
 labs(x="Macrophage ligand - FAP matrix rank correlation",y=NULL,tag="c")+theme(legend.position="top")
savefig(s1/s2/s3+plot_layout(heights=c(1.1,1.5,.8)),"SupplementaryFigure6",h=7)

human<-fread(file.path(root,"input/human_age_phase.tsv"))
hp<-ggplot(human,aes(coarse_cell_type,reorder(donor,as.integer(gsub("HM","",donor))),fill=log10(n_cells+1)))+
 geom_tile(colour="white",linewidth=.4)+geom_text(aes(label=n_cells),size=2.3)+facet_grid(age_class~.,scales="free_y",space="free_y")+
 scale_fill_gradient(low="white",high="#2878B5",name="log10(n + 1)")+
 labs(x=NULL,y="Donor",title="Recovered cells or nuclei by donor and lineage")+
 theme(axis.text.x=element_text(angle=35,hjust=1),legend.title=element_text(size=7))
savefig(hp,"SupplementaryFigure3",h=4.8)

# Schematic: programme positions are illustrative, not inferred timings.
pdf(file.path(fig,"Figure8_model.pdf"),width=7.2,height=6.5,family="Helvetica",useDingbats=FALSE)
grid.newpage()
txt<-function(x,y,label,size=8,font=1,just="left",col="#20252C")grid.text(label,x,y,just=just,gp=gpar(fontsize=size,fontface=font,col=col))
box<-function(x,y,w,h,label,fill="#EDF2F6",size=7.8,col="#617A8D"){
 grid.roundrect(x,y,w,h,r=unit(.018,"snpc"),gp=gpar(fill=fill,col=col,lwd=.6));txt(x,y,label,size,just="centre")
}
link<-function(x1,y1,x2,y2,dashed=FALSE,arrow=FALSE)grid.lines(c(x1,x2),c(y1,y2),arrow=if(arrow)arrow(length=unit(1.5,"mm"),type="closed")else NULL,gp=gpar(col="#617A8D",lwd=.65,lty=if(dashed)2 else 1))
txt(.025,.975,"a",11,2);txt(.068,.975,"Overlapping programmes during normal repair",9,2)
txt(.25,.935,"Early response",7.7,2);txt(.64,.935,"Repair and remodelling",7.7,2)
link(.26,.915,.94,.915,arrow=TRUE)
lanes<-c("Myogenic","Macrophage","FAP","Endothelial","Other immune")
states<-c("Activation / proliferation   >   commitment / fusion","Inflammatory recruitment   >   restorative states","Transient expansion   >   matrix turnover / clearance","Vascular response   >   niche reconstruction","Innate recruitment   >   immune regulation")
for(i in seq_along(lanes)){
 y<-.878-(i-1)*.056;txt(.04,y,lanes[i],8,2)
 box(.61,y,.68,.041,states[i],fill=if(i%in%c(1,2))"#EAF1F6"else"#F3F5F6",size=7.1)
}
txt(.27,.588,"Positions are schematic; programmes overlap across lineages.",6.8,col="#5A6771")
txt(.025,.545,"b",11,2);txt(.068,.545,"Observed state relationships",9,2)
box(.265,.449,.43,.14,"Age comparison\nCross-fitted RSI lower at 3 of 6 times\nNo significant time-specific age contrast\nMixed timing and programme-strength effects",fill="#EAF1F6",size=7.6)
box(.745,.449,.43,.14,"Macrophage depletion\nFive eligible animals; exact P = 0.20\nFour-lineage contrast: -0.354\nWithout macrophages: +0.072 (three lineages)",fill="#EAF1F6",size=7.6)
txt(.5,.352,"State alignment is reference dependent and is not a tissue-recovery score.",7.7,2,just="centre")
txt(.025,.31,"c",11,2);txt(.068,.31,"Mechanistic context and testable links",9,2)
cx<-c(.145,.38,.615,.85)
top<-c("mTORC1","Notch","Wnt / WISP1","Inflammatory signalling")
bot<-c("MuSC alerting / translation\nAutophagy balance","MuSC quiescence\nVascular niche","Myogenic / stromal\nprogramme support","Macrophage transitions\nTNF / TGF-beta - FAP")
for(i in seq_along(cx)){
 box(cx[i],.261,.215,.04,top[i],fill="#F3F5F6",size=7.7)
 link(cx[i],.241,cx[i],.216,arrow=TRUE)
 box(cx[i],.18,.215,.063,bot[i],fill="white",size=7)
 link(cx[i],.148,cx[i],.111,dashed=TRUE,arrow=TRUE)
}
box(.5,.081,.92,.05,"Expression-derived phase and programme strength",fill="#EAF1F6",size=8)
link(.045,.025,.095,.025);txt(.11,.025,"Supported in cited experimental literature",6.7)
link(.58,.025,.63,.025,dashed=TRUE);txt(.645,.025,"Proposed link to present observations",6.7)
dev.off()
