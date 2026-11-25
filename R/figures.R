suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(TRUE);root<-if(length(args))args[1]else"."
res<-file.path(root,"results");fig<-Sys.getenv("MUSCLE_FIGURE_DIR",file.path(root,"figures"));dir.create(fig,recursive=TRUE,showWarnings=FALSE)
theme_set(theme_classic(base_size=8,base_family="Helvetica")+
 theme(plot.title=element_text(size=9,face="bold"),plot.tag=element_text(face="bold",size=11),
  strip.background=element_blank(),strip.text=element_text(face="bold",size=8),
  legend.title=element_blank(),axis.text=element_text(colour="black",size=7)))
pdf_device<-function(filename,width,height,...) {
 if(capabilities("aqua")) quartz(type="pdf",file=filename,width=width,height=height,family="Helvetica")
 else cairo_pdf(filename=filename,width=width,height=height,family="Helvetica")
}
savefig<-function(p,name,w=7.2,h=5)ggsave(file.path(fig,paste0(name,".pdf")),p,width=w,height=h,device=pdf_device)
ages<-c(young="#376B91",old="#C49A48",geriatric="#BA6048")
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
 geom_point(size=1.5,colour="#2878B5")+xlim(0,1)+labs(x="Rank correlation with injury anchor before cross-fitting",y=NULL,tag="a")+theme(axis.text.y=element_text(size=6))
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
pdf_device(file.path(fig,"Figure8_model.pdf"),width=7.2,height=6.5)
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


# Frozen source measurements preserve the figure inputs without large matrices.
local({
out<-fig
sdout<-function(x,name)invisible(NULL)
source_tables<-list(
"audited_cell_counts.tsv"="dataset\tfinal_cells\nreused_core\t220854\nGSE159500\t6039\nGSE162172\t76133\nGSE205690\t9883\nGSE291782\t43113\nGSE214544\t44242\nGSE167186\t116369\n",
"complete_dynamics/all_validation_metrics.tsv"="validation\theld_out\tmodel\tn_samples\tMAE_tau\tRMSE_tau\nleave_study_out\tGSE138826\tabsolute_time_model\t7\t0.326695409924687\t0.366985949307812\nleave_study_out\tGSE138826\tcell_state_model\t7\t0.281528669852217\t0.3457700600648\nleave_study_out\tGSE138826\tpseudotime_model\t7\t0.391482215759864\t0.438891471867057\nleave_study_out\tGSE138826\tindependent_celltype_clocks\t7\t0.269886406509354\t0.306272080913165\nleave_study_out\tGSE138826\tcoupled_regeneration_clock\t7\t0.239203517211606\t0.297952860700845\nleave_study_out\tGSE143437\tabsolute_time_model\t10\t0.352832083357372\t0.401304406223996\nleave_study_out\tGSE143437\tcell_state_model\t10\t0.286370716394518\t0.356629617586606\nleave_study_out\tGSE143437\tpseudotime_model\t10\t0.350128540619084\t0.401674462081914\nleave_study_out\tGSE143437\tindependent_celltype_clocks\t10\t0.293614865168456\t0.365002456772885\nleave_study_out\tGSE143437\tcoupled_regeneration_clock\t10\t0.261358178563817\t0.332454544939944\nleave_study_out\tGSE159500\tabsolute_time_model\t2\t0.235839583213141\t0.304368679833822\nleave_study_out\tGSE159500\tcell_state_model\t2\t0.15893110709777\t0.17178119795901\nleave_study_out\tGSE159500\tpseudotime_model\t2\t0.226274561689089\t0.295741259427446\nleave_study_out\tGSE159500\tindependent_celltype_clocks\t2\t0.183592537743201\t0.22891310443415\nleave_study_out\tGSE159500\tcoupled_regeneration_clock\t2\t0.216742063766404\t0.220719227568831\nleave_study_out\tGSE205690\tabsolute_time_model\t6\t0.397439387029863\t0.42829985241144\nleave_study_out\tGSE205690\tcell_state_model\t6\t0.389441694802323\t0.460953342459878\nleave_study_out\tGSE205690\tpseudotime_model\t6\t0.399803911908376\t0.430484208616222\nleave_study_out\tGSE205690\tindependent_celltype_clocks\t6\t0.407131136818283\t0.430869672267542\nleave_study_out\tGSE205690\tcoupled_regeneration_clock\t6\t0.37220339734079\t0.40219051975007\nleave_study_out\tGSE232106\tabsolute_time_model\t8\t0.194987500240385\t0.20810829598903\nleave_study_out\tGSE232106\tcell_state_model\t8\t0.0865326354809992\t0.0950860496899139\nleave_study_out\tGSE232106\tpseudotime_model\t8\t0.19498982534619\t0.207631758997896\nleave_study_out\tGSE232106\tindependent_celltype_clocks\t8\t0.157897132433368\t0.185021923669608\nleave_study_out\tGSE232106\tcoupled_regeneration_clock\t8\t0.323395132219366\t0.383595165795417\nleave_time_point_out\t0.5\tabsolute_time_model\t1\t0.400566470621652\t0.400566470621652\nleave_time_point_out\t0.5\tcell_state_model\t1\t0.263036147138886\t0.263036147138886\nleave_time_point_out\t0.5\tpseudotime_model\t1\t0.538027991433227\t0.538027991433227\nleave_time_point_out\t0.5\tindependent_celltype_clocks\t1\t0.382012280145079\t0.382012280145079\nleave_time_point_out\t0.5\tcoupled_regeneration_clock\t1\t0.298386682771241\t0.298386682771241\nleave_time_point_out\t10\tabsolute_time_model\t1\t0.429602669755451\t0.429602669755451\nleave_time_point_out\t10\tcell_state_model\t1\t0.221046099400756\t0.221046099400756\nleave_time_point_out\t10\tpseudotime_model\t1\t0.431851860991716\t0.431851860991716\nleave_time_point_out\t10\tindependent_celltype_clocks\t1\t0.325797533453726\t0.325797533453726\nleave_time_point_out\t10\tcoupled_regeneration_clock\t1\t0.0523145567011588\t0.0523145567011588\nleave_time_point_out\t2\tabsolute_time_model\t4\t0.0626940365480291\t0.0626940365480291\nleave_time_point_out\t2\tcell_state_model\t4\t0.171151127433812\t0.196364626403674\nleave_time_point_out\t2\tpseudotime_model\t4\t0.0584138306281223\t0.0678602445956348\nleave_time_point_out\t2\tindependent_celltype_clocks\t4\t0.0618490661975432\t0.0687922962180651\nleave_time_point_out\t2\tcoupled_regeneration_clock\t4\t0.12392623137884\t0.132243821136341\nleave_time_point_out\t21\tabsolute_time_model\t1\t0.429602669755451\t0.429602669755451\nleave_time_point_out\t21\tcell_state_model\t1\t0.504493725888172\t0.504493725888172\nleave_time_point_out\t21\tpseudotime_model\t1\t0.427831429046039\t0.427831429046039\nleave_time_point_out\t21\tindependent_celltype_clocks\t1\t0.438642366202728\t0.438642366202728\nleave_time_point_out\t21\tcoupled_regeneration_clock\t1\t0.526909367469673\t0.526909367469673\nleave_time_point_out\t3.5\tabsolute_time_model\t5\t0.164873587429995\t0.164873587429995\nleave_time_point_out\t3.5\tcell_state_model\t5\t0.0676891480234598\t0.0761423656157474\nleave_time_point_out\t3.5\tpseudotime_model\t5\t0.163015626986555\t0.163147040677077\nleave_time_point_out\t3.5\tindependent_celltype_clocks\t5\t0.097522531526032\t0.125813794310321\nleave_time_point_out\t3.5\tcoupled_regeneration_clock\t5\t0.136992466537759\t0.175374009492336\nleave_time_point_out\t5\tabsolute_time_model\t4\t0.316616308279557\t0.316616308279557\nleave_time_point_out\t5\tcell_state_model\t4\t0.236174253354891\t0.247128429927825\nleave_time_point_out\t5\tpseudotime_model\t4\t0.322559814763439\t0.322579142388439\nleave_time_point_out\t5\tindependent_celltype_clocks\t4\t0.160487761341075\t0.161838722521762\nleave_time_point_out\t5\tcoupled_regeneration_clock\t4\t0.115860526728337\t0.124867621421259\nleave_time_point_out\t0\tabsolute_time_model\t6\t0.713063502512058\t0.713063502512058\nleave_time_point_out\t0\tcell_state_model\t6\t0.779651777441278\t0.786502037068283\nleave_time_point_out\t0\tpseudotime_model\t6\t0.715063796090247\t0.715069740293546\nleave_time_point_out\t0\tindependent_celltype_clocks\t6\t0.723311545453045\t0.723400198195499\nleave_time_point_out\t0\tcoupled_regeneration_clock\t6\t0.813364663101677\t0.854599604451211\nleave_time_point_out\t7\tabsolute_time_model\t5\t0.490974479720516\t0.490974479720516\nleave_time_point_out\t7\tcell_state_model\t5\t0.489563939297816\t0.537658669095866\nleave_time_point_out\t7\tpseudotime_model\t5\t0.487246302718081\t0.487392529972543\nleave_time_point_out\t7\tindependent_celltype_clocks\t5\t0.4861172517606\t0.493351716380184\nleave_time_point_out\t7\tcoupled_regeneration_clock\t5\t0.484235291591212\t0.509456928715064\nleave_time_point_out\t4\tabsolute_time_model\t2\t0.202854660514046\t0.202854660514046\nleave_time_point_out\t4\tcell_state_model\t2\t0.0479644964229571\t0.0531127536282112\nleave_time_point_out\t4\tpseudotime_model\t2\t0.189925946913622\t0.189940780047289\nleave_time_point_out\t4\tindependent_celltype_clocks\t2\t0.228825167004365\t0.22883778303393\nleave_time_point_out\t4\tcoupled_regeneration_clock\t2\t0.246108117131211\t0.255064225684129\nleave_time_point_out\t1\tabsolute_time_model\t4\t0.284576364407778\t0.284576364407778\nleave_time_point_out\t1\tcell_state_model\t4\t0.0637377471616234\t0.0695597555275734\nleave_time_point_out\t1\tpseudotime_model\t4\t0.283217182581758\t0.283232371758495\nleave_time_point_out\t1\tindependent_celltype_clocks\t4\t0.267956229467059\t0.268291135545403\nleave_time_point_out\t1\tcoupled_regeneration_clock\t4\t0.631738973468277\t0.643119955625886\nleave_injury_model_out\tCTX\tabsolute_time_model\t7\t0.326695409924687\t0.366985949307812\nleave_injury_model_out\tCTX\tcell_state_model\t7\t0.281528669852217\t0.3457700600648\nleave_injury_model_out\tCTX\tpseudotime_model\t7\t0.391482215759864\t0.438891471867057\nleave_injury_model_out\tCTX\tindependent_celltype_clocks\t7\t0.269886406509354\t0.306272080913165\nleave_injury_model_out\tCTX\tcoupled_regeneration_clock\t7\t0.239203517211606\t0.297952860700845\nleave_injury_model_out\tNTX\tabsolute_time_model\t20\t0.277995000096154\t0.326470586689628\nleave_injury_model_out\tNTX\tcell_state_model\t20\t0.285370572801503\t0.368097802334326\nleave_injury_model_out\tNTX\tpseudotime_model\t20\t0.276412317900938\t0.32515779792146\nleave_injury_model_out\tNTX\tindependent_celltype_clocks\t20\t0.226222670496746\t0.304027573080484\nleave_injury_model_out\tNTX\tcoupled_regeneration_clock\t20\t0.276494521007222\t0.366118320347915\nleave_injury_model_out\tBaCl2\tabsolute_time_model\t6\t0.397439387029863\t0.42829985241144\nleave_injury_model_out\tBaCl2\tcell_state_model\t6\t0.389441694802323\t0.460953342459878\nleave_injury_model_out\tBaCl2\tpseudotime_model\t6\t0.399803911908376\t0.430484208616222\nleave_injury_model_out\tBaCl2\tindependent_celltype_clocks\t6\t0.407131136818283\t0.430869672267542\nleave_injury_model_out\tBaCl2\tcoupled_regeneration_clock\t6\t0.37220339734079\t0.40219051975007\nleave_age_out\told\tabsolute_time_model\t21\t0.291891646924985\t0.340346852757483\nleave_age_out\tgeriatric\tabsolute_time_model\t24\t0.287218055635684\t0.335689187827858\nleave_age_out\told\tcell_state_model\t21\t0.44837194702311\t0.559337146468721\nleave_age_out\tgeriatric\tcell_state_model\t24\t0.263540709417756\t0.323815617497756\nleave_age_out\told\tpseudotime_model\t21\t0.292019110310651\t0.340703357166813\nleave_age_out\tgeriatric\tpseudotime_model\t24\t0.287346152430587\t0.335935783471376\nleave_age_out\told\tindependent_celltype_clocks\t21\t0.250410996991935\t0.309632791573543\nleave_age_out\tgeriatric\tindependent_celltype_clocks\t24\t0.25023674084813\t0.292953436918401\nleave_age_out\told\tcoupled_regeneration_clock\t21\t0.371590586233208\t0.43643058050962\nleave_age_out\tgeriatric\tcoupled_regeneration_clock\t24\t0.219531238648937\t0.266490003333783\n",
"trajectory_inference/leave_study_out_sample_predictions.tsv"="held_out_study\tsample_uid\ttarget_tau\ttrajectory_tau\tabsolute_error\nGSE138826\tGSE138826::0.5 DPI\t0.194987500240385\t0.45139321521426\t0.256405714973875\nGSE138826\tGSE138826::10 DPI\t1\t0.712246053263185\t0.287753946736815\nGSE138826\tGSE138826::2 DPI\t0.528320833573719\t0.537309972120723\t0.00898913854700356\nGSE138826\tGSE138826::21 DPI\t1\t0.710895415046942\t0.289104584953058\nGSE138826\tGSE138826::3.5 DPI\t0.723308333814104\t0.692278986890168\t0.0310293469239361\nGSE138826\tGSE138826::5 DPI\t0.861654166907052\t0.662033318109482\t0.19962084879757\nGSE138826\tGSE138826::Noninjured\t0\t0.632927859811037\t0.632927859811037\nGSE143437\tGSE143437::D0_A\t0\t0.315730121820912\t0.315730121820912\nGSE143437\tGSE143437::D0_B\t0\t0.355310476908175\t0.355310476908175\nGSE143437\tGSE143437::D0_Cv3\t0\t0.352877640650724\t0.352877640650724\nGSE143437\tGSE143437::D2_C\t0.528320833573719\t0.188761659165643\t0.339559174408076\nGSE143437\tGSE143437::D2_D\t0.528320833573719\t0.207367666446551\t0.320953167127168\nGSE143437\tGSE143437::D5_A\t0.861654166907052\t0.434953705020069\t0.426700461886983\nGSE143437\tGSE143437::D5_B\t0.861654166907052\t0.433053340971882\t0.42860082593517\nGSE143437\tGSE143437::D5_C\t0.861654166907052\t0.464982119454321\t0.396672047452731\nGSE143437\tGSE143437::D7_C\t1\t0.406847370254669\t0.593152629745331\nGSE143437\tGSE143437::D7_D\t1\t0.512230700407581\t0.487769299592419\nGSE232106\tGSE232106::GSM7328463_Yng_D1_A\t0.333333333333333\t0.617842068260791\t0.284508734927458\nGSE232106\tGSE232106::GSM7328464_Yng_D1_B\t0.333333333333333\t0.669259166311285\t0.335925832977952\nGSE232106\tGSE232106::GSM7328465_Yng_D1_C\t0.333333333333333\t0.675895090264128\t0.342561756930795\nGSE232106\tGSE232106::GSM7328466_Yng_D1_D\t0.333333333333333\t0.683777738717163\t0.35044440538383\nGSE232106\tGSE232106::GSM7328467_Yng_D3_5_A\t0.723308333814104\t0.705222886864137\t0.0180854469499671\nGSE232106\tGSE232106::GSM7328468_Yng_D3_5_B\t0.723308333814104\t0.653828883730091\t0.0694794500840126\nGSE232106\tGSE232106::GSM7328469_Yng_D3_5_C\t0.723308333814104\t0.651001893584278\t0.0723064402298259\nGSE232106\tGSE232106::GSM7328470_Yng_D3_5_D\t0.723308333814104\t0.653828810259491\t0.0694795235546131\nGSE159500\tGSE159500::D2_Ev3\t0.528320833573719\t0.574147390489031\t0.0458265569153119\nGSE159500\tGSE159500::D7_Ev3\t1\t0.636235193163064\t0.363764806836936\nGSE205690\tGSE205690::animal_1\t0.773976031629121\t0.679507662610717\t0.0944683690184042\nGSE205690\tGSE205690::animal_2\t0.773976031629121\t0.671454350804544\t0.102521680824577\nGSE205690\tGSE205690::animal_3\t1\t0.692615751257893\t0.307384248742107\nGSE205690\tGSE205690::animal_4\t1\t0.692117685333272\t0.307882314666728\nGSE205690\tGSE205690::animal_6\t0\t0.621620852617703\t0.621620852617703\nGSE205690\tGSE205690::animal_5\t0\t0.595987916880639\t0.595987916880639\n",
"trajectory_inference/leave_study_out_trajectory_metrics.tsv"="held_out_study\tn_samples\tMAE_tau\tphase_order_correlation\nGSE138826\t7\t0.243690205820471\t0.846881214933825\nGSE143437\t10\t0.401732584552769\t0.687835939428795\nGSE232106\t8\t0.192848948879807\t-0.109108945117996\nGSE159500\t2\t0.204795681876124\t1\nGSE205690\t6\t0.338310897125026\t0.956182887467515\n",
"complete_dynamics/Kuramoto_coupling_by_age.tsv"="age_class\ttarget\tsource\tK\nyoung\tMuSC_myogenic\tmacrophage\t0.00938013944154654\nyoung\tMuSC_myogenic\tFAP\t-0.00425666425311504\nyoung\tMuSC_myogenic\tendothelial\t-0.0272926598025447\nyoung\tMuSC_myogenic\tmyonucleus\t-0.0143098611705874\nyoung\tMuSC_myogenic\tneutrophil\t-0.00537700825417604\nyoung\tMuSC_myogenic\tT_NK\t-0.00866463207572522\nyoung\tmacrophage\tMuSC_myogenic\t0.00269394788309208\nyoung\tmacrophage\tFAP\t-0.00309079863043509\nyoung\tmacrophage\tendothelial\t0.0304084051771009\nyoung\tmacrophage\tmyonucleus\t0.0801661273707027\nyoung\tmacrophage\tneutrophil\t0.0360690041512767\nyoung\tmacrophage\tT_NK\t0.0654543765148319\nyoung\tFAP\tMuSC_myogenic\t-0.0137435610626397\nyoung\tFAP\tmacrophage\t0.0483504268512643\nyoung\tFAP\tendothelial\t-0.0917685488296243\nyoung\tFAP\tmyonucleus\t-0.0953921591959769\nyoung\tFAP\tneutrophil\t-0.0355529153442037\nyoung\tFAP\tT_NK\t-0.0463297087145587\nyoung\tendothelial\tMuSC_myogenic\t0.0459124527222261\nyoung\tendothelial\tmacrophage\t0.0289580874806598\nyoung\tendothelial\tFAP\t0.0290464791088752\nyoung\tendothelial\tmyonucleus\t-0.00734676924927483\nyoung\tendothelial\tneutrophil\t0.0187053164495559\nyoung\tendothelial\tT_NK\t-0.0238818639090788\nyoung\tmyonucleus\tMuSC_myogenic\t0.00852409299884013\nyoung\tmyonucleus\tmacrophage\t-0.0222701725642923\nyoung\tmyonucleus\tFAP\t0.0361587588501557\nyoung\tmyonucleus\tendothelial\t0.00945687077746109\nyoung\tmyonucleus\tneutrophil\t0.0077809357527002\nyoung\tmyonucleus\tT_NK\t-0.016343878680363\nyoung\tneutrophil\tMuSC_myogenic\t0.00972739271441558\nyoung\tneutrophil\tmacrophage\t0.0179153706005112\nyoung\tneutrophil\tFAP\t0.0173860076097945\nyoung\tneutrophil\tendothelial\t-0.0182277256183512\nyoung\tneutrophil\tmyonucleus\t-0.0220236030611819\nyoung\tneutrophil\tT_NK\t-0.00635872281063753\nyoung\tT_NK\tMuSC_myogenic\t0.00977833548436706\nyoung\tT_NK\tmacrophage\t-0.0321054792811522\nyoung\tT_NK\tFAP\t-0.0389694205146017\nyoung\tT_NK\tendothelial\t0.0190996472625639\nyoung\tT_NK\tmyonucleus\t0.0166541522922107\nyoung\tT_NK\tneutrophil\t-0.0181966666136143\naged\tMuSC_myogenic\tmacrophage\t0.00507467226555366\naged\tMuSC_myogenic\tFAP\t0.048372533210672\naged\tMuSC_myogenic\tendothelial\t-0.063459885412421\naged\tMuSC_myogenic\tmyonucleus\t-0.0223011766312592\naged\tMuSC_myogenic\tneutrophil\t0.0491964513487464\naged\tMuSC_myogenic\tT_NK\t-0.0345300290778988\naged\tmacrophage\tMuSC_myogenic\t-0.00168707212447905\naged\tmacrophage\tFAP\t-0.0164342542410547\naged\tmacrophage\tendothelial\t0.0021640698998468\naged\tmacrophage\tmyonucleus\t0.0357735868588983\naged\tmacrophage\tneutrophil\t0.0194982235173792\naged\tmacrophage\tT_NK\t0.0255436733555464\naged\tFAP\tMuSC_myogenic\t-0.00549588650449057\naged\tFAP\tmacrophage\t0.00984118914738331\naged\tFAP\tendothelial\t0.0429062016507347\naged\tFAP\tmyonucleus\t0.0318457430700665\naged\tFAP\tneutrophil\t0.0353432994529808\naged\tFAP\tT_NK\t0.0084147317233609\naged\tendothelial\tMuSC_myogenic\t0.079796348016042\naged\tendothelial\tmacrophage\t0.0659952913705787\naged\tendothelial\tFAP\t0.0939728714557015\naged\tendothelial\tmyonucleus\t0.0325626466715988\naged\tendothelial\tneutrophil\t0.077088649023762\naged\tendothelial\tT_NK\t0.0530348156498358\naged\tmyonucleus\tMuSC_myogenic\t0.0116082138943197\naged\tmyonucleus\tmacrophage\t0.0217951510455048\naged\tmyonucleus\tFAP\t0.0248973583616096\naged\tmyonucleus\tendothelial\t-0.00669176811395261\naged\tmyonucleus\tneutrophil\t0.0127233946930515\naged\tmyonucleus\tT_NK\t0.00167399421018711\naged\tneutrophil\tMuSC_myogenic\t-0.0502503032797186\naged\tneutrophil\tmacrophage\t-0.0155380221214626\naged\tneutrophil\tFAP\t-0.040404366605183\naged\tneutrophil\tendothelial\t0.00505830366227607\naged\tneutrophil\tmyonucleus\t0.0081096487679337\naged\tneutrophil\tT_NK\t0.0155957613740532\naged\tT_NK\tMuSC_myogenic\t0.0092721311674226\naged\tT_NK\tmacrophage\t-0.0137424464597592\naged\tT_NK\tFAP\t-0.00421104246570606\naged\tT_NK\tendothelial\t-0.00268727788331958\naged\tT_NK\tmyonucleus\t-0.000509211700587572\naged\tT_NK\tneutrophil\t-0.0069629676810104\n")
rd<-function(f)fread(text=source_tables[[f]])
ink <- '#252B33'; muted <- '#69727C'; blue <- '#376B91'; teal <- '#338B83'
rust <- '#BA6048'; gold <- '#C49A48'; purple <- '#89719D'; light <- '#E9EDF0'
agepal <- c(young=blue, old=gold, geriatric=rust)
ctlab <- c(MuSC_myogenic='MuSC / myogenic',macrophage='Macrophage',FAP='FAP',
  endothelial='Endothelial',myonucleus='Myonucleus',neutrophil='Neutrophil',
  T_NK='T / NK',pericyte='Pericyte',other_major='Other major',other_immune='Other immune')
ctorder <- names(ctlab)
modelkeys <- c('cell_state_model','independent_celltype_clocks','coupled_regeneration_clock',
               'absolute_time_model','pseudotime_model')
modellab <- c('Cell state','Independent','Coupled','Training mean','Pseudotime')
names(modellab) <- modelkeys
modelpal <- setNames(c(teal,blue,rust,'#9CA5AC',purple),modelkeys)
base <- theme_classic(base_size=7, base_family='Helvetica') + theme(
  text=element_text(colour=ink),axis.text=element_text(size=6,colour=ink),
  axis.title=element_text(size=6.5),axis.line=element_line(linewidth=.3,colour=ink),
  axis.ticks=element_line(linewidth=.25,colour=ink),axis.ticks.length=unit(1.2,'mm'),
  plot.title=element_text(size=7,face='bold',margin=margin(b=5)),
  plot.subtitle=element_text(size=6,colour=muted,margin=margin(b=6)),
  plot.caption=element_text(size=5.5,colour=muted,hjust=0,margin=margin(t=5)),
  plot.tag=element_text(size=8,face='bold'),legend.title=element_text(size=6),
  legend.text=element_text(size=5.8),legend.key.size=unit(3,'mm'),
  legend.spacing.x=unit(1.5,'mm'),strip.background=element_blank(),
  strip.text=element_text(size=6,face='bold'),plot.margin=margin(6,7,5,6),
  panel.spacing=unit(3,'mm'))
theme_set(base)
heat <- theme(axis.line=element_blank(),axis.ticks=element_blank())
divscale <- function(lim, name) scale_fill_gradient2(low=blue,mid='#F5F5F3',high=rust,
  midpoint=0,limits=c(-lim,lim),breaks=c(-lim,0,lim),labels=scales::label_number(accuracy=.01),na.value='#E7E9EB',name=name,
  guide=guide_colourbar(barwidth=unit(23,'mm'),barheight=unit(2,'mm'),title.position='top'))
savefig <- function(p, stem, height=155, width=183){
 p<-p+plot_annotation(tag_levels='a') & theme(plot.tag=element_text(face='bold',size=8))
 pdf_device(file.path(out,paste0(stem,'.pdf')),width=width/25.4,height=height/25.4)
 print(p);dev.off()
}
factct <- function(x) factor(x,levels=rev(ctorder),labels=rev(unname(ctlab)))

# FIGURE 1 — observed coverage, counts and a conceptual workflow.
reg <- fread(file.path(root,'config/datasets.tsv'))
flow <- data.table(x=1:4, label=c('Public data','Lineage states','Expression phase','Relative timing'),
  detail=c('Mouse injury + age\nHuman + spatial transfer','Myogenic / immune\nStromal / vascular',
           'Six expression programmes\nExpected young phase','Phase deviation\nRegenerative synchrony'))
p1a <- ggplot(flow,aes(x,1))+
  geom_segment(data=data.table(x=1:3),aes(x=x+.37,xend=x+.63,y=1,yend=1),inherit.aes=FALSE,
               arrow=arrow(length=unit(1.5,'mm')),linewidth=.35,colour=muted)+
  geom_label(aes(label=label),size=2.45,fontface='bold',fill='#EFF3F5',linewidth=.3,
             label.padding=unit(2,'mm'))+
  geom_text(aes(y=.63,label=detail),size=2.15,lineheight=1.2)+
  coord_cartesian(xlim=c(.5,4.5),ylim=c(.35,1.24),clip='off')+theme_void()+
  labs(title='A study-aware framework for relative regenerative timing',
       subtitle='516,633 non-overlapping QC profiles  |  42 condition-time combinations  |  80 dynamic samples')+
  theme(plot.title=element_text(size=7,face='bold'),plot.subtitle=element_text(size=6),
        plot.margin=margin(8,8,4,8))
cov <- reg[grepl('^GSE',accession)]
cov <- rbindlist(lapply(seq_len(nrow(cov)),function(i){
 d<-cov[i]; t<-strsplit(d$timepoints_dpi,';',fixed=TRUE)[[1]]
 if(d$accession=='GSE214544')t<-c('Pre','3 h')
 data.table(accession=d$accession,time=t,
   domain=ifelse(d$organism=='human','Human',ifelse(grepl('spatial',d$role),'Spatial','Mouse')))
}))
tl <- c('0','0.5','1','2','3.5','4','5','7','10','14','21','Pre','3 h','baseline')
cov[,time:=factor(time,levels=tl)]
cov[,accession:=factor(accession,levels=rev(reg[grepl('^GSE',accession)]$accession))]
p1b<-ggplot(cov,aes(time,accession,fill=domain))+geom_tile(width=.8,height=.72)+
  scale_fill_manual(values=c(Mouse=blue,Spatial=purple,Human=rust))+
  labs(title='Public study coverage',x='Days after injury / human sampling',y=NULL,fill=NULL)+
  heat+theme(legend.position='bottom',axis.text.x=element_text(size=5.8,angle=45,hjust=1))
cnt<-rd('audited_cell_counts.tsv');cnt[,label:=ifelse(dataset=='reused_core','Core atlas (unique)',dataset)]
cnt[,label:=factor(label,levels=label[order(final_cells)])]
p1c<-ggplot(cnt,aes(final_cells/1000,label))+geom_segment(aes(x=0,xend=final_cells/1000,yend=label),
    linewidth=2.5,colour='#D8E1E8')+geom_point(size=2.2,colour=blue)+
  geom_text(aes(label=format(final_cells,big.mark=',')),hjust=-.18,size=2.05)+
  scale_x_continuous(expand=expansion(mult=c(.02,.3)))+
  labs(title='Audited profiles without double counting',x='Cells / nuclei (thousands)',y=NULL,
    caption='Accession counts are not independent study counts.')+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
stages<-data.table(i=1:6,phase=c('Quiescence','Activation','Proliferation','Differentiation','Fusion','Remodelling'))
p1d<-ggplot(stages,aes(i,1))+annotate('segment',x=1,xend=6,y=1,yend=1,colour='#CBD2D7',linewidth=.7)+
 geom_point(aes(fill=phase),shape=21,size=3.2,stroke=.3)+
 geom_text(aes(y=.76,label=phase),size=2.05)+
 annotate('text',x=3.5,y=.37,label='MuSC / myogenic   /   immune   /   FAP   /   vascular   /   myonuclei',size=2.3)+
 scale_fill_manual(values=setNames(c('#D6DFE6',blue,teal,gold,rust,purple),stages$phase),guide='none')+
 coord_cartesian(xlim=c(.6,6.4),ylim=c(.2,1.25),clip='off')+theme_void()+
 labs(title='Biological coverage of the regeneration ontology',subtitle='Ordered stages; no periodicity or inferred peak locations implied')+
 theme(plot.title=element_text(size=7,face='bold'),plot.subtitle=element_text(size=6))
savefig(p1a / (p1b|p1c) / p1d + plot_layout(heights=c(.75,2.1,.7)), 'Figure1_overview',165)
sdout(cov,'Figure1b_coverage');sdout(cnt,'Figure1c_counts')

# FIGURE 2 — error landscape, study-level gain and aggregate ranking.
m<-rd('complete_dynamics/all_validation_metrics.tsv')
m[,model:=factor(model,levels=modelkeys)]
ls<-m[validation=='leave_study_out']
ls[,study:=sprintf('%s\n(n = %d)',held_out,n_samples)]
p2a<-ggplot(ls,aes(model,study,fill=MAE_tau))+
 geom_tile(colour='white',linewidth=1)+geom_text(aes(label=sprintf('%.3f',MAE_tau)),size=2.4,colour=ink)+
 scale_x_discrete(labels=modellab)+scale_fill_gradient(low='#F3F6F7',high='#6C9CB5',limits=c(0,.55),name='MAE',
   guide=guide_colourbar(barheight=unit(2,'mm'),barwidth=unit(23,'mm'),title.position='top'))+
 labs(title='Prediction error across held-out studies',x=NULL,y=NULL,
   subtitle='Numbers are MAE in normalized injury-time target; lower is better')+heat+
 theme(legend.position='bottom',axis.text.x=element_text(angle=30,hjust=1))
w<-dcast(ls,held_out~model,value.var='MAE_tau');w[,gain:=independent_celltype_clocks-coupled_regeneration_clock]
w[,held_out:=factor(held_out,levels=held_out[order(gain)])]
p2b<-ggplot(w,aes(gain,held_out))+geom_vline(xintercept=0,colour='#AEB7BD',linetype=2,linewidth=.3)+
 geom_segment(aes(x=0,xend=gain,yend=held_out),colour='#CBD3D9',linewidth=1.1)+
 geom_point(aes(fill=gain>0),shape=21,size=2.6,stroke=.3)+
 geom_text(aes(label=sprintf('%+.3f',gain)),hjust=ifelse(w$gain>0,-.3,1.3),size=2.15)+
 scale_fill_manual(values=c('TRUE'=teal,'FALSE'=rust),guide='none')+
 scale_x_continuous(limits=c(-.23,.10),breaks=c(-.2,-.1,0,.05))+
 labs(title='Coupling helps in 3 of 5 studies',x='Independent MAE - coupled MAE',y=NULL,
   caption='Positive values favour coupling.')
rank<-m[validation=='leave_study_out',.(sample_weighted_MAE=weighted.mean(MAE_tau,n_samples)),by=model][order(sample_weighted_MAE)];rank[,model:=factor(model,levels=rev(modelkeys))]
p2c<-ggplot(rank,aes(sample_weighted_MAE,model,fill=model))+
 geom_segment(aes(x=.20,xend=sample_weighted_MAE,yend=model),colour='#D6DDE1',linewidth=1)+
 geom_point(shape=21,size=2.6,stroke=.25)+
 geom_text(aes(label=sprintf('%.3f',sample_weighted_MAE)),hjust=-.25,size=2.15)+
 scale_y_discrete(labels=modellab)+scale_fill_manual(values=modelpal,guide='none')+
 scale_x_continuous(limits=c(.20,.35),breaks=c(.20,.25,.30,.35))+
 labs(title='Cell state has the lowest pooled error',x='Sample-weighted MAE',y=NULL)
stress<-m[validation!='leave_study_out',.(MAE=weighted.mean(MAE_tau,n_samples)),by=.(validation,model)]
stress[,domain:=factor(validation,levels=c('leave_age_out','leave_injury_model_out','leave_time_point_out'),
                       labels=c('Age','Injury model','Time point'))]
p2d<-ggplot(stress,aes(domain,MAE,colour=model,group=model))+
 geom_line(linewidth=.4)+geom_point(aes(shape=model),size=1.7)+
 scale_colour_manual(values=modelpal,labels=modellab)+scale_shape_manual(values=15:19,labels=modellab)+
 labs(title='Sensitivity to the held-out domain',x='Held-out group',y='Sample-weighted MAE',colour=NULL,shape=NULL)+
 theme(legend.position='bottom',legend.text=element_text(size=5.5))+guides(colour=guide_legend(nrow=2),shape=guide_legend(nrow=2))
savefig((p2a|p2b)/(p2c|p2d)+plot_layout(heights=c(1.1,1)), 'Figure2_validation',158)
sdout(ls,'Figure2a_leave_study');sdout(w,'Figure2b_coupling_gain');sdout(rank,'Figure2c_ranking');sdout(stress,'Figure2d_stress')

tr<-rd('trajectory_inference/leave_study_out_sample_predictions.tsv')
sw<-rd('trajectory_inference/leave_study_out_trajectory_metrics.tsv')
sdout(tr,'Supplementary2_trajectory_predictions');sdout(sw,'Supplementary2_trajectory_metrics')
dskeys<-sort(unique(tr$held_out_study));dspal<-setNames(c(blue,teal,gold,rust,purple)[seq_along(dskeys)],dskeys)
ps2a<-ggplot(tr,aes(target_tau,trajectory_tau,colour=held_out_study))+
 geom_abline(slope=1,intercept=0,linetype=2,linewidth=.3,colour=muted)+geom_point(size=1.8)+
 facet_wrap(~held_out_study,nrow=1)+coord_fixed(xlim=c(0,1),ylim=c(0,1))+
 scale_colour_manual(values=dspal,guide='none')+scale_x_continuous(breaks=c(0,1))+
 scale_y_continuous(breaks=c(0,.5,1))+
 labs(title='Fold-refitted Slingshot transfer',x='Anchored injury-time target',y='Held-out trajectory phase',
      subtitle='33 held-out samples; trajectories learned inside each training fold')
terr<-tr[,.(MAE=mean(absolute_error),n_samples=.N),by=held_out_study]
ps2b<-ggplot(terr,aes(MAE,reorder(held_out_study,MAE),fill=held_out_study))+
 geom_segment(aes(x=0,xend=MAE,yend=reorder(held_out_study,MAE)),colour=light,linewidth=1.5)+
 geom_point(shape=21,size=2.5,stroke=.3)+geom_text(aes(label=sprintf('%.3f (n=%d)',MAE,n_samples)),hjust=-.2,size=2.2)+
 scale_fill_manual(values=dspal,guide='none')+scale_x_continuous(limits=c(0,.55))+
 labs(title='Transfer error varies by study',x='MAE',y=NULL)
kk<-rd('complete_dynamics/Kuramoto_coupling_by_age.tsv')
ka<-kk[,.(mean_abs_K=mean(abs(K))),by=age_class]
ps2c<-ggplot(kk,aes(age_class,abs(K),fill=age_class))+
 geom_boxplot(width=.42,outlier.shape=NA,linewidth=.3,alpha=.5)+
 geom_point(shape=21,size=1.25,stroke=.2,position=position_jitter(width=.13,seed=22),alpha=.6)+
 geom_text(data=ka,aes(y=.13,label=sprintf('Mean |K| = %.4f',mean_abs_K)),size=2.2,inherit.aes=TRUE)+
 scale_fill_manual(values=c(young=blue,aged=rust),guide='none')+
 labs(title='No global young-coupling advantage',x=NULL,y='Absolute fitted coupling |K|',
      caption='Directed lineage pairs are coefficients, not biological replicates.')
savefig(ps2a/(ps2b|ps2c)+plot_layout(heights=c(.9,1.1)),'SupplementaryFigure2',142)
sdout(kk,'Supplementary2c_coupling')

})
