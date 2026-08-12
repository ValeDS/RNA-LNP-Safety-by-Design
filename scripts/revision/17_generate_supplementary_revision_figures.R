suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(scales)})
source(file.path("scripts", "config.R"))
inp <- revision_output_dir
out <- figure_output_dir
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cols <- c(Low="#2B8CBE",Intermediate="#F4A340",High="#C0392B")
theme_pub <- theme_classic(base_size=11)+theme(plot.title=element_text(face="bold",size=11),plot.subtitle=element_text(size=9),legend.position="bottom",strip.background=element_rect(fill="#F2F2F2",colour=NA),strip.text=element_text(face="bold"))

# Figure S1: SbD weight sensitivity.
w <- read.csv(file.path(inp,"10_SbD_weight_set_summary.csv"))
cand <- read.csv(file.path(inp,"10_SbD_candidate_robustness.csv")) %>% arrange(original_rank) %>% slice_head(n=12) %>% mutate(lnp_id=factor(lnp_id,levels=rev(lnp_id)))
p1a <- ggplot(w,aes(spearman_rank_vs_original))+geom_histogram(binwidth=.02,boundary=0,fill="#3973AC",colour="white")+geom_vline(xintercept=median(w$spearman_rank_vs_original),linetype=2)+labs(title="A  Global rank stability",x="Spearman correlation with original ranking",y="Weight sets")+theme_pub
p1b <- ggplot(cand,aes(y=lnp_id,x=median_rank,xmin=min_rank,xmax=max_rank))+geom_errorbar(width=.18,colour="#888888",orientation="y")+geom_point(aes(size=top10_frequency,colour=top1_frequency))+scale_colour_viridis_c(labels=percent)+scale_size_continuous(labels=percent,range=c(2,6))+labs(title="B  Candidate rank robustness",x="Rank across 3,876 weight sets",y=NULL,colour="Top-1 frequency",size="Top-10 frequency")+theme_pub+theme(legend.position="right")
fig1 <- p1a+p1b+plot_annotation(title="Supplementary Figure S1. Safety-by-Design weight sensitivity",theme=theme(plot.title=element_text(margin=margin(4,4,10,4))))
ggsave(file.path(out,"Figure_S1_SbD_weight_sensitivity.png"),fig1,width=12,height=5.8,dpi=300,bg="white")
ggsave(file.path(out,"Figure_S1_SbD_weight_sensitivity.pdf"),fig1,width=12,height=5.8)

# Figure S2: balanced risk analysis.
rs <- read.csv(file.path(inp,"16_balanced_risk_endpoint_summary.csv")) %>% mutate(inflammatory_risk_class=factor(inflammatory_risk_class,levels=c("Low","Intermediate","High")))
rl <- rs %>% select(inflammatory_risk_class,delivery,adaptive,innate,cytokine,off_target,sbd) %>% pivot_longer(-inflammatory_risk_class,names_to="endpoint",values_to="value") %>% group_by(endpoint) %>% mutate(relative_to_low=value/value[inflammatory_risk_class=="Low"]) %>% ungroup() %>% mutate(endpoint=recode(endpoint,delivery="Delivery",adaptive="Adaptive",innate="Innate",cytokine="Cytokine",off_target="Off-target",sbd="SbD"))
p2a <- ggplot(rl,aes(inflammatory_risk_class,relative_to_low,group=endpoint,colour=endpoint))+geom_hline(yintercept=1,colour="#BBBBBB")+geom_line(linewidth=.7)+geom_point(size=2)+labs(title="A  Endpoint change relative to Low risk",x=NULL,y="Ratio to Low-risk mean",colour="Endpoint")+theme_pub
rr <- read.csv(file.path(inp,"16_balanced_risk_formulation_rankings_pareto.csv")) %>% filter(sbd_rank<=10) %>% mutate(inflammatory_risk_class=factor(inflammatory_risk_class,levels=c("Low","Intermediate","High")))
p2b <- ggplot(rr,aes(inflammatory_risk_class,mean_safety_score,group=lnp_id,colour=lnp_id))+geom_line(alpha=.65)+geom_point(size=1.8)+geom_line(data=filter(rr,lnp_id=="LNP_0462"),linewidth=1.4,colour="black")+geom_point(data=filter(rr,lnp_id=="LNP_0462"),size=3,colour="black")+labs(title="B  Leading formulations across risk strata",subtitle="Black line: LNP_0462",x=NULL,y="Mean Safety-by-Design score",colour="Formulation")+theme_pub+theme(legend.position="none")
fig2 <- p2a+p2b+plot_annotation(title="Supplementary Figure S2. Balanced inflammatory-risk robustness",theme=theme(plot.title=element_text(margin=margin(4,4,10,4))))
ggsave(file.path(out,"Figure_S2_balanced_risk_robustness.png"),fig2,width=12,height=6.2,dpi=300,bg="white",limitsize=FALSE)
ggsave(file.path(out,"Figure_S2_balanced_risk_robustness.pdf"),fig2,width=12,height=6.2)

# Figure S3: BPCI null and donor bootstrap uncertainty.
null <- read.csv(file.path(inp,"09_BPCI_exact_null_summary.csv"))
boot <- read.csv(file.path(inp,"09_BPCI_donor_bootstrap_summary.csv"))
p3a <- ggplot(null,aes(factor(day),observed_BPCI))+geom_linerange(aes(ymin=null_q025,ymax=null_q975),linewidth=4,colour="#D9D9D9")+geom_point(size=2.5,colour="#2166AC")+geom_point(aes(y=null_mean),shape=4,size=2.5,colour="#B2182B")+labs(title="A  Observed BPCI versus exact shuffled-label null",subtitle="Grey: central 95% null interval; red cross: null mean",x="Post-vaccination day",y="BPCI")+coord_cartesian(ylim=c(0,1))+theme_pub
p3b <- ggplot(boot,aes(factor(day),bootstrap_median_BPCI))+geom_linerange(aes(ymin=bootstrap_q025_BPCI,ymax=bootstrap_q975_BPCI),linewidth=.8,colour="#555555")+geom_point(aes(size=n_donors),colour="#2166AC")+geom_point(aes(y=original_unpaired_BPCI),shape=1,size=3,colour="#B2182B")+scale_size_continuous(breaks=c(3,6))+labs(title="B  Paired-donor bootstrap uncertainty",subtitle="Red open point: original unpaired estimate",x="Post-vaccination day",y="BPCI",size="Donors")+coord_cartesian(ylim=c(0,1))+theme_pub
fig3 <- p3a+p3b+plot_annotation(title="Supplementary Figure S3. BPCI null comparison and uncertainty")
ggsave(file.path(out,"Figure_S3_BPCI_null_bootstrap.png"),fig3,width=11,height=5.7,dpi=300,bg="white")
ggsave(file.path(out,"Figure_S3_BPCI_null_bootstrap.pdf"),fig3,width=11,height=5.7)
message("Supplementary figures generated.")
