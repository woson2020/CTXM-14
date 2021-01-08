library(reshape2)
library(plyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(readr);
library(scales);
library(RColorBrewer);
library("scales")
##Fig.2A
Figure.2A<- fitness_caz_5_1mic_ctxm %>%
  mutate(grp = ceiling(pos / 792)) %>%
  mutate(is100=ifelse(effect_size>100 & gt!="mean","A",NA))%>%
  mutate(new_pos=pos+81)%>%
  mutate(s=10^effect_size - 1) %>%
  ggplot(aes(x=gt,y=new_pos,fill=log10(effect_size)))+
  geom_tile(width=1,height=1)+#"brown","white","green" "blue","white","red""#FFFF00"huang
  scale_fill_gradientn(name="Relative growth",colours = c("blue","white","red"),
                       na.value = "white",limits=c(-1,1),breaks=c(-1,0,1),oob=squish,
                       labels=c(TeX("$10^{-1}$"),TeX("$10^{0}$"),TeX("$10^{1}$")))+
  xlab("")+ylab("")+theme_bw()+
  scale_y_reverse(expand = c(0,0),breaks = seq(1+81,792+81,49),position = "left")+
  scale_x_discrete(expand = c(0,0))+ 
  facet_wrap(~grp,ncol=1,scales = "free_y")+
  geom_point( aes( x = gt, y =new_pos, alpha=is.na(is100),shape=is100),size=1,show.legend=T)+
  guides(fill = guide_colorbar(order = 0,override.aes = list(size = 1),barwidth=12.5,keywidth=6.5,title.position ="top"))+
  scale_shape_manual("",values=c(21,10),label=c("",""),guide="none")+
  scale_alpha_manual(values=c("TRUE"=0,"FALSE"=1),guide="none") +
  theme(legend.position = "top",
        legend.title = element_text(size = 20,hjust=0.5), 
        legend.text = element_text(size = 18),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))



##Fig.2B
Figure.2B<- fitness_ctx_5_1mic_ctxm %>%
  mutate(grp = ceiling(pos / 792)) %>%
  mutate(is100=ifelse(effect_size>100 & gt!="mean","A",NA))%>%
  mutate(new_pos=pos+81)%>%
  mutate(s=10^effect_size - 1) %>%
  ggplot(aes(x=gt,y=new_pos,fill=log10(effect_size)))+
  geom_tile(width=1,height=1)+#"brown","white","green" "blue","white","red""#FFFF00"huang
  scale_fill_gradientn(name="Relative growth",colours = c("#8B658B","white","green"),
                       na.value = "white",limits=c(-1,1),breaks=c(-1,0,1),oob=squish,
                       labels=c(TeX("$10^{-1}$"),TeX("$10^{0}$"),TeX("$10^{1}$")))+
  xlab("")+ylab("")+theme_bw()+
  scale_y_reverse(expand = c(0,0),breaks = seq(1+81,792+81,49),position = "right")+
  scale_x_discrete(expand = c(0,0))+ 
  facet_wrap(~grp,ncol=1,scales = "free_y")+
  geom_point( aes( x = gt, y =new_pos, alpha=is.na(is100),shape=is100),size=1,show.legend=T)+
  guides(fill = guide_colorbar(order = 0,override.aes = list(size = 1),barwidth=12.5,keywidth=6.5,title.position ="top"))+
  scale_shape_manual("",values=c(21,10),label=c("",""),guide="none")+
  scale_alpha_manual(values=c("TRUE"=0,"FALSE"=1),guide="none") +
  theme(legend.position = "top",
        legend.title = element_text(size = 20,hjust=0.5), 
        legend.text = element_text(size = 18),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))


##Fig.2F
Figure.2F <- growth_plot3 %>%
  filter(sd!=0) %>%
  ggplot(aes(x=od,y=log10(mean),color=mic,lty=drug))+
  geom_point(aes(group=drug,shape=drug),stat="identity",size=1.5)+
  geom_line()+
  geom_errorbar (aes(ymin = log10(mean - sd),ymax =log10(mean + sd)), width=.05,linetype="solid")+
  scale_x_continuous(breaks = c(1,2,3),labels=c("T1 (0.23)","T2 (0.43)","T3 (0.75)"),limits=c(0.7,3.3))+
  scale_y_continuous(breaks = c(0.5,1,1.5),labels = c(TeX("$10^{0.5}$"),TeX("$10^{1.0}$"),TeX("$10^{1.5}$")),limits=c(0,1.8))+
  theme_classic()+
  scale_linetype_manual("Antibiotic",values=c("solid","dotted"),labels=c("Ceftazidime","Cefotaxime"))+
  scale_color_manual("Concentration",values=fill_snr)+
  scale_shape_manual("Antibiotic",values=c(16,17),labels=c("Ceftazidime","Cefotaxime"))+
  guides(color= guide_legend(order = 2),shape= guide_legend(order =1),lty= guide_legend(order =1),nrow = 1)+
  ylab("Relative growth")+
  xlab("Time (OD)")+
  theme( legend.title = element_text(face="bold.italic",size=18),
         legend.text = element_text(size=16),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.text.x=element_text(size=20),
         axis.text.y=element_text(size=20),
         axis.title.x=element_text(size=24),
         axis.title.y=element_text(size=24))

##Fig.2G
Figure.2G <- ggplot(clinic_enrich,aes(x=cut,y=odd))+
  geom_bar(position=position_dodge(0.9), width=.5,stat="identity",fill="#6699CC")+
  geom_errorbar(aes(ymin = odd-pvalue/sqrt(1),ymax = odd + pvalue/sqrt(1)), 
                width=.2,linetype="solid",alpha=1,
                position = position_dodge(0.9),show.legend = F)+
  theme_classic()+
  geom_hline(yintercept = c(1),color="gray",linetype="longdash",size=1)+
  ylab("Fold enrichment for \nclinical isolates")+
  geom_text(aes(label=p,y=odd+pvalue+1),size=8)+
  xlab("Threshold for relative growth (≥ x)")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)),breaks=c(0,1,5,10,15))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    axis.title.x=element_text(size=21),
    axis.title.y=element_text(size=21))

save(Figure.2A,Figure.2B,Figure.2F,Figure.2G,file = "/mnt/data2/disk/smrtanalysis/pacbio_data_new/new_illumina_data_2-14/20202015_sample_data/gfp_plot_data/CTXM_Figure2.RData")


