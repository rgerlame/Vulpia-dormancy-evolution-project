require(ggplot2)
require(car)
require(lmerTest)
require(visreg)
require(dplyr)
require(lmodel2)
require(glmmTMB)
require(lsmeans)
require(cowplot)
require(DHARMa)
require(ggeffects)
require(lsmeans)

#data

envt<-read.csv("./data/env.means.all.csv")
germ1<-read.csv("./data/germ from lambda expt.csv")
germ2<-read.csv("./data/germ from trait expt.csv")
awn<-read.csv("./data/awn_data.csv")
hair<-read.csv("./data/seedhair_zeroes_R.csv")
adhere<-read.csv("./data/deerdata_final.csv")


###########
# Figure 2 
###########

# Panel A (i)

germ1.envt<-merge(germ1,envt,by.x="pop",by.y="site")
germ1.envt$n<-7

germ1.envt<-na.omit(germ1.envt)

lm.dorm_prod<-glmmTMB(prop.dorm~productivity+(1|pop),data=germ1.envt,weights=n,family="binomial")
Anova(lm.dorm_prod)
germ1.envt$fit<-fitted.values(lm.dorm_prod)

vis<-ggpredict(lm.dorm_prod, terms = c("productivity[all]"), type="fe")


#fig
ggplot(vis,aes(x=x,y=predicted)) + 
  geom_ribbon(data=vis,aes(ymin=conf.low, ymax=conf.high, x=x), fill="orange", alpha=0.2) +
  geom_boxplot(data=germ1.envt, aes(y=prop.dorm,x=productivity,group=pop), fill="lightgrey",color="grey") + 
  theme_classic() + 
  geom_line(aes(x=x,y=predicted), color="orange", size=1) + 
  labs(x = "Site productivity", y = "Dormancy rate") +
  geom_point(data=germ1.envt, aes(y=fit,x=productivity,group=pop),pch=21,fill="cornflowerblue",size=3) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) 

ggsave(height=5,width=5,"./figures/fig 1 - dormancy.pdf")

# Panel A (ii)

lm.dorm_envt<-glmmTMB(prop.dorm~percent.sm+ph+(1|pop),data=germ1.envt,weights=n,family="binomial")
Anova(lm.dorm_envt)

#fig
visreg2d(fit=lm.dorm_envt,xvar="percent.sm",yvar="ph",plot.type="persp",scale="response",xlab=" ",ylab=" ",zlab=" ",col="cornflowerblue",zlim=c(0.4,1))


# Panel B (i)

adhere.envt<-merge(adhere,envt,by.x="pop",by.y="site")
adhere.envt<-na.omit(adhere.envt)

adhere.envt<-merge(awn,adhere.envt,by="id")
dispersal.data<-merge(adhere.envt,hair,by="id")
dispersal.data$total.adhere1<-1-dispersal.data$total.adhere

lm.disp_prod<-glmmTMB(total.adhere1~ productivity + (1|pop),data=dispersal.data,family="binomial",weights=total_seed)
Anova(lm.disp_prod)

dispersal.data$fit<-fitted.values(lm.disp_prod)

vis<-ggpredict(lm.disp_prod, terms = c("productivity[all]"), type="fe")


#fig
ggplot(vis,aes(x=x,y=predicted)) + 
  geom_ribbon(data=vis,aes(ymin=conf.low, ymax=conf.high, x=x), fill="orange", alpha=0.2) +
  geom_boxplot(data=dispersal.data, aes(y=total.adhere1,x=productivity,group=pop), fill="lightgrey",color="grey") + 
  theme_classic() + 
  geom_line(aes(x=x,y=predicted), color="orange", size=1) + 
  labs(x = "Site productivity", y = "Dispersal ability") +
  geom_point(data=dispersal.data, aes(y=fit,x=productivity,group=pop),pch=21,fill="cornflowerblue",size=3) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) 


ggsave(height=5,width=5,"./figures/fig 1 - disp ability.pdf")

#AIC comparison
lm.pop<-glmer(total.adhere1~productivity + (1|pop),data=dispersal.data,family="binomial",weights=total_seed)
lm.nopop<-glm(total.adhere1~productivity,data=dispersal.data,family="binomial",weights=total_seed)
anova(lm.pop,lm.nopop)


# Panel B (ii)

#stats below: glmTMB converges but is not compatible with visreg. glmer is used for plotting only, glmmtmb for actual stats

lm.disp_envt<-glmmTMB(total.adhere1~percent.sm+ph + (1|pop),data=dispersal.data,family="binomial",weights=total_seed)
Anova(lm.disp_envt)

lm<-glmer(total.adhere1~percent.sm+ph + (1|pop),data=dispersal.data,family="binomial",weights=total_seed)
Anova(lm)

#fig
visreg2d(fit=lm,xvar="percent.sm",yvar="ph",plot.type="persp",scale="response",xlab=" ",ylab=" ",zlab=" ",col="cornflowerblue")


##########
#Figure 3
##########

micro.site<-read.csv("./data/microsite availability.csv")
micro.site$weight<-25
micro.site$prop<-micro.site$prop.bad/25

lm.site<-glmer(prop~treatment+(1|site),data=micro.site,weights=weight,family="binomial")
Anova(lm.site)

lsmeans(lm.site,"treatment",type="response")
0.3921009/0.1205309

vis<-ggpredict(lm.site, terms = c("treatment"), type="fe")

vis$x<-factor(vis$x,levels=c("harsh","benign"))

fig3a<-ggplot(data=vis,aes(x=x,y=predicted)) + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), colour="black", width=.05)+
  geom_point(fill="cornflowerblue", pch=21, size = 5) +
  theme_classic() + labs(x = "Site productivity", y = "Proportion of microsites unsuitable to germination") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) + 
  scale_x_discrete(labels=c("low","high")) 

ggsave(height=5,width=5,"./figures/microsite availability.pdf")


############
#Figure 4
############

sim<-read.csv("./data/plastic dormancy, 2 seeds.csv")


df.means<-sim %>% group_by(dormancy.rate, env) %>%
  filter(time>90) %>%
  summarise(mean=mean(n))

ggplot(data=df.means,aes(x=1-env,y=mean,color=dormancy.rate, group=dormancy.rate)) +
  geom_point() + 
  geom_line() +
  theme_classic() +
  labs(y = "Avg. N at equilibrium", x = "Proportion of microsites unsuitable to germination", color="Dorm. rate") +
  geom_vline(xintercept = c(0.12,0.39), linetype="dashed", color="grey") +
  theme(legend.position = c(0.9,0.6)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
  annotate("text", x = 0.12, y = 0, label = "high") +
  annotate("text", x = 0.39, y = 0, label = "low") +
  scale_color_continuous(limits = c(0,0.9), breaks = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

ggsave(height=4,width=8,"./figures/simulation.pdf")



###########
# Table 1
##########

###dormancy covariates

mean.dorm1<-germ1 %>%
  group_by(pop) %>%
  summarise(mean.dorm1=mean(prop.dorm))

mean.hair<-hair %>%
  group_by(pop) %>%
  summarise(mean.hair=mean(avg_hair, na.rm = TRUE))

mean.awn<-awn %>%
  group_by(pop) %>%
  summarise(mean.awn=mean(avg_awn, na.rm = TRUE),mean.seed=mean(avg_seed, na.rm = TRUE))

germ1.traits<-merge(mean.dorm1,mean.hair,by="pop")
germ1.traits<-merge(germ1.traits,mean.awn,by="pop")

#ggplot(germ1.traits,aes(y=mean.dorm1,x=mean.hair,color=pop)) + geom_point() 
#ggplot(germ1.traits,aes(y=mean.dorm1,x=mean.awn,color=pop)) + geom_point() 
#ggplot(germ1.traits,aes(y=mean.dorm1,x=mean.seed,color=pop)) + geom_point() 

lm.dorm_traits<-lm(mean.dorm1~mean.hair+mean.awn+mean.seed,data=germ1.traits)
#Anova(lm.dorm_traits)
summary(lm.dorm_traits)

###dispersal covariates

hist(dispersal.data$ratio1) 
hist(dispersal.data$ratio2)

dispersal.data$avg_seed<-dispersal.data$avg_seed-dispersal.data$avg_awn

#need to adjust # seeds for ratio2
dispersal.data$total.attached<-dispersal.data$total_seed-(dispersal.data$total_seed*dispersal.data$total.adhere)

#ratio 1 (i.e., getting initially attached): both seed traits affect it, but hair most strongly
lm<-glmmTMB(ratio1~avg_hair+avg_awn + avg_seed + (1|pop),data=dispersal.data,family="binomial")
#Anova(lm)
summary(lm)
visreg2d(fit=lm,xvar="avg_hair",yvar="avg_awn",plot.type="persp",scale="response")

#ratio 2 (i.e., staying attached): only hair
lm<-glmmTMB(ratio2~avg_hair+avg_awn+avg_seed + (1|pop),data=dispersal.data,family="binomial",weights=total.attached)
#Anova(lm)
summary(lm)
visreg2d(fit=lm,xvar="avg_hair",yvar="avg_awn",plot.type="persp",scale="response")

#total adhere: mostly hair, marginal awn
lm<-glmmTMB(total.adhere~avg_hair+avg_awn+avg_seed + (1|pop),data=dispersal.data,family="binomial",weights=total.attached)
#Anova(lm)
summary(lm)
visreg2d(fit=lm,xvar="avg_hair",yvar="avg_awn",plot.type="persp",scale="response")



##########
#Figure 5
##########

require(lmodel2)

mean.dorm1<-germ1.envt %>%
  group_by(pop) %>%
  summarise(mean.dorm1=mean(fit))

mean.disp<-na.omit(dispersal.data)%>%
  group_by(pop) %>%
  summarise(mean.disp=mean(fit))

disp.dorm<-merge(mean.dorm1,mean.disp,by="pop")
disp.dorm<-merge(disp.dorm,envt,by.x="pop",by.y="site")


lm.dorm_disp<-lmodel2(mean.dorm1~mean.disp,data=disp.dorm, nperm=999)
plot(lm.dorm_disp, "MA",pch=21,cex=1.5,bg="cornflowerblue",xlab="Avg. dispersal ability", ylab="Avg. dormancy rate",main="",las=1,col="orange")
points(y=disp.dorm$mean.dorm1,x=disp.dorm$mean.disp,pch=21,cex=1.5,bg="cornflowerblue")

print(lm.dorm_disp)

#ggsave(height=5,width=5,"./figures/covariance.pdf")
