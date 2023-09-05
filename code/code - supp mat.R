#Important: run main script first since this script uses data frames copiled there

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

###################
#climate over time
###################

df<-read.csv("./data/climate.csv")
head(df)
df$time<-1:nrow(df)


dte_formatter <- function(x) { 
  mth <- substr(format(x, "%b"),1,1) 
  mth 
} 

df$date<-as.Date(df$date, "%Y-%d-%m")

#get annual totals
ggplot(data=df,aes(x=date, y=prec_tot)) + 
  geom_bar(aes(x=date, y=annual),stat="identity", fill="blue", alpha=0.1) +
  geom_point(color="blue") + geom_line(color="blue") + theme_classic() +
  geom_line(aes(y = temp_ave/5),color="red") +
  geom_point(aes(y = temp_ave/5),color="red") +
  scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Temperature °C")) +
  scale_x_date(date_breaks="1 months", date_labels = dte_formatter(df$date)) +
  labs(x="Month",y="Precipitation (mm)", color="blue") +
  theme(axis.title.y = element_text(color = "blue", size = 13),
        axis.title.y.right = element_text(color = "red", size = 13))

ggsave("./figures/climate.pdf")


##############################################################################################
#figure S1, Assumption 1:  Vulpia does best in more productive sites, but only in the absence of competition
##############################################################################################

#test: is fitness of Vulpia higher in productive sites, in the absence of competition only
#data used: transplant fitness + NDVI per site, of seeds that germinated

#data prep
df<-read.csv("./data/cleaned_Fitnessdata.csv")
df$seed<-as.numeric(as.character(df$seed))
df$grid<-as.factor(df$grid)

df$treatment<-ifelse(df$treatment=="A","no","yes")

env<-read.csv("./data/ndvi_szojka_buffer25m2.csv")
head(df); head(env)
df.env<-merge(df,env,by.x="grid",by.y="ï..grid"); head(df.env)

df.seed<-subset(df.env,status=="normal"|status=="no fruit")


#super overdispersed
lm1<-glmmTMB(seed~poly(mean_ndvi,2)*treatment+(1|grid),data=df.seed,family="poisson")
lm2<-glmmTMB(seed~poly(mean_ndvi,2)*treatment+(1|grid),data=df.seed,family="nbinom1")
lm3<-glmmTMB(seed~poly(mean_ndvi,2)*treatment+(1|grid),data=df.seed,family="nbinom2")

lm4<-glmmTMB(seed~poly(mean_ndvi,2)*treatment+(1|grid),data=df.seed,family="poisson", ziformula = ~1)
lm5<-glmmTMB(seed~poly(mean_ndvi,2)*treatment+(1|grid),data=df.seed,family="nbinom1", ziformula = ~1)
lm6<-glmmTMB(seed~poly(mean_ndvi,2)*treatment+(1|grid),data=df.seed,family="nbinom2", ziformula = ~1)

lm7<-glmmTMB(seed~mean_ndvi*treatment+(1|grid),data=df.seed,family="poisson")
lm8<-glmmTMB(seed~mean_ndvi*treatment+(1|grid),data=df.seed,family="nbinom1")
lm9<-glmmTMB(seed~mean_ndvi*treatment+(1|grid),data=df.seed,family="nbinom2")

lm10<-glmmTMB(seed~mean_ndvi*treatment+(1|grid),data=df.seed,family="poisson", ziformula = ~1)
lm11<-glmmTMB(seed~mean_ndvi*treatment+(1|grid),data=df.seed,family="nbinom1", ziformula = ~1)
lm12<-glmmTMB(seed~mean_ndvi*treatment+(1|grid),data=df.seed,family="nbinom2", ziformula = ~1)


#no zero inflation needed
anova(lm1,lm2,lm3,lm4,lm5,lm6,lm7,lm8,lm9,lm10,lm11,lm12) 

#lm9 = lowest AIC, DHARMa looks good

simulationOutput<-simulateResiduals(lm9)
plot(simulationOutput)
testZeroInflation(simulationOutput)


Anova(lm9,type=3) #model with lowest AIC
summary(lm9)
visreg(lm9,xvar="mean_ndvi",by="treatment",scale="response", overlay=TRUE)

#Answer: yes, negligible effect of productivity with competition, huge effect without it
#Answer: should expect sites with lower mean NDVI to be more variable due to variability in competition

#testing by treatment: with vs without competitors
lm.yes<-glmmTMB(seed~mean_ndvi+(1|grid),data=subset(df.seed,treatment=="yes"),family="nbinom2")
lm.no<-glmmTMB(seed~mean_ndvi+(1|grid),data=subset(df.seed,treatment=="no"),family="nbinom2")

#both increase with prod, just A (abiotic only) way more 
Anova(lm.yes)
Anova(lm.no)

visreg(fit=lm.yes, xvar="mean_ndvi", scale="response")
visreg(fit=lm.no, xvar="mean_ndvi", scale="response")

#####plotting

vis<-ggpredict(lm9, terms = c("mean_ndvi[all]","treatment"), type="fe")

ggplot(vis, aes(x = x, y = predicted, group = group, fill = group, color = group)) + 
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha = 0.2, color=NA) +
  geom_line(size=1) +
  labs(x = "Site productivity (NDVI)", y = "Seeds produced per individual") +
  theme_classic()  + 
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA),legend.position = c(0.15,0.75), strip.text.x = element_text(size = 12, face = "bold")) +
  labs(color = "Competitors", fill = "Competitors") 

ggsave("./figures/assumption 1 - transplants.pdf", width = 5, height = 5)


##############################################################################################
#Figure S3, Assumption 2:  low productivity sites are more variable
##############################################################################################

env<-read.csv("./data/env.means.all.csv")
df<-read.csv("./data/percent green.csv"); head(df)
df$area<-df$X.area/100
df$Year<-as.factor(df$Year)
df$Plot.<-as.factor(df$Plot.)

df.env<-merge(df, env, by.x = "Plot.", by.y = "site")

lm1<-glmmTMB(area ~ productivity * Year + (1|Plot.), family=beta_family(link="logit"), data=df.env)
lm2<-glmmTMB(area ~ productivity * Year + (1|Plot.), family=beta_family(link="logit"), data=subset(df.env,Plot.!=5))


simulationOutput<-simulateResiduals(lm1)
plot(simulationOutput)

simulationOutput<-simulateResiduals(lm2)
plot(simulationOutput)


Anova(lm1, type=2) #interaction not sig so type 2
Anova(lm2, type=3)


vis<-ggpredict(lm1, terms = c("productivity","Year"), type="fe")

figS3a<-ggplot(vis, aes(x = x, y = predicted*100, group = group, fill = group, color = group)) + 
  geom_ribbon(aes(ymax=conf.high*100, ymin=conf.low*100), alpha = 0.2, color=NA) +
  geom_line(size=1) +
  labs(x = "Site productivity (NDVI)", y = "Percent greenness in plot") +
  theme_classic()  + 
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA),legend.position = c(0.2,0.75), strip.text.x = element_text(size = 12, face = "bold")) +
  labs(color = "Competitors", fill = "Competitors") 

vis<-ggpredict(lm2, terms = c("productivity","Year"), type="fe")

figS3b<-ggplot(vis, aes(x = x, y = predicted*100, group = group, fill = group, color = group)) + 
  geom_ribbon(aes(ymax=conf.high*100, ymin=conf.low*100), alpha = 0.2, color=NA) +
  geom_line(size=1) +
  labs(x = "Site productivity (NDVI)", y = "Percent greenness in plot") +
  theme_classic()  + 
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA),legend.position = "none", strip.text.x = element_text(size = 12, face = "bold")) +
  labs(color = "Competitors", fill = "Competitors") 


plot_grid(figS3a,figS3b, labels = c("A","B"))
ggsave("./figures/prod vs time.png", width = 8, height = 4)


##############################################################################################
#Figure S4, spatial autocorrelation
##############################################################################################

df1<-read.csv("./data/raw env data.csv")
df2<-read.csv("./data/env.means.all, lat lon.csv")

df<-merge(df1,df2[,c(1,15)])

head(df)

require(ggplot2); require(cowplot)

#order by productivity
df$site = factor(df$site, levels=c(4,28,2,19,20,1,26,27,5,29,16,11,18,10,14,8,13,17,25,6,12,24,7,15,3,9,23,30,22,21)) 

#cut out sites with no population
pops<-c(2,3,4,5,6,9,10,11,13,15,16,17,18,19,20,23,25,26,27,29)

df<-subset(df, site %in% pops)

df$soil.depth<-as.numeric(df$soil.depth)


#soil depth
figa<-ggplot(data=df, aes(x=site,y=soil.depth, color=productivity)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5),linetype="dashed") +
  labs(x="",y="Soil depth (cm)")

#soil moisture
figb<-ggplot(data=df, aes(x=site,y=sm1, color=productivity)) +
  geom_text(aes(label=plot),position=position_jitterdodge()) +
  geom_text(aes(x=site,y=sm2,label=plot),position=position_jitterdodge()) +
  geom_text(aes(x=site,y=sm3,label=plot),position=position_jitterdodge()) +
  theme_classic() +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5),linetype="dashed") +
  labs(x="Site",y="Soil moisture (%)")

plot_grid(figa,figb,labels="auto",nrow=2)

ggsave("./figures/autocorrelation of env.png", width = 8, height = 4)


################
#plot distances
################

geog<-read.csv("./data/env.means.all, lat lon.csv")

#germ1.envt is generated in main code file
merge(geog, germ1.envt)

geog<-geog %>%
  filter(site %in% germ1.envt$pop) %>%
  select(site,lat,long)

ggplot(data=geog,aes(x=long,y=lat)) +
  geom_text(aes(label=site))

expand.grid(geog$site,geog$site)

sites<-combn(geog$site, m=2)

library(geosphere)
distm(c(geog$long[1], geog$lat[1]), c(geog$long[2], geog$lat[2]), fun = distHaversine)


distances<-c()

for(i in 1:ncol(sites)) {
  
  site1<-as.numeric(sites[1,i])
  site2<-as.numeric(sites[2,i])
  
  tmp1<-subset(geog,site==site1)
  tmp2<-subset(geog,site==site2)
  
  tmp<-distm(c(tmp1$long, tmp1$lat), c(tmp2$long, tmp2$lat), fun = distHaversine)
  
  distances<-c(distances, tmp)
  
}

mean(distances)
min(distances)
max(distances)


#################################################
# Figure S5 - testing how dormancy varies among populations when the trial is performed under suitable conditions
#################################################

germ1.envt<-merge(germ1,envt,by.x="pop",by.y="site")
germ1.envt$n<-7

germ1.envt<-na.omit(germ1.envt)

lm.dorm_prod<-glmmTMB(prop.dorm~productivity+(1|pop),data=germ1.envt,weights=n,family="binomial")
Anova(lm.dorm_prod)
germ1.envt$fit<-fitted.values(lm.dorm_prod)


germ2.envt<-merge(germ2,envt,by.x="pop",by.y="site")
germ2.envt$n<-7

lm<-glmmTMB(prop.dorm~productivity+(1|pop),data=germ2.envt,weights=n,family="binomial")
Anova(lm)
germ2.envt$fit<-fitted.values(lm)

#plotting

ggplot(germ1.envt,aes(y=prop.dorm,x=productivity,group=pop)) + geom_boxplot(fill="lightgrey",color="grey") +  theme_classic()+
  geom_point(aes(x=productivity,y=fit),color="black",size=2) + labs(x = "Site productivity", y = "Dormancy rate") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) + 
  geom_boxplot(data=germ2.envt,fill="pink",color="pink") + geom_point(data=germ2.envt,aes(x=productivity,y=fit),color="red",size=2) 

ggsave(height=5,width=5,"./figures/fig S7 - dorm rate, both trials.pdf")

