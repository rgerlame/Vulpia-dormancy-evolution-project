#setwd("C:/Users/16047/Dropbox/2 - Manuscripts in prep/3 - in revision - Kelley - evolution of dispersal and dormancy/AmNat resubmission/data")


#add in variability to each in fitness due to climate + competition, spatial decay in competition
#currently not working: # individuals not going down as 0 are produced
#death or burial rate
#

library(tidyr)
library(dplyr)
library(ggplot2)

###################
####model
###################

df.all<-data.frame(dormancy.rate=vector(), time=vector(), env=vector(), n=vector(), n_dorm=vector())

#initial conditions
df<-data.frame(env=vector(), dorm.rate=vector(), cells=vector())

env<-seq(from=0, to = 1, by = 0.1)
d<-seq(from=0, to=0.9,by=0.1)

#for each environment j (0 = bad) and each dormancy level k
#create a 100 x 100 grid; 

##saved data frame to be called if model crashes mid-run
#df.all<-read.csv("C:/Users/Rachel/Documents/plastic dormancy, 2 seeds.csv")

for(j in 1:length(env)) {
  
  for(k in 1:length(d)) {
    
    inds           <- array(data = 0, dim = c(10000, 9))
    colnames(inds) <- c("env","dormancy.rate", "dispersal.rate", "x_loc", "y_loc","dorm","birth","loc", "ind");
    
    rows <- sample(nrow(expand.grid(1:100,1:100)))
    grid <- grid[rows, ]
    
    #initial placement
    #all in unique random cell
    inds[, 1]<-rbinom(n=10000, prob=env[j], size=1)
    inds[, 2]<-rep(c(d[k]),times=10000)
    inds[, 4]<-grid[,1]
    inds[, 5]<-grid[,2]
    
    inds.orig<-as.data.frame(inds)
    inds.orig$loc<-paste(inds.orig$x_loc,inds.orig$y_loc,sep=",")
    
    #change to 100
    for(t in 1:100){ tryCatch({
      
      
      inds<-as.data.frame(inds)
      
      inds<-inds %>% relocate(ind, .after = loc)
      
      colnames(inds)[1:9] <- c("env","dormancy.rate", "dispersal.rate", "x_loc", "y_loc","dorm","birth","loc","ind"); 
      
      #inds[, 6] <- rbinom(n=nrow(inds),prob=inds[, 2], size=1) #non-plastic version
      inds[, 6] <- ifelse(inds[, 1]==0,rbinom(n=nrow(inds),prob=inds[, 2], size=1),0) #plastic version
      
      inds<-as.data.frame(inds)
      
      inds$loc<-paste(inds$x_loc,inds$y_loc,sep=",")
      inds<-inds %>% relocate(ind, .after = loc)
      
      inds<-as.data.frame(inds)
      
      #separates the dormant and non-dormant outcomes
      #if more than 1 non-dormant, select 1
      x1 <- inds %>% filter(dorm==1)
      x2 <- inds %>% filter(dorm==0) %>% group_by(loc) %>% slice_sample(n = 1)
      
      x1<-as.data.frame(x1)
      x2<-as.data.frame(x2)
      
      inds <- rbind(x1,x2) 
      
      print(table(inds$dormancy.rate))
      
      tmp1<-data.frame(dormancy.rate=c(dormancy.rate=seq(from=0,to=0.9,by=0.1)), time=t, env=env[j])
      tmp2<-inds %>% group_by(dormancy.rate) %>%
        summarise(n=n(), n_dorm = sum(dorm == 1))
      
      tmp3<-left_join(tmp1,tmp2)
      
      ###WHAT IS THIS???? I think it's if I had 2 species comparison
      #tmp3<-subset(tmp3,dormancy.rate==0.0|dormancy.rate==d[k])
      tmp3<-subset(tmp3,dormancy.rate==d[k])
      
      df.all<-rbind(df.all,tmp3)
      
      df.all[is.na(df.all)] <- 0
      
      print(tail(df.all))
      
      for( i in 1: dim(inds)[1]){
        
        inds[i, 9] <- i
        
        #with dormancy, 1 is carried forward regardless of environment
        #without dormancy, if the environment is bad, 0 are produced
        #without dormancy, if the environment is good, 2 are produced        
        inds[i, 7] <- ifelse(inds[i, 6]==1,1,inds[i, 7]) 
        inds[i, 7] <- ifelse(inds[i, 6]!=1 & inds[i, 1]==0, 0, inds[i, 7])
        inds[i, 7] <- ifelse(inds[i, 6]!=1 & inds[i, 1]==1, 2, inds[i, 7])
        
        #each of those offspring take on the dormancy of the parent
        #make them disperse
        
      }   
      
      inds<-as.data.frame(inds)
      
      x<-subset(inds,birth>0)
      
      #error here because 0 rows

      #expands each birth to its own row, which can then be dispersed, with same dispersal rate as parent
      if (nrow(x)>0) {
        x<-x %>% group_by(ind, env, dormancy.rate, dispersal.rate, x_loc, y_loc, dorm) %>%
          complete(birth = full_seq(1:birth, 1))
      }  else {
        x<-x
      }
      
      
      #"disperses" each birth by re-assigning location
      #torus, uniform from all over landscape 
      x$x_loc<-sample(x=1:100, replace=TRUE, size=dim(x)[1])
      x$y_loc<-sample(x=1:100, replace=TRUE, size=dim(x)[1])

      x$loc<-paste(x$x_loc,x$y_loc,sep=",")
      
      #randomly select 1 individual per location
      
      #####IS THIS CORRECT?
      
      if (nrow(x)>0) {
        x$ind <- 1:nrow(x)
      }  else{
        inds['ind'] <- vector()
      }
      
      x<- x %>% relocate(ind, .after = loc)
      
      inds<-x
      inds<-as.data.frame(inds)
      
      inds<- inds %>% relocate(ind, .after = loc)
      
      inds<-merge(inds,inds.orig[c(1,8)],by="loc")
      inds<-inds[,c(1,10,3:9)]
      inds<-inds %>% relocate(loc, .after = ind)
      colnames(inds)[1]<-"env"
      head(inds)
      
      
    }, error=function(e){}) 
    }}
  write.csv(df.all, "plastic dormancy, 2 seeds.csv")
  }
    
  

#without plasticity, delays the inevitable
ggplot(data=df.all,aes(x=time,y=n,color=dormancy.rate)) +
  geom_point() +
  facet_wrap(~env, scales="free")

write.csv(df.all, "plastic dormancy, 2 seeds.csv")

##########
#plotting
###########

df.all<-read.csv("C:/Users/Rachel/Desktop/plastic dormancy, 2 seeds.csv")


df.means<-df.all %>% group_by(dormancy.rate, env) %>%
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
  scale_color_continuous(limits = c(0,0.9), breaks = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  scale_x_continuous(breaks = seq(from=0, to=1, by = 0.1), limits=c(0,1))

ggsave(height=4,width=8,"simulation.pdf")

