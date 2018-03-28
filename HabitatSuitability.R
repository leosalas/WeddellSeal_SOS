# TODO: explore effects of di_drygals, dist_eit, dist_250
# 
# Author: lsalas
###############################################################################


#load the data
library(XLConnect)
library(ggplot2)

#we will change this path to connect directly to GitHub when we are done, so any user can just copy-paste the code in R and the data will
#be read from the GitHub repo
basepth<-"//prbo.org/Data/Home/Petaluma/lsalas/Documents/lsalas/Antarctica/SealsFromSpace/HabitatSuitability"

data<-try(readWorksheetFromFile(paste(basepth,"5kmgrid_RossSeaAll_geodesic_20180212.xlsx",sep="/"),sheet="5kmgrid_RS2011_geodesic"))

if(inherits(data,"try-error"))stop("Could not read the data. Please check the path, name of file and name of sheet")

data$slope<-data$Max_GRID_1-data$Min_GRID_1
minedge<-min(subset(data,di_edge>0)$di_edge);data$lgdi_edge<-log((data$di_edge+minedge)/1000)		#add min distance to avoid log(0)
#minshlf<-min(subset(data,dist_shlf>0)$dist_shlf);data$lgdi_shlf<-log((data$dist_shlf+minshlf)/1000)	#add min distance
mincont<-min(subset(data,dist_cont>0)$dist_cont);data$lgdi_cont<-log((data$dist_cont+mincont)/1000)		#add min distance to avoid log(0)
minglac<-min(subset(data,dist_glac>0)$dist_glac);data$lgdi_glac<-log((data$dist_glac+minglac)/1000)	#add min distance
minisld<-min(subset(data,dist_isld>0)$dist_isld);data$lgdi_isld<-log((data$dist_isld+minisld)/1000)		#add min distance to avoid log(0)
mindadpe<-min(subset(data,dist_adpe>0)$dist_adpe);data$lgdi_adpe<-log((data$dist_adpe+mindadpe)/1000)	#add min distance
mindempe<-min(subset(data,dist_empe>0)$dist_empe);data$lgdi_empe<-log((data$dist_empe+mindempe)/1000)		#add min distance to avoid log(0)
mincadpe<-min(subset(data,near_adpe>0)$near_adpe);data$lgnear_adpe<-log((data$near_adpe+mincadpe))	#add min distance
mincempe<-min(subset(data,near_empe>0)$near_empe);data$lgnear_empe<-log((data$near_empe+mincempe))	#add min distance
#data$lgdi_drygals<-log((data$di_drygals/1000)+1);data$lgdi_eit<-log((data$dist_eit/1000)+1)
data$lgdi_250<-log((data$dist_250/1000)+1);data$lgdi_deep<-log((data$dist_deep/100)+1)
minicew<-min(subset(data,icewidth>0)$icewidth);data$lgicewidth<-log((data$icewidth+minicew))
#data$bindi_shlf<-ifelse(data$lgdi_shlf<5.849,0,1) #using the median to create the dychotomous dummy
data$ratio_edge<-ifelse(data$icewidth==0,0,
		ifelse(data$di_edge==0,minedge/data$icewidth,data$di_edge/data$icewidth))
minratio<-min(subset(data,ratio_edge>0)$ratio_edge)
data$lgratio_edge<-log(data$ratio_edge+minratio)

#create a variable lgdi_shore<-min of lgdi_cont or lgdi_isld
data$lgdi_shore<-ifelse(data$lgdi_cont<data$lgdi_isld,data$lgdi_cont,data$lgdi_isld)

#let's plot the values and see how they cluster - this may show justification for quadratic forms
covars<-c("slope","lgdi_edge","lgdi_cont","lgdi_glac","lgdi_isld","lgdi_adpe","lgdi_empe",
		"lgnear_adpe","lgnear_empe","Min_GRID_1","SD_GRIDC_1","lgdi_250","lgdi_deep","lgicewidth","lgratio_edge")	#"lgdi_shlf", "lgdi_shore",

#Add here all the covariates described below!!
idvar<-c("Id","sealspres","tomnodaoi","crack","Year")	#,"highres"
df<-data[,c(idvar,covars)]
dfr<-reshape(data=df,idvar=idvar,varying=list(6:19),direction="long",times=names(df)[6:19],timevar="Parameter",v.names="Value",new.row.names=as.character(1:(nrow(df)*15)))
#ggplot(data=dfr,aes(x=Value)) + geom_histogram() + facet_wrap(~Parameter,ncol=4,scales="free")

#let's fit the first model - let's call this the "full" model
allcovars<-c("tomnodaoi","crack","Year","mlsearch",covars)
intcovars<-c("lgdi_cont*drygalski","lgdi_empe*drygalski","lgdi_edge*drygalski","drygalski")	#"lgdi_shlf*lgdi_edge","bindi_shlf*lgdi_edge",,"lgdi_shlf:lgdi_edge*drygalski"
#fmlfull<-as.formula(paste("sealspres~",paste(allcovars,collapse="+"),"+I(",paste(covars,collapse="^2)+I("),"^2)+",paste(intcovars,collapse="+"),sep=""))
fmlfull<-as.formula(paste("sealspres~",paste(allcovars,collapse="+"),"+I(",paste(covars,collapse="^2)+I("),"^2)",sep=""))

fullmodel<-glm(formula=fmlfull,data=data,family="binomial")
summary(fullmodel)

# Find top model via stepwise, using BIC as criterion to avoid AIC's over-parametrization issues
k<-log(nrow(data))
stepm<-step(fullmodel,k=k)	#using BIC as the selection criterion 
summary(stepm)

#Not an ideal fit...

#Our top model is really:
topfml<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgdi_edge + lgdi_glac + lgdi_deep + 
    lgicewidth + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgdi_empe^2) + I(lgdi_deep^2)")
topmdl<-glm(formula=topfml,data=data,family="binomial");summary(topmdl)

topfml_l<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgdi_edge + lgdi_glac + lgdi_deep + lgdi_cont + lgdi_adpe +
    lgnear_empe + lgicewidth + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_l<-glm(formula=topfml_l,data=data,family="binomial");summary(topmdl_l)

#not worth adding the linear effects to those quadratics that do not heve them

# Remove ERS and redo
datawrs<-subset(data,Region=="WRS")
fullmodelwrs<-glm(formula=fmlfull,data=datawrs,family="binomial")
summary(fullmodelwrs)
k<-log(nrow(datawrs))
stepmwrs<-step(fullmodelwrs,k=k)	#using BIC as the selection criterion 
summary(stepmwrs)

topmdl_wrs<-glm(formula=topfml,data=data,family="binomial");summary(topmdl_wrs)

#we tried using shore on both the overal and the WRS models - not best.

#2) rearrange the factor for AllRegion to use WRSS as the ref and do the shotgun below
data$AllRegion<-ifelse(data$Region=="ERS","ERS",ifelse(data$drygalski==0,"WRSS","WRSN"))
data$regionOrder<-ifelse(data$Region=="ERS",2,ifelse(data$drygalski==0,1,3))
data$AllRegion<-as.factor(data$AllRegion)
data$AllRegion<-reorder(data$AllRegion,data$regionOrder)

#add some interactions:
topfml_ra<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + AllRegion*lgdi_edge + lgdi_glac + AllRegion*lgdi_deep + 
				lgicewidth + AllRegion*lgdi_adpe + AllRegion*lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_ra<-glm(formula=topfml_ra,data=data,family="binomial");summary(topmdl_ra)

#explore using lgratio_edge:
topfml_ra<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + AllRegion*lgratio_edge + lgdi_glac + AllRegion*lgdi_deep + 
				lgicewidth + AllRegion*lgdi_adpe + AllRegion*lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_ra<-glm(formula=topfml_ra,data=data,family="binomial");summary(topmdl_ra)

#using logratio but no interactions:
topfml_nint<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgratio_edge + lgdi_glac + lgdi_deep + 
				lgicewidth + lgdi_adpe + lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_nint<-glm(formula=topfml_nint,data=data,family="binomial");summary(topmdl_nint)

##PLOTS: 
# Partial dependence of individual effects

# will need these functions...
makePredictions<-function(df,mdlfit,varnam,othervars=NULL){
	preddf<-predict(mdlfit,df,se=T)
	df$lgpred<-preddf$fit; df$sepred<-preddf$se.fit
	df$lgymin<-df$lgpred-(1.96*df$sepred)
	df$lgymax<-df$lgpred+(1.96*df$sepred)
	df$predicted<-exp(df$lgpred)/(1+exp(df$lgpred))
	df$ymin<-exp(df$lgymin)/(1+exp(df$lgymin))
	df$ymax<-exp(df$lgymax)/(1+exp(df$lgymax))
	if(!is.null(othervars)){
		df<-df[,c(varnam,"predicted","ymin","ymax",othervars)]
	}else{
		df<-df[,c(varnam,"predicted","ymin","ymax")]
	}
	return(df)
}

makedf<-function(reg,var,data){
	#based on the region, make the means and sequences
	# will need these data...
	listseq<-list();listmn<-list()
	mn_lgratio_edge<-mean(subset(data,AllRegion %in% reg)$lgratio_edge,na.rm=T); listmn[["lgratio_edge"]]<-mn_lgratio_edge
	seq_lgratio_edge<-seq(from=min(data$lgratio_edge,na.rm=T), to=max(data$lgratio_edge,na.rm=T),length.out=100);listseq[["lgratio_edge"]]<-seq_lgratio_edge
	mn_lgdi_glac<-mean(subset(data,AllRegion %in% reg)$lgdi_glac,na.rm=T); listmn[["lgdi_glac"]]<-mn_lgdi_glac
	seq_lgdi_glac<-seq(from=min(data$lgdi_glac,na.rm=T), to=max(data$lgdi_glac,na.rm=T),length.out=100);listseq[["lgdi_glac"]]<-seq_lgdi_glac
	mn_lgdi_deep<-mean(subset(data,AllRegion %in% reg)$lgdi_deep,na.rm=T); listmn[["lgdi_deep"]]<-mn_lgdi_deep
	seq_lgdi_deep<-seq(from=min(data$lgdi_deep,na.rm=T), to=max(data$lgdi_deep,na.rm=T),length.out=100);listseq[["lgdi_deep"]]<-seq_lgdi_deep
	mn_lgicewidth<-mean(subset(data,AllRegion %in% reg)$lgicewidth,na.rm=T); listmn[["lgicewidth"]]<-mn_lgicewidth
	seq_lgicewidth<-seq(from=min(data$lgicewidth,na.rm=T), to=max(data$lgicewidth,na.rm=T),length.out=100);listseq[["lgicewidth"]]<-seq_lgicewidth
	mn_lgdi_adpe<-mean(subset(data,AllRegion %in% reg)$lgdi_adpe,na.rm=T); listmn[["lgdi_adpe"]]<-mn_lgdi_adpe
	seq_lgdi_adpe<-seq(from=min(data$lgdi_adpe,na.rm=T), to=max(data$lgdi_adpe,na.rm=T),length.out=100);listseq[["lgdi_adpe"]]<-seq_lgdi_adpe
	mn_lgdi_empe<-mean(subset(data,AllRegion %in% reg)$lgdi_empe,na.rm=T); listmn[["lgdi_empe"]]<-mn_lgdi_empe
	seq_lgdi_empe<-seq(from=min(data$lgdi_empe,na.rm=T), to=max(data$lgdi_empe,na.rm=T),length.out=100);listseq[["lgdi_empe"]]<-seq_lgdi_empe
	mn_lgdi_cont<-mean(subset(data,AllRegion %in% reg)$lgdi_cont,na.rm=T); listmn[["lgdi_cont"]]<-mn_lgdi_cont
	seq_lgdi_cont<-seq(from=min(data$lgdi_cont,na.rm=T), to=max(data$lgdi_cont,na.rm=T),length.out=100);listseq[["lgdi_cont"]]<-seq_lgdi_cont
	mn_lgnear_empe<-mean(subset(data,AllRegion %in% reg)$lgnear_empe,na.rm=T); listmn[["lgnear_empe"]]<-mn_lgnear_empe
	seq_lgnear_empe<-seq(from=min(data$lgnear_empe,na.rm=T), to=max(data$lgnear_empe,na.rm=T),length.out=100);listseq[["lgnear_empe"]]<-seq_lgnear_empe
	
	newdf<-data.frame(tomnodaoi=rep(1,100),crack=rep(1,100),mlsearch=rep(1,100),Year=rep(2011,100))
	newdf[,var]<-listseq[[var]]
	for(nn in names(listmn)){
		if(nn != var){
			newdf[,nn]<-listmn[[nn]]
		}
	}
	if(NROW(reg)==1){
		newdf$AllRegion<-reg
	}else{
		newdf$AllRegion<-"All regions"
	}
	
	return(newdf)
}

#DO this first overall then by region

#distance to deep
pdf<-makedf(reg="WRSN",var="lgdi_deep",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_deep");plotdfn$Region<-"WRSN"
pdf<-makedf(reg="WRSS",var="lgdi_deep",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_deep");plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgdi_deep",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_deep");plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgdi_deep");nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
p<-ggplot(data=plotdf,aes(x=lgdi_deep,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Distance to deep waters)")
dev.new();print(p)
names(plotdf)<-gsub("lgdi_deep","CovariateValue",names(plotdf))
plotdf$CovarName<-"Log(Dist. to deep waters)"
regionsdf<-plotdf
names(nreg)<-gsub("lgdi_deep","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to deep waters)"
regionsdf<-plotdf
nregdf<-nreg

#lgratio_edge
pdf<-makedf(reg="WRSN",var="lgratio_edge",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgratio_edge");plotdfn$Region<-"WRSN"
pdf<-makedf(reg="WRSS",var="lgratio_edge",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgratio_edge");plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgratio_edge",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgratio_edge");plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgratio_edge");nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
p<-ggplot(data=plotdf,aes(x=lgratio_edge,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Dist. to ice edge/ice width)")
dev.new();print(p)
names(plotdf)<-gsub("lgratio_edge","CovariateValue",names(plotdf))
plotdf$CovarName<-"Log(Dist. to ice edge/ice width)"
regionsdf<-rbind(regionsdf,plotdf)
names(nreg)<-gsub("lgratio_edge","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to ice edge/ice width)"
nregdf<-rbind(nregdf,nreg)

#distance to emperor penguins
pdf<-makedf(reg="WRSN",var="lgdi_empe",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_empe");plotdfn$Region<-"WRSN"
pdf<-makedf(reg="WRSS",var="lgdi_empe",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_empe");plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgdi_empe",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_empe");plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgdi_empe");nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
p<-ggplot(data=plotdf,aes(x=lgdi_empe,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Distance to EMPE)")
dev.new();print(p)
names(plotdf)<-gsub("lgdi_empe","CovariateValue",names(plotdf))
plotdf$CovarName<-"Log(Dist. to EMPE)"
regionsdf<-rbind(regionsdf,plotdf)
names(nreg)<-gsub("lgdi_empe","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to EMPE)"
nregdf<-rbind(nregdf,nreg)

#distance to Adelie penguins
pdf<-makedf(reg="WRSN",var="lgdi_adpe",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_adpe");plotdfn$Region<-"WRSN"
pdf<-makedf(reg="WRSS",var="lgdi_adpe",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_adpe");plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgdi_adpe",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_adpe");plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgdi_adpe");nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
p<-ggplot(data=plotdf,aes(x=lgdi_adpe,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Distance to ADPE)")
dev.new();print(p)
names(plotdf)<-gsub("lgdi_adpe","CovariateValue",names(plotdf))
plotdf$CovarName<-"Log(Dist. to ADPE)"
regionsdf<-rbind(regionsdf,plotdf)
names(nreg)<-gsub("lgdi_adpe","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to ADPE)"
nregdf<-rbind(nregdf,nreg)
nregdf$RegionName<-"All regions"

regionsdf$RegionName<-ifelse(regionsdf$Region=="ERS","Eastern Ross Sea",
		ifelse(regionsdf$Region=="WRSS","Western Ross Sea South","Western Ross Sea North"))

p<-ggplot(data=regionsdf,aes(x=CovariateValue,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_grid(RegionName~CovarName,scales="free") +
		labs(y="Probability of seal presence", x="Covariate Value")

alldata<-rbind(regionsdf,nregdf)
alldata$CovarLetter<-ifelse(alldata$CovarName=="Log(Dist. to ADPE)","A",
		ifelse(alldata$CovarName=="Log(Dist. to EMPE)","B",
				ifelse(alldata$CovarName=="Log(Dist. to deep waters)","C","D")))
p<-ggplot(data=alldata,aes(x=CovariateValue,y=predicted)) + 
		geom_line(color="black",size=2) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_grid(RegionName~CovarLetter,scales="free") +
		labs(y="Probability of seal presence", x="Covariate Value") +
		theme(axis.text.x=element_text(size=24, color="black")) + theme(axis.title.x=element_text(size=28,face="bold")) +
		theme(axis.text.y=element_text(size=24, color="black")) + theme(axis.title.y=element_text(size=28,face="bold")) +
		theme(strip.text.x=element_text(size=24)) + theme(strip.text.y=element_text(size=24))

jpeg(filename=paste(basepth,"/paper/plots/Plasticity.jpg",sep=""),width=1500,height=1400,quality=100)
print(p)
dev.off()


#then the non-interacting effects: lgicewidth, lgdi_glac, lgdi_cont, lgnear_empe
#lgicewidth
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgicewidth",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgicewidth");plotdfe$Region<-"All Regions"

p<-ggplot(data=plotdf,aes(x=lgicewidth,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability seal presence", x="Log(Fast ice width)") +
		theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.title.x=element_text(size=22,face="bold")) +
		theme(axis.text.y=element_text(size=20, color="black")) + theme(axis.title.y=element_text(size=24,face="bold"))

jpeg(filename=paste(basepth,"/paper/plots/FastIceWidth.jpg",sep=""),width=400,height=400,quality=100)
print(p)
dev.off()


#lgdi_glac
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgdi_glac",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_glac");plotdfe$Region<-"All Regions"

p<-ggplot(data=plotdf,aes(x=lgdi_glac,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability of seal presence", x="Log(Distance to glaciers)") +
		theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.title.x=element_text(size=22,face="bold")) +
		theme(axis.text.y=element_text(size=20, color="black")) + theme(axis.title.y=element_text(size=24,face="bold"))

jpeg(filename=paste(basepth,"/paper/plots/DistToGlaciers.jpg",sep=""),width=400,height=400,quality=100)
print(p)
dev.off()


#lgdi_cont
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgdi_cont",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_cont");plotdfe$Region<-"All Regions"

p<-ggplot(data=plotdf,aes(x=lgdi_cont,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability of seal presence", x="Log(Distance to shore)") +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold"))

jpeg(filename=paste(basepth,"/paper/plots/DistToContinent.jpg",sep=""),width=400,height=400,quality=100)
print(p)
dev.off()


#lgnear_empe
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgnear_empe",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgnear_empe");plotdfe$Region<-"All Regions"

p<-ggplot(data=plotdf,aes(x=lgnear_empe,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability of seal presence", x="Log(Size of EMPE colony)") +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold"))

jpeg(filename=paste(basepth,"/paper/plots/AbundNearestEMPE.jpg",sep=""),width=400,height=400,quality=100)
print(p)
dev.off()

###Let's examine the relationship with distsnce to deep waters for close distances
data_shll<-subset(data,lgdi_deep<=2)
topmdl_shll<-glm(formula=topfml_ra,data=data_shll,family="binomial");summary(topmdl_shll)
topfml_lin<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + AllRegion*lgratio_edge + lgdi_glac + AllRegion*lgdi_deep + 
				lgicewidth + AllRegion*lgdi_adpe + AllRegion*lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2)")
topmdl_shll_lin<-glm(formula=topfml_lin,data=data_shll,family="binomial");summary(topmdl_shll_lin)
anova(topmdl_shll_lin, topmdl_shll)
1-pchisq(3.6115,1)  #Almost significant! but only because the quadrature is NOT what we thought.

topfml_nint<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgratio_edge + lgdi_glac + lgdi_deep + 
				lgicewidth + lgdi_adpe + lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_nint<-glm(formula=topfml_nint,data=data_shll,family="binomial");summary(topmdl_nint)
nreg<-makePredictions(df=data_shll,mdlfit=topmdl_nint,varnam="lgdi_deep");nreg$Region<-"ALL"

p<-ggplot(nreg,aes(x=lgdi_deep,y=predicted)) + geom_point() +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Prob. seal presence", x="Log(Dist. to deep waters)") +
		theme(axis.text.x=element_text(size=22)) + theme(axis.title.x=element_text(size=24)) +
		theme(axis.text.y=element_text(size=22)) + theme(axis.title.y=element_text(size=26))

