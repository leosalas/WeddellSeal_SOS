# TODO: explore effects of di_drygals, dist_eit, dist_250
# 
# Author: lsalas
###############################################################################


#load the data
library(XLConnect)
library(ggplot2)
library(generalhoslem)
library(plyr)

# will need these functions...
makePredictions<-function(df,mdlfit,varnam,othervars=NULL,intsC,reg,conf=95){
	preddf<-predict(mdlfit,df,se=T)
	adjM<-intsC[reg]
	curmn<-mean(preddf$fit);adjC<-curmn-adjM
	df$lgpred<-preddf$fit-adjC; df$sepred<-preddf$se.fit
	if(conf==95){
		df$lgymin<-df$lgpred-(1.96*df$sepred)	
		df$lgymax<-df$lgpred+(1.96*df$sepred)	
	}else if(conf==0){
		df$lgymin<-df$lgpred
		df$lgymax<-df$lgpred
	}else{
		df$lgymin<-df$lgpred-(df$sepred)	
		df$lgymax<-df$lgpred+(df$sepred)	
	}
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
	seq_lgratio_edge<-seq(from=min(subset(data,AllRegion %in% reg)$lgratio_edge,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgratio_edge,na.rm=T),length.out=100);listseq[["lgratio_edge"]]<-seq_lgratio_edge
	mn_lgdi_glac<-mean(subset(data,AllRegion %in% reg)$lgdi_glac,na.rm=T); listmn[["lgdi_glac"]]<-mn_lgdi_glac
	seq_lgdi_glac<-seq(from=min(subset(data,AllRegion %in% reg)$lgdi_glac,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgdi_glac,na.rm=T),length.out=100);listseq[["lgdi_glac"]]<-seq_lgdi_glac
	mn_lgdi_deep<-mean(subset(data,AllRegion %in% reg)$lgdi_deep,na.rm=T); listmn[["lgdi_deep"]]<-mn_lgdi_deep
	seq_lgdi_deep<-seq(from=min(subset(data,AllRegion %in% reg)$lgdi_deep,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgdi_deep,na.rm=T),length.out=100);listseq[["lgdi_deep"]]<-seq_lgdi_deep
	mn_lgicewidth<-mean(subset(data,AllRegion %in% reg)$lgicewidth,na.rm=T); listmn[["lgicewidth"]]<-mn_lgicewidth
	seq_lgicewidth<-seq(from=min(subset(data,AllRegion %in% reg)$lgicewidth,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgicewidth,na.rm=T),length.out=100);listseq[["lgicewidth"]]<-seq_lgicewidth
	mn_lgdi_adpe<-mean(subset(data,AllRegion %in% reg)$lgdi_adpe,na.rm=T); listmn[["lgdi_adpe"]]<-mn_lgdi_adpe
	seq_lgdi_adpe<-seq(from=min(subset(data,AllRegion %in% reg)$lgdi_adpe,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgdi_adpe,na.rm=T),length.out=100);listseq[["lgdi_adpe"]]<-seq_lgdi_adpe
	mn_lgdi_empe<-mean(subset(data,AllRegion %in% reg)$lgdi_empe,na.rm=T); listmn[["lgdi_empe"]]<-mn_lgdi_empe+1.832581 #this is -log(mindempe/1000)
	seq_lgdi_empe<-seq(from=min(subset(data,AllRegion %in% reg)$lgdi_empe,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgdi_empe,na.rm=T),length.out=100);listseq[["lgdi_empe"]]<-seq_lgdi_empe
	mn_lgdi_cont<-mean(subset(data,AllRegion %in% reg)$lgdi_cont,na.rm=T); listmn[["lgdi_cont"]]<-mn_lgdi_cont
	seq_lgdi_cont<-seq(from=min(subset(data,AllRegion %in% reg)$lgdi_cont,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgdi_cont,na.rm=T),length.out=100);listseq[["lgdi_cont"]]<-seq_lgdi_cont
	mn_lgnear_empe<-mean(subset(data,AllRegion %in% reg)$lgnear_empe,na.rm=T); listmn[["lgnear_empe"]]<-mn_lgnear_empe
	seq_lgnear_empe<-seq(from=min(subset(data,AllRegion %in% reg)$lgnear_empe,na.rm=T), to=max(subset(data,AllRegion %in% reg)$lgnear_empe,na.rm=T),length.out=100);listseq[["lgnear_empe"]]<-seq_lgnear_empe
	
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

estExplaniedDeviance<-function(topmdl,vars,data){
	devdf<-ldply(.data=vars, .fun=function(x,vars,data,topmdl){
				redvars<-subset(vars,!vars %in% x);
				redfml<-paste("sealspres ~ ",paste(redvars,collapse=" + "));
				redmdl<-glm(formula=redfml,data=data,family="binomial");
				devexp<-redmdl$deviance-topmdl_ra$deviance;
				tdf<-data.frame(covariate=x,deviance_exp=devexp);
				return(tdf)
			},vars=vars,data=data,topmdl=topmdl)
	return(devdf)
}


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

#let's fit the first model - let's call this the "full" model
allcovars<-c("tomnodaoi","crack","Year","mlsearch",covars)
fmlfull<-as.formula(paste("sealspres~",paste(allcovars,collapse="+"),"+I(",paste(covars,collapse="^2)+I("),"^2)",sep=""))

fullmodel<-glm(formula=fmlfull,data=data,family="binomial")
summary(fullmodel)

# Find top model via stepwise, using BIC as criterion to avoid AIC's over-parametrization issues
k<-log(nrow(data))
stepm<-step(fullmodel,k=k)	#using BIC as the selection criterion 
summary(stepm)

#Stepm is not an ideal fit... And...
#Rearrange the factor for AllRegion to use WRSS as the ref and do the shotgun below
data$AllRegion<-ifelse(data$Region=="ERS","ERS",ifelse(data$drygalski==0,"WRSS","WRSN"))
data$regionOrder<-ifelse(data$Region=="ERS",2,ifelse(data$drygalski==0,1,3))
data$AllRegion<-as.factor(data$AllRegion)
data$AllRegion<-reorder(data$AllRegion,data$regionOrder)

#Our top model is really:
topfml<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgdi_edge + lgdi_glac + lgdi_deep + 
				lgicewidth + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgdi_empe^2) + I(lgdi_deep^2) + AllRegion")
topmdl<-glm(formula=topfml,data=data,family="binomial");summary(topmdl)

#Adding linear covars for those quadratics that did not include them to see if it's worth adding them 
#(i.e., to see if it makes sense not to force the parabola to have min/max at 0)
topfml_l<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgdi_edge + lgdi_glac + lgdi_deep + lgdi_cont + lgdi_adpe +
				lgnear_empe + lgicewidth + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2) + AllRegion") 
topmdl_l<-glm(formula=topfml_l,data=data,family="binomial");summary(topmdl_l)
#RESULT: not worth adding the linear effects to those quadratics that do not heve them

# Remove ERS and redo, to see if ERS is different from WRS
datawrs<-subset(data,Region=="WRS")
topmdl_wrs<-glm(formula=topfml,data=datawrs,family="binomial");summary(topmdl_wrs)

#we tried using shore on both the overal and the WRS models - not best.

#add some interactions:
topfml_rint<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgdi_adpe + lgdi_glac + lgicewidth +  
				AllRegion*lgdi_edge + AllRegion*lgdi_deep + AllRegion*lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_rint<-glm(formula=topfml_rint,data=data,family="binomial");summary(topmdl_ra)

#explore using lgratio_edge:
topfml_ra<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgdi_glac + lgicewidth +  
				AllRegion*lgratio_edge + AllRegion*lgdi_deep + AllRegion*lgdi_empe + I(lgdi_cont^2) + I(lgdi_adpe^2) + I(lgnear_empe^2) + I(lgdi_deep^2)")
topmdl_ra<-glm(formula=topfml_ra,data=data,family="binomial");summary(topmdl_ra)

#using logratio but no interactions - this is the top model, but with the non-important linear effects for lgdi_adpe and lgdi_empe, and using lgratio_edge:
topfml_nint<-as.formula("sealspres ~ tomnodaoi + crack + Year + mlsearch + lgratio_edge + lgdi_glac + lgicewidth + lgdi_empe + 
				lgdi_cont + I(lgdi_cont^2) + lgdi_adpe + I(lgdi_adpe^2) + lgnear_empe + I(lgnear_empe^2) + lgdi_deep + I(lgdi_deep^2) + AllRegion")
topmdl_nint<-glm(formula=topfml_nint,data=data,family="binomial");summary(topmdl_nint)

#remove the non-important linear effects (attention: here it is only lgdi_adpe and lgdi_cont):
topfml_nint<-as.formula("sealspres ~ tomnodaoi + crack + mlsearch + Year + lgdi_glac + I(lgdi_cont^2) + lgicewidth + lgdi_deep + I(lgdi_deep^2) + I(lgdi_adpe^2) +  
				 lgratio_edge + lgdi_empe + lgnear_empe + I(lgnear_empe^2) + AllRegion") 
topmdl_nint<-glm(formula=topfml_nint,data=data,family="binomial");summary(topmdl_nint)

#Is there an interaction between distance and size of penguin colonies? A: NO
topfml_nint_peng<-as.formula("sealspres ~ tomnodaoi + crack + mlsearch + Year + lgdi_glac + I(lgdi_cont^2) + lgicewidth + lgdi_deep + I(lgdi_deep^2) + I(lgdi_adpe^2) +  
				lgratio_edge + lgdi_empe + lgnear_empe + I(lgnear_empe^2) + AllRegion + lgdi_adpe*lgnear_adpe + lgdi_empe*lgnear_empe") 
topmdl_nint_peng<-glm(formula=topfml_nint_peng,data=data,family="binomial");summary(topmdl_nint_peng)

#Nadav wants LRTs for the quadratics that also have linear effects, independently. So:
library(lmtest)
deepreffml<-as.formula("sealspres ~ tomnodaoi + crack + mlsearch + Year + lgdi_glac + I(lgdi_cont^2) + lgicewidth + I(lgdi_adpe^2) +  
				lgratio_edge + lgdi_empe + lgnear_empe + I(lgnear_empe^2) + AllRegion")
deeprefmdl<-glm(formula=deepreffml,data=data,family="binomial"); lrtest(deeprefmdl,topmdl_nint)

nempreffml<-as.formula("sealspres ~ tomnodaoi + crack + mlsearch + Year + lgdi_glac + I(lgdi_cont^2) + lgicewidth + lgdi_deep + I(lgdi_deep^2) + I(lgdi_adpe^2) +  
				lgratio_edge + lgdi_empe + AllRegion")
nemprefmdl<-glm(formula=nempreffml,data=data,family="binomial"); lrtest(nemprefmdl,topmdl_nint)

#BEST MODEL: region interactions with lgratio_edge, and compare to model without regions (topmdl_nint) to show how much of an improvement 
summary(topmdl_ra);summary(topmdl_nint)
#How important is to add the interaction?
(spin<-anova(topmdl_nint,topmdl_ra))
1-pchisq(spin[2,4],9)  

#GOF for the models
logitgof(data$sealspres,fitted(topmdl_ra))
logitgof(data$sealspres,fitted(topmdl_nint))
#In both cases there is an interval with no data, so, though not as good a test...
logitgof(data$sealspres,fitted(topmdl_ra),g=5)
logitgof(data$sealspres,fitted(topmdl_nint),g=6)

write.csv(summary(topmdl_nint)$coefficients,file=paste(basepth,"/paper/table3.csv",sep=""))
write.csv(summary(topmdl_ra)$coefficients,file=paste(basepth,"/paper/table4.csv",sep=""))

#calculate the deviance explained by each covariate in topmdl_ra and topmdl_nint
nint_vars<-c("tomnodaoi","crack","mlsearch","Year","lgdi_glac","I(lgdi_cont^2)","lgicewidth","lgdi_deep+I(lgdi_deep^2)",
		"I(lgdi_adpe^2)","lgratio_edge","lgdi_empe","lgnear_empe+I(lgnear_empe^2)","AllRegion")

ra_vars<-c("tomnodaoi","crack","Year","mlsearch","lgdi_glac","lgicewidth","AllRegion:lgratio_edge","AllRegion:lgdi_deep","AllRegion:lgdi_empe",
		"AllRegion","I(lgdi_cont^2)","I(lgdi_adpe^2)","I(lgnear_empe^2)","I(lgdi_deep^2)")

devexp_nint<-estExplaniedDeviance(topmdl=topmdl_nint,vars=nint_vars,data=data)
devexp_nint<-devexp_nint[order(devexp_nint$deviance_exp,decreasing=TRUE),]

devexp_ra<-estExplaniedDeviance(topmdl=topmdl_nint,vars=ra_vars,data=data)
devexp_ra<-devexp_ra[order(devexp_ra$deviance_exp,decreasing=TRUE),]

##PLOTS: 
#recentering results by using a model with only region effects
fml_region<-as.formula("sealspres ~ AllRegion"); mdl_region<-glm(formula=fml_region,data=data,family="binomial");summary(mdl_region)
qvr<-summary(mdl_region)$coefficients;
wrssInt<-qvr[1,1];ersInt<-qvr[2,1]+wrssInt;wrsnInt<-qvr[3,1]+wrssInt;intsC<-c(wrssInt,wrsnInt,ersInt);names(intsC)<-c("WRSS","WRSN","ERS")

# Partial dependence of individual effects
#distance to deep
pdf<-makedf(reg="WRSN",var="lgdi_deep",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_deep",intsC=intsC,reg="WRSN",conf=86);plotdfn$Region<-"WRSN"		#conf=68 or any other value will yield 1SD, 0 means no SD
pdf<-makedf(reg="WRSS",var="lgdi_deep",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_deep",intsC=intsC,reg="WRSS",conf=86);plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgdi_deep",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_deep",intsC=intsC,reg="ERS",conf=86);plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgdi_deep",intsC=intsC,reg="ALL",conf=86);nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
plotdf$lgdi_deep<-plotdf$lgdi_deep-1
plotdf<-subset(plotdf,lgdi_deep>0)
p<-ggplot(data=plotdf,aes(x=lgdi_deep,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Distance to deep waters)")
dev.new();print(p)
names(plotdf)<-gsub("lgdi_deep","CovariateValue",names(plotdf))
plotdf$CovarName<-"LGDI_DEEP"
regionsdf<-plotdf
names(nreg)<-gsub("lgdi_deep","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to deep waters)"
regionsdf<-plotdf
nregdf<-nreg

#lgratio_edge
pdf<-makedf(reg="WRSN",var="lgratio_edge",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgratio_edge",intsC=intsC,reg="WRSN",conf=86);plotdfn$Region<-"WRSN"
pdf<-makedf(reg="WRSS",var="lgratio_edge",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgratio_edge",intsC=intsC,reg="WRSS",conf=86);plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgratio_edge",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgratio_edge",intsC=intsC,reg="ERS",conf=86);plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgratio_edge",intsC=intsC,reg="ALL",conf=86);nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
plotdf<-subset(plotdf,lgratio_edge>0)
p<-ggplot(data=plotdf,aes(x=lgratio_edge,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Dist. to ice edge/ice width)")
dev.new();print(p)
names(plotdf)<-gsub("lgratio_edge","CovariateValue",names(plotdf))
plotdf$CovarName<-"LGICE_RATIO"
regionsdf<-rbind(regionsdf,plotdf)
names(nreg)<-gsub("lgratio_edge","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to ice edge/ice width)"
nregdf<-rbind(nregdf,nreg)

#distance to emperor penguins
pdf<-makedf(reg="WRSN",var="lgdi_empe",data=data)
plotdfn<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_empe",intsC=intsC,reg="WRSN",conf=86);plotdfn$Region<-"WRSN"
pdf<-makedf(reg="WRSS",var="lgdi_empe",data=data)
plotdfs<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_empe",intsC=intsC,reg="WRSS",conf=86);plotdfs$Region<-"WRSS"
pdf<-makedf(reg="ERS",var="lgdi_empe",data=data)
plotdfe<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_empe",intsC=intsC,reg="ERS",conf=86);plotdfe$Region<-"ERS"
nreg<-makePredictions(df=pdf,mdlfit=topmdl_nint,varnam="lgdi_empe",intsC=intsC,reg="ALL",conf=86);nreg$Region<-"ALL"

plotdf<-rbind(plotdfn,plotdfs);plotdf<-rbind(plotdf,plotdfe)
plotdf<-subset(plotdf,lgdi_empe>0)
p<-ggplot(data=plotdf,aes(x=lgdi_empe,y=predicted)) + 
		geom_line(color="magenta",size=1) + geom_errorbar(aes(ymin=ymin,ymax=ymax)) +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_wrap(~Region,ncol=3) +
		labs(y="Prob. seal presence", x="Log(Distance to EMPE)")
dev.new();print(p)
names(plotdf)<-gsub("lgdi_empe","CovariateValue",names(plotdf))
plotdf$CovarName<-"LGDI_EMPE"
regionsdf<-rbind(regionsdf,plotdf)
names(nreg)<-gsub("lgdi_empe","CovariateValue",names(nreg))
nreg$CovarName<-"Log(Dist. to EMPE)"
nregdf<-rbind(nregdf,nreg)


##combining
regionsdf$RegionName<-ifelse(regionsdf$Region=="ERS","Eastern Ross Sea",
		ifelse(regionsdf$Region=="WRSS","Western Ross Sea South","Western Ross Sea North"))

regionsdf$Covar<-ifelse(regionsdf$CovarName=="LGDI_DEEP","Log(Distance to deep waters)",
		ifelse(regionsdf$CovarName=="LGDI_EMPE","Log(Distance to EMPE colony)","Log(Dist. to edge/ice width)"))

regionsdf$CovOrder<-ifelse(regionsdf$CovarName=="LGDI_DEEP",1,ifelse(regionsdf$CovarName=="LGDI_EMPE",2,3))

regionsdf<-within(regionsdf,{Covar<-reorder(Covar,CovOrder)})

p<-ggplot(data=regionsdf,aes(x=CovariateValue,y=predicted)) + 
		geom_line(color="black",size=1.3) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_grid(RegionName~Covar,scales="free") +
		labs(y="Probability of seal presence", x="Predictor Value") +
		scale_x_continuous(breaks=c(0:8)) +
		theme(axis.text.x=element_text(size=24, color="black")) + theme(axis.title.x=element_text(size=28,face="bold")) +
		theme(axis.text.y=element_text(size=24, color="black")) + theme(axis.title.y=element_text(size=28,face="bold")) +
		theme(strip.text.x=element_text(size=24)) + theme(strip.text.y=element_text(size=24))

jpeg(filename=paste(basepth,"/paper/plots/Plasticity.jpg",sep=""),width=1200,height=1100,quality=100)
print(p)
dev.off()

### DO NOT USE This:
nregdf<-subset(nregdf,CovariateValue>0);nregdf$RegionName<-"All Regions"
alldata<-rbind(regionsdf,nregdf)
alldata$CovarLetter<-ifelse(alldata$CovarName=="Log(Dist. to EMPE)","A",
				ifelse(alldata$CovarName=="Log(Dist. to deep waters)","B","C"))
p<-ggplot(data=alldata,aes(x=CovariateValue,y=predicted)) + 
		geom_line(color="black",size=2) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		facet_grid(RegionName~CovarLetter,scales="free_x") +
		labs(y="Probability of seal presence", x="Predictor Value") +
		theme(axis.text.x=element_text(size=24, color="black")) + theme(axis.title.x=element_text(size=28,face="bold")) +
		theme(axis.text.y=element_text(size=24, color="black")) + theme(axis.title.y=element_text(size=28,face="bold")) +
		theme(strip.text.x=element_text(size=24)) + theme(strip.text.y=element_text(size=24))
###################


#then the non-interacting effects: lgdi_adpe, lgdi_glac, lgicewidth, lgdi_cont, lgnear_empe
#lgdi_cont
#recentering... but this time to the entire region
fml_region<-as.formula("sealspres ~ 1"); mdl_region<-glm(formula=fml_region,data=data,family="binomial");summary(mdl_region)
qvr<-summary(mdl_region)$coefficients;
wrssInt<-qvr[1,1];intsC<-c(wrssInt);names(intsC)<-c("WRSS")

pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgdi_cont",data=data);pdf$AllRegion<-"WRSS"	#reg=c("WRSN","WRSS","ERS")
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_cont",intsC=intsC,reg="WRSS",conf=68);plotdf$Region<-"All Regions"	#conf=68 or any other value will yield 1SD, 0 means no SD
nintpdf<-plotdf;nintpdf$Variable<-"A";names(nintpdf)<-gsub("lgdi_cont","value",names(nintpdf))


p<-ggplot(data=subset(plotdf,lgdi_cont>=0),aes(x=lgdi_cont,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability of seal presence", x="Log(Distance to shore)") +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold"))

jpeg(filename=paste(basepth,"/paper/plots/DistToContinent.jpg",sep=""),width=400,height=400,quality=100)
print(p)
dev.off()

#lgdi_adpe
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgdi_adpe",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_adpe",intsC=intsC,reg="WRSS",conf=68);plotdf$Region<-"All Regions"
tdf<-plotdf;tdf$Variable<-"B";names(tdf)<-gsub("lgdi_adpe","value",names(tdf))
tdf<-subset(tdf,value>0)
nintpdf<-rbind(nintpdf,tdf)

p<-ggplot(data=subset(plotdf,lgdi_adpe>=0),aes(x=lgdi_adpe,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability of seal presence", x="Log(Distance to ADPE colony)") +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold"))

jpeg(filename=paste(basepth,"/paper/plots/DistToNearestADPE.jpg",sep=""),width=420,height=400,quality=100)
print(p)
dev.off()

#lgdi_glac
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgdi_glac",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgdi_glac",intsC=intsC,reg="WRSS",conf=68);plotdf$Region<-"All Regions"
tdf<-plotdf;tdf$Variable<-"C";names(tdf)<-gsub("lgdi_glac","value",names(tdf))
nintpdf<-rbind(nintpdf,tdf)

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

#lgicewidth
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgicewidth",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgicewidth",intsC=intsC,reg="WRSS",conf=68);plotdf$Region<-"All Regions"
tdf<-plotdf;tdf$Variable<-"D";names(tdf)<-gsub("lgicewidth","value",names(tdf))
nintpdf<-rbind(nintpdf,tdf)

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

#lgnear_empe
pdf<-makedf(reg=c("WRSN","WRSS","ERS"),var="lgnear_empe",data=data);pdf$AllRegion<-"WRSS"
plotdf<-makePredictions(df=pdf,mdlfit=topmdl_ra,varnam="lgnear_empe",intsC=intsC,reg="WRSS",conf=68);plotdf$Region<-"All Regions"
tdf<-plotdf;tdf$Variable<-"E";names(tdf)<-gsub("lgnear_empe","value",names(tdf))
nintpdf<-rbind(nintpdf,tdf)

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

#combined
p<-ggplot(data=nintpdf,aes(x=value,y=predicted)) + 
		geom_line(color="black",size=1) + geom_ribbon(aes(ymin=ymin,ymax=ymax),color="light gray",alpha=0.5) +
		theme_bw() +
		scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
		labs(y="Probability of seal presence", x="Predictor value") +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold")) +
		facet_wrap(~Variable,scales="free_x",ncol=3) +
		theme(strip.text.x=element_text(size=24)) + theme(strip.background=element_blank(),strip.text.x=element_blank())

jpeg(filename=paste(basepth,"/paper/plots/Non_interacting_effects.jpg",sep=""),width=1100,height=750,quality=100)
print(p)
dev.off()



###Let's examine the relationship with distance to deep waters 
data0<-subset(data,sealspres==0);data1<-subset(data,sealspres==1)
cut0<-table(cut(data0$lgdi_deep,breaks=seq(0,5.6,0.7)));cut1<-table(cut(data1$lgdi_deep,breaks=seq(0,5.6,0.7)))
deepbar<-data.frame(LgIntDeep=rep(names(cut0),2),count=c(as.numeric(cut0),count1=as.numeric(cut1)),Presence=c(rep("No",8),rep("Yes",8)))
p<-ggplot(deepbar,aes(x=LgIntDeep,y=count)) + geom_bar(stat="identity",aes(fill=Presence),position="stack") +
		scale_fill_manual(values=c("gray","black")) + theme_bw() +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold")) +
		theme(legend.text=element_text(size=22, color="black"), legend.title=element_text(size=22, color="black")) +
		labs(x="Intervals of Log(Distance to deep waters)",y="Cell count",fill="Seal presence")

jpeg(filename=paste(basepth,"/paper/plots/Pattern_deep.jpg",sep=""),width=950,height=550,quality=100)
print(p)
dev.off()


#Similarly, examine the lgdi_adpe quadrature for lgdi_adpe<0 for WRSS
data0s<-subset(data0,AllRegion=="WRSS");data1s<-subset(data1,AllRegion=="WRSS");
cut0<-table(cut(data0s$lgdi_adpe,breaks=c(seq(-3.6,4.5,by=0.9),Inf)));cut1<-table(cut(data1s$lgdi_adpe,breaks=c(seq(-3.6,4.5,by=0.9),Inf)))
adpebar<-data.frame(LgIntADPE=rep(names(cut0),2),count=c(as.numeric(cut0),count1=as.numeric(cut1)),Presence=c(rep("No",10),rep("Yes",10)))
p<-ggplot(adpebar,aes(x=LgIntADPE,y=count)) + geom_bar(stat="identity",aes(fill=Presence),position="stack") +
		scale_fill_manual(values=c("gray","black")) + theme_bw() +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold")) +
		theme(legend.text=element_text(size=22, color="black"), legend.title=element_text(size=22, color="black")) +
		labs(x="Intervals of Log(Distance to deep waters)",y="Cell count",fill="Seal presence")

jpeg(filename=paste(basepth,"/paper/plots/Pattern_adpe.jpg",sep=""),width=1200,height=550,quality=100)
print(p)
dev.off()

# Known biology: do a bar graph of percent cells with seals for crack and no-crack from raw data, or predict from top model with all else constant (should be very similar)
crdata<-aggregate(as.formula("sealspres~crack+AllRegion"),data=data,FUN=sum,na.rm=T)
nrdata<-aggregate(as.formula("sealspres~crack+AllRegion"),data=data,FUN=NROW);names(nrdata)<-c("crack","AllRegion","nrecs")
crackpct<-merge(crdata,nrdata,by=c("crack","AllRegion"),all.x=T)
crackpct$percent<-round(crackpct$sealspres*100/crackpct$nrecs)
crackpct$CrackPres<-ifelse(crackpct$crack==0,"No crack","Crack")
crackpct$Region<-ifelse(crackpct$AllRegion=="ERS","Eastern Ross Sea",ifelse(crackpct$AllRegion=="WRSN","W. Ross Sea North","W. Ross Sea South"))
p<-ggplot(data=crackpct,aes(x=Region,y=percent)) + 
		geom_bar(stat="identity",position="dodge",aes(fill=CrackPres)) +
		scale_fill_grey() + theme_bw() +
		theme(axis.text.x=element_text(size=22, color="black")) + theme(axis.title.x=element_text(size=24,face="bold")) +
		theme(axis.text.y=element_text(size=22, color="black")) + theme(axis.title.y=element_text(size=26,face="bold")) +
		theme(legend.text=element_text(size=22, color="black"), legend.title=element_text(size=22, color="black")) +
		labs(x="",y="Probability of Seal Presence",fill="Presence of cracks")

jpeg(filename=paste(basepth,"/paper/plots/Crack_by_region.jpg",sep=""),width=1000,height=550,quality=100)
print(p)
dev.off()
