library(COSTcore)
library(tidyr)
library(dplyr)
library(sf)
library(ggplot2)
library(mapdata)
library(maptools)
library(raster)
library(rasterVis)

#listspp
codespp<-read.table("../../data/refTables/ISIH-19926-espece_fao-2018.txt",header=T,sep=";",fileEncoding="Latin1") 
codespp1<-codespp%>%filter(grepl("Raie",ESPF_FRA_LIB)|substr(ESPF_COD,1,2)=="RJ"|ESPF_COD%in%c("RAJ","SRX","SKA","JDP"))%>%arrange(ESPF_COD)
listparent<-unique(codespp1$ESPF_PARENT_COD)
codespp<-codespp%>%filter(grepl("ray",ESPF_ENG_LIB)|grepl("skate",ESPF_ENG_LIB)|grepl("Raie",ESPF_FRA_LIB)|substr(ESPF_COD,1,2)=="RJ"|ESPF_COD%in%c("RAJ","SRX","SKA","JDP")|ESPF_COD%in%listparent)%>%arrange(ESPF_COD)
scival<-unique(codespp$NOM_VALIDE)
scival<-scival[scival!=""]
#getaphiaid
library("taxizesoap")
codespp$worms<-""
for(i in 1:nrow(codespp)){
		   print(paste("###########################################"))
		   print(paste(i,codespp$ESPF_SCI_LIB[i]))
		   codespp$worms[i]<-as.numeric(get_wormsid(searchterm=codespp$ESPF_SCI_LIB[i])[1])
		   if(is.na(codespp$worms[i])&codespp$NOM_VALIDE[i]!=""){
				print("NOM VALIDE")
		   		codespp$worms[i]<-as.numeric(get_wormsid(searchterm=codespp$NOM_VALIDE[i])[1])
			}
		   print(worms_name(ids=codespp$worms[i]))
}
#add TSN
tsn<-read.csv("./TSN/longnames.csv")
tsn<-tsn%>%transmute(ESPF_SCI_LIB=completename,tsn=tsn)
codespp<-left_join(codespp,tsn)
save(codespp,file="codespp.rdata")
#format codespp
listspp<-codespp%>%transmute(Name=ESPF_ENG_LIB,FAO=ESPF_COD,WORMS=worms)%>%distinct()

load("codespp.rdata")

if(F){

	#ENG CS
	pipo<-data.frame(line=readLines("./ENG_COST/cs.txt"))
	tr<-pipo%>%filter(grepl("TR,",line))%>%
		tidyr::separate(line,into=as.character(1:17),sep=",")%>%
		select(-1)
	hh<-pipo%>%filter(grepl("HH,",line))%>%
		tidyr::separate(line,into=as.character(1:29),sep=",")%>%
		select(-1)
	sl<-pipo%>%filter(grepl("SL,",line))%>%
		tidyr::separate(line,into=as.character(1:18),sep=",")%>%
		select(-1)
	hl<-pipo%>%filter(grepl("HL,",line))%>%
		tidyr::separate(line,into=as.character(1:18),sep=",")%>%
		select(-1)

	#cost and quick correction
	library(COSTcore)
	#tr
	names(tr)<-names(csData()@tr)
	#hh
	hh<-cbind(hh[,1:7],data.frame(uu=""),hh[,8:ncol(hh)])
	names(hh)<-names(csData()@hh)
	hh$latIni<-as.numeric(hh$latIni)
	hh$latFin<-as.numeric(hh$latFin)
	hh$lonIni<-as.numeric(hh$lonIni)
	hh$lonFin<-as.numeric(hh$lonFin)
	#hh%>%group_by(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum)%>%summarise(n=n())%>%filter(n>1)#select(1:7)%>%distinct()
	hhdouble<-hh%>%select(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum)%>%data.frame()%>%duplicated()
	hh<-hh%>%filter(!hhdouble)

	#sl
	names(sl)<-names(csData()@sl)
	sl0<-sl
	sl0$sex<-sl$catchCat
	sl0$catchCat<-sl$landCat
	sl0$landCat<-sl$commCatScl
	sl0$commCatScl<-"EU"
	sl0$wt<-as.numeric(sl0$wt)
	sl0$subSampWt<-as.numeric(sl0$subSampWt)
	sl<-sl0
	#sldouble<-sl%>%select(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum,spp,sex,catchCat)%>%data.frame()%>%duplicated()
	#sl[sldouble,]
	#sl%>%filter(trpCode==62744&staNum==999&grepl("Raja clav",spp))
	sl<-sl%>%tidyr::separate(spp,into=c("genre","esp"),sep=" ",remove=FALSE)%>%tidyr::unite(spp,c("genre","esp"),sep=" ")
	sl<-sl%>%group_by(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum,spp,catchCat,landCat,commCatScl,commCat,subSampCat,sex)%>%
			summarise(wt=sum(wt,na.rm=T),subSampwt=sum(subSampWt,na.rm=T),lenCode="mm")%>%ungroup()
	#remove no ray
	listspuk<-sort(unique(sl$spp))[-c(3,4,5,8,14,15,16)]
	sl<-sl%>%filter(spp%in%listspuk)

	#hl
	hl<-hl[,-9]
	names(hl)<-names(csData()@hl)
	hl$lenNum<-as.numeric(hl$lenNum)
	hl<-hl%>%tidyr::separate(spp,into=c("genre","esp"),sep=" ",remove=FALSE)%>%tidyr::unite(spp,c("genre","esp"),sep=" ")
	hl$commCatScl<-"EU"
	hl$commCat<-"English"
	#hldouble<-hl%>%select(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum,catchCat,spp,landCat,commCatScl,sex,lenCls)%>%data.frame()%>%duplicated()
	#hl[hldouble,]%>%head()
	#hl%>%filter(trpCode==62744&staNum==999&grepl("Raja clav",spp))%>%arrange(sex,lenCls)
	hl<-hl%>%group_by(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum,spp,catchCat,landCat,commCatScl,commCat,subSampCat,sex,lenCls)%>%
			summarise(lenNum=sum(lenNum,na.rm=T))%>%ungroup()

	hl<-hl%>%filter(spp%in%listspuk)

	cs<-csData(tr=tr,hh=hh,sl=data.frame(sl),hl=data.frame(hl))
	saveRDS(cs,file="csuk.rds")



	#CL
	dat1<-data.frame()
	for(i in 2017:2000){
		year<-as.character(i)
		print(paste(i,year))
		#print(load(paste0("../../../../dataCOST/sacrois_v335_",year,"_navday_newmet6.rdata")))
		dat<-readRDS(paste0("../../dataCOST/cl",year,".rds"))
		dat<-dat@cl
		#selection des lignes avec des crevettes dans le 27.7.d
		dat<-tbl_df(dat)%>%filter(area%in%c("27.7.d","27.7.e","27.4.c","27.4.b")&taxon%in%codespp$ESPF_COD) 
		dat1<-rbind(dat1,data.frame(dat))
	}
	saveRDS(dat1,file="cl.rds")
	#CE
	dat1<-data.frame()
	for(i in 2017:2000){
		year<-as.character(i)
		print(paste(i,year))
		#print(load(paste0("../../../../dataCOST/sacrois_v335_",year,"_navday_newmet6.rdata")))
		dat<-readRDS(paste0("../../dataCOST/ce",year,".rds"))
		dat<-dat@ce
		#selection des lignes avec des crevettes dans le 27.7.d
		dat<-tbl_df(dat)%>%filter(area%in%c("27.7.d","27.7.e","27.4.c","27.4.b"))#&taxon%in%codespp$ESPF_COD) 
		dat1<-rbind(dat1,data.frame(dat))
	}
	saveRDS(dat1,file="ce.rds")
	#CS
	dat1<-csData()
	for(i in 2017:2000){
		year<-as.character(i)
		print(paste(i,year))
		#print(load(paste0("../../../../dataCOST/sacrois_v335_",year,"_navday_newmet6.rdata")))
		dat<-readRDS(paste0("../../dataCOST/cs",year,".rds"))
		#dat<-dat@ce
		#selection des lignes avec des crevettes dans le 27.7.d
		dat<-subset(dat,area%in%c("27.7.d","27.7.e","27.4.c","27.4.b"),table="hh")
		dat<-subset(dat,spp%in%c(unique(codespp$ESPF_SCI_LIB),scival),table="sl")
		#dat<-tbl_df(dat)%>%filter(area%in%c("27.7.d","27.7.e","27.4.c","27.4.b"))#&taxon%in%codespp$ESPF_COD) 
		dat1<-rbind2(dat1,dat)
	}
	saveRDS(dat1,file="cs.rds")

	#campagne
	path<-"./CGFS/"
	fich<-paste0(path,"data.csv")
	pipo<-readLines(fich,n=-1)
	id<-which(substr(pipo,1,3)=="Rec")
	hh<-pipo[substr(pipo,1,2)=="HH"]
	hhname<-strsplit(pipo[id[1]],split=",")%>%unlist()
	hh1<-tidyr::separate(data.frame(nom=hh),nom,into=hhname,sep=",")
	hl<-pipo[substr(pipo,1,2)=="HL"]
	hlname<-strsplit(pipo[id[2]],split=",")%>%unlist()
	hl1<-tidyr::separate(data.frame(nom=hl),nom,into=hlname,sep=",")
	caname<-strsplit(pipo[id[3]],split=",")%>%unlist()
	ca<-pipo[substr(pipo,1,2)=="CA"]
	ca1<-tidyr::separate(data.frame(nom=ca),nom,into=caname,sep=",")
	hl1<-hl1%>%subset(ValidAphiaID%in%codespp$worms)
	cgfs<-left_join(hl1%>%dplyr::select(-RecordType),hh1%>%dplyr::select(-RecordType))
	cgfs0<-hh1

	#BTS
	path<-"./BTS/"
	fich<-paste0(path,"data.csv")
	pipo<-readLines(fich,n=-1)
	id<-which(substr(pipo,1,3)=="Rec")
	hh<-pipo[substr(pipo,1,2)=="HH"]
	hhname<-strsplit(pipo[id[1]],split=",")%>%unlist()
	hh1<-tidyr::separate(data.frame(nom=hh),nom,into=hhname,sep=",")
	hl<-pipo[substr(pipo,1,2)=="HL"]
	hlname<-strsplit(pipo[id[2]],split=",")%>%unlist()
	hl1<-tidyr::separate(data.frame(nom=hl),nom,into=hlname,sep=",")
	caname<-strsplit(pipo[id[3]],split=",")%>%unlist()
	ca<-pipo[substr(pipo,1,2)=="CA"]
	ca1<-tidyr::separate(data.frame(nom=ca),nom,into=caname,sep=",")
	hl1<-hl1%>%filter(SpecCodeType=="W")%>%subset(ValidAphiaID%in%codespp$worms)
	bts<-left_join(hl1%>%dplyr::select(-RecordType),hh1%>%dplyr::select(-RecordType))
	bts0<-hh1

	#NSBTS
	path<-"./NSBTS/"
	fich<-paste0(path,"data.csv")
	pipo<-readLines(fich,n=-1)
	id<-which(substr(pipo,1,3)=="Rec")
	hh<-pipo[substr(pipo,1,2)=="HH"]
	hhname<-strsplit(pipo[id[1]],split=",")%>%unlist()
	hh1<-tidyr::separate(data.frame(nom=hh),nom,into=hhname,sep=",")
	hl<-pipo[substr(pipo,1,2)=="HL"]
	hlname<-strsplit(pipo[id[2]],split=",")%>%unlist()
	hl1<-tidyr::separate(data.frame(nom=hl),nom,into=hlname,sep=",")
	caname<-strsplit(pipo[id[3]],split=",")%>%unlist()
	ca<-pipo[substr(pipo,1,2)=="CA"]
	ca1<-tidyr::separate(data.frame(nom=ca),nom,into=caname,sep=",")
	hl1<-hl1%>%subset(ValidAphiaID%in%codespp$worms)
	nsbts<-left_join(hl1%>%dplyr::select(-RecordType),hh1%>%dplyr::select(-RecordType))
	nsbts0<-hh1

	camp<-rbind(cgfs,bts,nsbts)%>%filter(Year>=2000)
	camp0<-rbind(cgfs0,bts0,nsbts0)%>%filter(Year>=2000)

	#filtering space
	#add rect loc
	load("../../data/refTables/geosihsextant/rect.rdata")
	rect<-rect%>%filter(F_DIVISION%in%c("27.7.d","27.7.e","27.4.c","27.4.b"))
	camp<-camp%>%filter(StatRec%in%rect$ICESNAME)
	camp0<-camp0%>%filter(StatRec%in%rect$ICESNAME)
	saveRDS(camp,file="camp.rds")	
	saveRDS(camp0,file="camp0.rds")	

stop("ici")

	


	nomfich<-"CL_FRA_2000-2017.csv"
	fwrite(dat1,file="CL_FRA_2000-2017.csv")
	pipo<-readr::format_csv(dat1)
	readr::write_file(pipo,nomfich)
	write.csv(dat1,file=nomfich)

}

metcl<-dat1%>%transmute(met=foCatEu6)%>%distinct()%>%
		tidyr::separate("met",c("gear","spp","m1","m2"),
				remove=F,sep="_")%>%
		#remove metier wich selective device
		mutate(
		m2=ifelse(grepl("<",m1),m1,m2),
		m1=ifelse(grepl("<",m1),0,m1),
		m2=ifelse(grepl(">",m1),99999,m2),
		m1=gsub(">","",m1),m1=gsub("=","",m1),
		m2=gsub(">","",m2),m2=gsub("=","",m2),
		m1=gsub("<","",m1),m1=gsub("=","",m1),
		m2=gsub("<","",m2),m2=gsub("=","",m2),
		m1=ifelse(is.na(m1),0,m1),
		m2=ifelse(is.na(m2),0,m2)
		)

#read met list ices
metices<-read.csv("metices.csv")%>%
		transmute(met=Metier.Fleet.Name)%>%distinct()%>%
		filter(!grepl("FDF",met))%>%
		tidyr::separate("met",c("gear","spp","mesh","m2"),
				remove=F,sep="_") %>%
		#remove metier wich selective device
		filter(m2==0)%>%dplyr::select(-m2)%>%
		tidyr::separate("mesh",c("m1","m2"),
				remove=F,sep="-") %>%
		mutate(
		m2=ifelse(grepl("<",mesh),m1,m2),
		m1=ifelse(grepl("<",mesh),0,m1),
		m2=ifelse(grepl(">",mesh),99999,m2),
		m1=gsub(">","",m1),m1=gsub("=","",m1),
		m2=gsub(">","",m2),m2=gsub("=","",m2),
		m1=gsub("<","",m1),m1=gsub("=","",m1),
		m2=gsub("<","",m2),m2=gsub("=","",m2),
		m1=ifelse(is.na(m1),0,m1),
		m2=ifelse(is.na(m2),0,m2)
		)


	metarea$ices<-""
	for(i in 1:nrow(metarea)){
		metarea$ices[i]<-findmet(metarea[i,],reftab)
	}
	rez<-metarea%>%select(area,met,ices)




load("sacrois.rdata")
#clean some info
dat1<-dat1%>%filter(QUANT_POIDS_VIF_SACROIS<9999)
plot(dat1$QUANT_POIDS_VIF_SACROIS,n=100)
dat1%>%filter(QUANT_POIDS_VIF_SACROIS>2000)%>%View()

#format data
dat2<-dat1%>%transmute(idmaree=MAREE_ID,
		    year=substr(MAREE_DATE_RET,7,10),month=substr(MAREE_DATE_RET,4,5),
		    div=SECT_COD_SACROIS_NIV3,rect=SECT_COD_SACROIS_NIV5,
		    nav=NAVS_COD,port=LIEU_COD_RET_SACROIS,
		    #met=METIER_DCF_6_COD,
		    w=QUANT_POIDS_VIF_SACROIS,spp=ESP_COD_FAO,
		    tps=TP_NAVIRE_SACROIS,prix=MONTANT_EUROS_SACROIS)%>%
		filter(spp=="CSH")%>%group_by(idmaree)%>%mutate(n=1/length(idmaree))%>%ungroup()
#add nav info
print(load("../../../../data/refTables/navires/navires.rdata"))
navires<-navires%>%transmute(nav=as.character(NAVS_COD),nom=CARN_NOM,year=as.character(fpcyear),
			  l=NAVP_LONGUEUR_PP,puissance=NAVP_PUISSANCE_AD)
dat2<-left_join(dat2,navires)
#add port info
print(load("../../../../data/refTables/geosihsextant/port.rdata"))
port<-port%>%transmute(port=LOC_LABEL,nomport=PORT_LIB,X,Y)%>%
	group_by(port)%>% summarise(nomport=first(nomport))
dat2<-left_join(dat2,port)

#some checks
dat2<-dat2%>%filter(w!=99999)

pipo<-dat2%>%ungroup()%>%mutate(dat=as.numeric(year)+as.numeric(month)/12)%>%group_by(nom,dat)%>%
	summarise(w=sum(w,na.rm=T),tps=sum(tps,na.rm=T),n=sum(n))%>%ungroup()
liste<-pipo%>%group_by(nom)%>%filter(max(w/n)>1000)%>%dplyr::select(nom)%>%distinct()
ggplot(pipo%>%filter(nom%in%liste$nom),aes(x=dat,y=w/n))+geom_point()+facet_wrap(~nom,scale="free")+geom_line()



n_distinct(dat1$MAREE_ID[dat1$ESP_COD_FAO=="CSH"])
sum(dat2$n)

#some checks
dat2%>%group_by(idmaree)%>%summarise(n=n_distinct(rect))%>%filter(n>2)
dat2%>%group_by(idmaree)%>%summarise(n=n_distinct(n))%>%filter(n>1)
dat2%>%filter(idmaree=="10674365")

dat2%>%filter(year==2013&month=="09"&rect=="27E9")

#monthly data
dat2<-dat2%>%group_by(year,month,div,rect,spp,port,nomport,nav,nom)%>%summarise(w=sum(w,na.rm=T),
								    prix=sum(prix,na.rm=T),
								    tps=sum(tps,na.rm=T),
								    nbmar=sum(n))
sum(dat2$nbmar)


#add baie de seine baie de somme id
dat2<-dat2%>%mutate(baie=ifelse(rect%in%c("28E8","28E9","28F0","27E8","27E9","27F0"),"seine",rect))
dat2<-dat2%>%mutate(baie=ifelse(rect%in%c("28F1","29E9","29F0","29F1","30E8","30F0","30F1"),"somme",baie))
dat2<-dat2%>%mutate(baie=ifelse(rect==""&port%in%c("ADP","FBL","DBL","XDK","XBL"),"somme",baie))
dat2<-dat2%>%mutate(baie=ifelse(rect==""&port%in%c("DCN","BCN","XLH","ACN","CCN"),"seine",baie))
dat2<-dat2%>%filter(baie!="")
table(dat2$nomport,dat2$baie)

pipo<-dat2%>%mutate(dat=as.numeric(year)+as.numeric(month)/12)%>%group_by(baie,dat)%>%
	summarise(w=sum(w,na.rm=T),tps=sum(tps,na.rm=T))
pipo<-dat2%>%mutate(dat=as.numeric(year))%>%group_by(baie,dat)%>%
	summarise(w=sum(w,na.rm=T),tps=sum(tps,na.rm=T),nbmar=sum(nbmar,na.rm=T),nbnav=n_distinct(nav))
ggplot(pipo,aes(x=dat,y=w/1000))+geom_point()+facet_grid(~baie,scale="free")+geom_line()
ggplot(pipo,aes(x=dat,y=nbmar))+geom_point()+facet_grid(~baie,scale="free")+geom_line()
ggplot(pipo,aes(x=dat,y=w/tps))+geom_point()+facet_grid(~baie,scale="free")+geom_line()
ggplot(pipo,aes(x=dat,y=nbnav))+geom_point()+facet_grid(~baie,scale="free")+geom_line()
#plot(dat2$w/dat2$tps)#,log="xy")

datpech<- dat2%>%mutate(dat=as.numeric(year))%>%group_by(baie,dat)%>%
	summarise(w=sum(w,na.rm=T),tps=sum(tps,na.rm=T),nbmar=sum(nbmar,na.rm=T),nbnav=n_distinct(nav))%>%
	filter(baie=="seine")%>%ungroup()

#build up env series for seine and somme
#poly for rect



#load strata
maprect<- rgdal::readOGR("/home/moi/ifremer/data/refTables/geosihsextant/Sextant-19926-15181899015/3", "IFR_F_RECTANGLE_CIEM")
maprect<-st_as_sf(maprect)%>%filter(F_DIVISION=="27.7.d"&SECT_COD%in%c("27F0"))

#build the coast
coast1<-map('worldHires',xlim=c(-1.4,.5),ylim=c(49.2,49.8),bg="light blue",col="grey20",lwd=0.1,fill=T)
IDs <- sapply(strsplit(coast1$names, ":"), function(x) x[1])
coast2<-fortify(coast1)
coast3<-map2SpatialPolygons(coast1,IDs=IDs)
#polymap3<-st_as_sfc(polymap2)
plot(maprect["ICESNAME"],axes=T)
plot(coast3,add=T)
#text( getSpPPolygonsLabptSlots(st_as_sfc(maprect)) , labels=as.character(mapstrate@data$Ifremer_id), cex=1)

#verif
#plot(dist2seine,add=T)
plot(mean(chl,na.rm=T))
points(0.2,49.45)
plot(maprect,add=T)


#chl
pipo<-stack("chl")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
pipo<-crop(pipo,extent(-0.6,0.5,49.2,49.7))
tschl<-cellStats(pipo,mean)
pipo<-stack("sst")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
pipo<-crop(pipo,extent(-0.6,0.5,49.2,49.7))
tssst<-cellStats(pipo,mean)
pipo<-stack("kd490")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
pipo<-crop(pipo,extent(-0.6,0.5,49.2,49.7))
tskd490<-cellStats(pipo,mean)
pipo<-stack("spm")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
pipo<-crop(pipo,extent(-0.6,0.5,49.2,49.7))
tsspm<-cellStats(pipo,mean)

datenv<-data.frame(t(rbind(tschl[4:21],c(tssst[6:20],NA,tssst[21:22]),tskd490[4:21],tsspm[4:21])))
names(datenv)<-c("chl","sst","kd490","spm")
#debit
load("../debit/datdeb.rdata")


datall<-cbind(data.frame(w=datpech$w,h=datpech$tps),datenv,disw=datdeb[6:23,4],year=2000:2017)
datall$sst[16]<-c(13.57993+13.03493)/2

#
plot(datall)
library(psych)
pairs.panels(datall,method="spearman",density=T,ellipses=T)


mod1<-(lm(w/h~chl+sst+kd490+spm+disw,data=datall))

mod1<-(lm(w~h+chl+sst+kd490+spm+disw,data=datall))
mod2<-(lm(w~h+sst+kd490+spm+disw,data=datall))

mod1<-lm(w~.,data=datall)
anova(mod1)
anova(lm(w~h+chl+sst+kd490+spm+disw,data=datall))
anova(lm(w~h+chl+kd490+spm+disw,data=datall))
anova(lm(w~h+kd490+spm+disw,data=datall))
anova(lm(w~h+spm,data=datall))
m<-lm(w~h+spm,data=datall)

summary(m)
confint(m,level=0.95)
influence(m)
par(mfrow = c(2,2))
plot(m)
par(mfrow = c(1,1))
plot(predict(m)~datall$w)
pipo <- data.frame(w= datall$w/datall$h, wpred = predict(m), year = 2000:2017)
ggplot(data = pipo, aes(y=wpred, x = w)) + geom_point() + geom_smooth(method="lm") + labs(x = "Débarquements observés", y = "Débarquements prédits") 





plot(mean(pipo))
plot(coast3,add=T)
plot(stackApply(pipo,rep(1,21),fun=mean))
chlstrata<-extract(chl,maprect,method="simple")#,fun=mean,na.rm=T)
chlstrata<-extract(chl,maprect,method="bilinear",fun=mean,na.rm=T)
chlstrata<-extract(chl,data.frame(x=0.2,y=49.45),buffer=30,fun=mean,na.rm=T)




chlstrata<-cbind(maprect$SECT_COD,chlstrata)
datchl<-tidyr::gather(chlstrata,"year","chl",3:ncol(chlstrata))
#kp490
pipo<-stack("kd490")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
levelplot(pipo,zscale=T,contour=T)
pipo<-extract(pipo,mapstrate,method="bilinear",fun=mean,na.rm=T)
pipo<-cbind(mapstrate@data,pipo)
datkd490<-tidyr::gather(pipo,"year","kd490",3:ncol(pipo))
#bbp443
pipo<-stack("bbp443")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
levelplot(pipo,zscale=T,contour=T)
pipo<-extract(pipo,mapstrate,method="bilinear",fun=mean,na.rm=T)
pipo<-cbind(mapstrate@data,pipo)
datbbp443<-tidyr::gather(pipo,"year","bbp443",3:ncol(pipo))
#cdm443
pipo<-stack("cdm443")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
levelplot(pipo,zscale=T,contour=T)
pipo<-extract(pipo,mapstrate,method="bilinear",fun=mean,na.rm=T)
pipo<-cbind(mapstrate@data,pipo)
datcdm443<-tidyr::gather(pipo,"year","cdm443",3:ncol(pipo))
#spm
pipo<-stack("spm")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
levelplot(pipo,zscale=T,contour=T)
pipo<-extract(pipo,mapstrate,method="bilinear",fun=mean,na.rm=T)
pipo<-cbind(mapstrate@data,pipo)
datspm<-tidyr::gather(pipo,"year","spm",3:ncol(pipo))
#sst
pipo<-stack("sst")
pipo<-stackApply(pipo,substr(names(pipo),7,10),mean,na.rm=T)
levelplot(pipo)
levelplot(pipo,zscale=T,contour=T)
pipo<-extract(pipo,mapstrate,method="bilinear",fun=mean,na.rm=T)
pipo<-cbind(mapstrate@data,pipo)
datsst<-tidyr::gather(pipo,"year","sst",3:ncol(pipo))

#datsat
datsat<-datsst%>%full_join(datchl)%>%full_join(datkd490)%>%
	full_join(datcdm443)%>%full_join(datbbp443)%>%full_join(datspm)%>%
	mutate(year=gsub("index_","",year))
save(datsat,file="datsat.rdata")
load("datsat.rdata")

#verif
pipo<-tidyr::gather(datsat,"var","val",4:9)
ggplot(pipo,aes(x=year,y=val,group=Ifremer_id,color=Ifremer_id))+
	geom_path()+facet_wrap(~var,scale="free")






