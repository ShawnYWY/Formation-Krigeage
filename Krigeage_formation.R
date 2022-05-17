
library(raster)
library(gstat)
library(fields)
library(sp)
library(rgdal)
library(proj4)

######################
# LECTURE DES DONNEES
######################
setwd('/mnt/dcappli/DEV_PERSO/WYa/tuto_krigeage/Formation-Krigeage/')
data = read.table(file='Troyes/data.csv',h=T)
coordinates(data) = ~X+Y 
proj4string(data) = CRS("+init=epsg:27572")
urbain=data[which(data$Typo=="U"),]

grid = read.table(file='Troyes/grid.csv',h=T)
coordinates(grid) = ~X+Y
proj4string(grid) = CRS("+init=epsg:27572")

# CREATION DU VARIOGRAMME EXPERIMENTAL
vexp=variogram(NO2~1,urbain)

# CREATION DU VARIOGRAMME THEORIQUE
v0=vgm(model='Sph',nugget=2,psill=15,range=3000) 	# Variogramme d initilasation
vfit=fit.variogram(vexp,model=v0)				# Iteratif

plot(vexp,vfit)

# KRIGEAGE ORDINAIRE
objg=gstat(id='obs',formula=NO2~1,data=urbain,model=vfit) #DEFINIT OBJET Gstat

#pt = data.frame(X=729500,Y=2367500)	#pt d export
#coordinates(pt)=~+ X+Y
#proj4string(pt) = CRS("+init=epsg:27572")

predict.gstat(object=objg,pt)		#Prediction au point
krig=predict.gstat(object=objg,grid)	#Prediction sur grille pas de vairaible auxiliaire grace a NO2~1

#CARTE
grille_fin=t(matrix(krig$obs.pred,nrow=length(unique(grid$y)),ncol=length(unique(grid$y))))
grille_fin=matrix(krig$obs.pred,byrow=T,nrow=length(unique(grid$x)),ncol=length(unique(grid$y)))

image.plot(grille_fin)				# sutilise avec une matrice 
spplot(krig)				# utilise avec des spatial point 

# VALIDATION CROISEE
vc=gstat.cv(objg,nfold=nrow(urbain))
plot(vc@data$obs.pred,vc@data$observed)
res=(vc@data$obs.pred)-(vc@data$observed)
hist(res)

#####################
# DERIVE EXTERNE
#####################
# REFAIRE Krigeage en fonction emission NOx
# DERIVE EXTERNE DE LA VARIABLE EXPLICATIVE NOX
# On suppose que NO2 depend des emissions de NOX de maniere lineaire
# On estime pour chaque point la regression lineaire entre NOX et NO2
# On cree le variogramme des residus 
# On applique le krigeage des residus pour obtenir une grille de residu 
# On 
lin=lm(NO2~NOx,data=urbain)


#FONCTION iterative : estime la regression lineaire entre NO2 et NOx et renvoi residus de la regression
reg_VC = function(formula, df) {
	res=NULL
	for (i in 1:nrow(df)) {
		reg=lm(formula, data=df[-i,])
		restmp=predict(reg,newdata=df[i,])-df$NO2[i]
		res=c(res, restmp)
	}
	return (res)
}

res=reg_VC(formula(lin),urbain)
residu=cbind(res,urbain$NOx)
residu=data.frame(residu)
names(residu)=c("res","NOx")
coordinates(residu)=coordinates(urbain)
proj4string(residu)=CRS("+init=epsg:27572")


# CREATION DU VARIOGRAMME EXPERIMENTAL derive externe
vdexp=variogram(res~1,data=residu)

# CREATION DU VARIOGRAMME THEORIQUE derive externe
vde0=vgm(model='Sph',nugget=0,psill=3.5,range=4000) 	# Variogramme d initilasation
vdefit=fit.variogram(vdexp,model=vde0)				# Iteratif

objgde=gstat(id='obs',formula=NO2~NOx,data=urbain,model=vdefit) #DEFINIT OBJET Gstat
krigde=predict.gstat(object=objgde,grid)	#Prediction sur grille pas de vairaible auxiliaire grace a NO2~1

# VALIDATION CROISEE
vcde=gstat.cv(objgde,nfold=nrow(residu))
plot(vcde@data$obs.pred,vcde@data$observed)
res=(vcde@data$obs.pred)-(vcde@data$observed)
hist(res)

######################
# KRIGEAGE DES RESIDUS
######################

lin=lm(NO2~NOx,data=urbain)
a=lin$coefficients[1]
b=lin$coefficients[2]

prevlin=lin$fitted.values
obs=urbain$NO2

reslin=prevlin-obs
urbain$reslin=reslin

# on va kriger les residus : krigeage ordinaire

# CREATION DU VARIOGRAMME EXPERIMENTAL des residus
vexp=variogram(reslin~1,urbain)

# variogramme theo
v0=vgm(model="Sph",range=4000, nugget=0,psill=3.5)
vres=fit.variogram(vexp,v0)

# krigeage des residus
objres=gstat(id='obs',formula=reslin~1,data=urbain,model=vres) #DEFINIT OBJET Gstat
krigres=predict.gstat(object=objres,grid)	#Prediction sur grille pas de vairaible auxiliaire grace a NO2~1


# Sommme a*grille_emission$NOX+b+krigres = NO2 a partir des emissions NOX en ajoutant les residu
no2pred=predict(lin,grid)+krigres@data$obs.pred
gridfinres=a*grid$NOx+b+krigres@data$obs.pred
gridfinres=matrix(gridfinres,byrow=T,nrow=length(unique(coordinates(grid)[,1])),ncol=length(unique(coordinates(grid)[,2])))
spplot(grid,"gridfinres")


################################
# MANIPULATION DES CARTES 
################################
shape=readOGR(dsn='C:/Users/LibreService/Desktop/FormationR/Troyes/shapefiles',layer='Troyes')	#LECTURE DU SHAPEFILE
contour = spTransform(shape,CRS("+init=epsg:27572"))

X = unique(coordinates(grid)[,1])
Y = unique(coordinates(grid)[,2])
rst_no2 =raster(list(x=X,y=Y,z=gridfinres))
proj4string(rst_no2)=proj4string(contour)
rst_no2_10m = disaggregate(rst_no2,fact=c(2),method='bilinear') # RESOLUTION A 10m
rst_no2_10m_decoupe=rasterize(x=contour,y=rst_no2_10m,mask=T)	#DECOUPAGE EN FONCTION DU CONTOUR commune	

###############################
# test de la fonction interp.surface pour COMBINE SIRANE + PREVALP
###############################
#LECTURE DONNEES CHIMERE
library(ncdf)
nc=open.ncdf('C:/Users/LibreService/Desktop/FormationR/DONNEE_JOUJOU/stat_ana_avg24.20130101_20131231_REG01KM.nc')
POL=get.var.ncdf(nc,varid='pm10_moy_an')
LAT=get.var.ncdf(nc,varid='lat')
LON=get.var.ncdf(nc,varid='lon')
close.ncdf(nc)
polluant=data.frame(x=as.vector(LON),y=as.vector(LAT),data=as.vector(POL))
coordinates(polluant)=~x+y
proj4string(polluant)=CRS("+init=epsg:4326")
polluant=spTransform(polluant,CRS("+init=epsg:32631"))
X=unique(coordinates(polluant)[,1])
Y=unique(coordinates(polluant)[,2])
Z=polluant@data$data

#LECTURE DONNEE SIRANE
rst_sirane=raster('Grille_PM_Moy.CUGDLyon.20130303.grd')
X=unique(coordinates(rst_sirane)[,1])
Y=unique(coordinates(rst_sirane)[,2])


#INTERPOLATION CHIMERE --> grille SIRANE
prevalp_10m=interp.surface()



