library(Momocs)
library(geomorph)
#Repeatability and measurement error
MTFTWCoords <- list.files("D:\\Coordinates\\Free Threshing\\FTW _Method_Test", full.names = TRUE)
MTFTWFrame <- read.csv("D:\\Coordinate_Frames\\MethodTest_FTW.csv", header = TRUE)
MTFTWTxt <- import_txt(MTFTWCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
MTFTWOut <- Out(MTFTWTxt, fac=MTFTWFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
MTFTWOut.l <- filter(MTFTWOut, View == "l")#creates subset of lateral views only
MTFTWOut.d <- filter(MTFTWOut, View == "d")#creates subset of lateral views only
MTFTWOut.l.efour <- efourier(MTFTWOut.l, nb.h=8, norm = FALSE, start = FALSE)
MTFTWOut.d.efour <- efourier(MTFTWOut.d, nb.h=8, norm = FALSE, start = FALSE)
MTFTWOut.d.l <- combine (MTFTWOut.d, MTFTWOut.l)
#Procrustes ANOVA
export(MTFTWOut.l.efour)
export(MTFTWOut.d.efour)
MTFTW_l <- read.csv("D:\\Coordinate_Frames\\MTFTW_l.csv")
MTFTW_d <- read.csv("D:\\Coordinate_Frames\\MTFTW_d.csv")
#Calculate PCA scores for above
PCAObservation_d<-prcomp(MTFTW_d)$x 
PCAObservation_l<-prcomp(MTFTW_l)$x 
#Creates 5 level session factor (consecutive)
sessionfactor<-as.factor(gl(5, 5))
sessionfactor
#Create 5 level individual factor (repeated)
individualfactor<-as.factor(rep(1:5, 5))
individualfactor
#Geomorph data frames for PCA computations of coefficients
gdf_d<-geomorph.data.frame(shape=PCAObservation_d) #ID: a vector containing the specimens ID 
gdf_l<-geomorph.data.frame(shape=PCAObservation_l) #ID: a vector containing the specimens ID 
#Procrustes ANOVA using session factor and then individual factor as dependent variables
summary(procD.lm(shape~sessionfactor, data=gdf_d))
summary(procD.lm(shape~sessionfactor, data=gdf_l))
#=session difference not significant for dorsal or lateral view
mod<-summary(procD.lm(shape~individualfactor, data= gdf_d))
mod2<-summary(procD.lm(shape~individualfactor, data= gdf_l))
mod
mod2
#=indiv difference significant for dorsal and lateral view
#dorsal measurement error (See Claude 2008)
s2within<-mswithin<-mod[[1]][2,3]
mod[[1]][2, 3]
MSamong<-mod[[1]][1, 3]
MSamong
s2among<-(MSamong-mswithin)/5
s2within/(s2within+s2among)*100#measurement error
#lateral measurement error
s2within<-mswithin<-mod2[[1]][2,3]
mod2[[1]][2, 3]
MSamong<-mod2[[1]][1, 3]
MSamong
s2among<-(MSamong-mswithin)/5
s2within/(s2within+s2among)*100#measurement error

#Stage 1 Analysis 
NeoRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\CharredNeo_Ref", full.names = TRUE)
NeoRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Neo_Ref.csv", header = TRUE)
NeoRefTxt <- import_txt(NeoRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut <- Out(NeoRefTxt, fac=NeoRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2 <-coo_scale (NeoRefOut)
#lateral
NeoRefOut.l <- filter(NeoRefOut2, View == "l")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoRefOut.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoRefOut.l.efour <- efourier(NeoRefOut.l, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut.d <- filter(NeoRefOut2, View == "d")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoRefOut.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoRefOut.d.efour <- efourier(NeoRefOut.d, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut.d.l <- combine (NeoRefOut.d, NeoRefOut.l)#dataset of just dorsal and lateral views
#?LDA of species
NeoRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoRefOut.d.l.lda <- NeoRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(NeoRefOut.d.l.lda, col=c("#24947A", "#0000FF", "#c44918", "#AA29E1"), points= TRUE , labelspoints= FALSE, zoom= 1.25, cex.labelsgroups= 0.95, rect.labelsgroups= TRUE, cex = 0.1)

#LDA of accession
NeoRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('Accession', scale=FALSE, center= TRUE) 
NeoRefOutA.d.l.lda <- NeoRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('Accession', scale=FALSE, center= TRUE)

###############Landrace only
NeoRefCoords_LR <- list.files("D:\\Coordinates\\Free Threshing\\CharredNeo_Ref_LR", full.names = TRUE)
NeoRefFrame_LR <- read.csv("D:\\Coordinate_Frames\\Charred_Neo_Ref_LR.csv", header = TRUE)
NeoRefTxt_LR <- import_txt(NeoRefCoords_LR, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut_LR <- Out(NeoRefTxt_LR, fac=NeoRefFrame_LR) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2_LR <-coo_scale (NeoRefOut_LR)
#lateral
NeoRefOut_LR.l <- filter(NeoRefOut2_LR, View == "l")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoRefOut_LR.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoRefOut_LR.l.efour <- efourier(NeoRefOut_LR.l, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut_LR.d <- filter(NeoRefOut2_LR, View == "d")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoRefOut_LR.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoRefOut_LR.d.efour <- efourier(NeoRefOut_LR.d, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut_LR.d.l <- combine (NeoRefOut_LR.d, NeoRefOut_LR.l)#dataset of just dorsal and lateral views
#?LDA
NeoRefOut_LR.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoRefOut_LR.d.l.lda <- NeoRefOut_LR.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
NeoRefOut_LR.d.l.lda
plot(NeoRefOut_LR.d.l.lda, col=c("#24947A", "#0000FF", "#629F53", "#AA29E1"), points= TRUE , labelspoints= FALSE, zoom= 1, cex.labelsgroups= 0.95, rect.labelsgroups= TRUE, cex = 0.1)
###

#Stage 2.1 comparison with archaeobotanical material (all material)
NeoCombCoords <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo", full.names = TRUE)
NeoCombFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Combined_Neo.csv", header = TRUE)
NeoCombTxt <- import_txt(NeoCombCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoCombOut <- Out(NeoCombTxt, fac=NeoCombFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoCombOut2 <-coo_scale (NeoCombOut)
#lateral
NeoCombOut.l <- filter(NeoCombOut2, View == "l")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoCombOut.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoCombOut.l.efour <- efourier(NeoCombOut.l, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoCombOut.d <- filter(NeoCombOut2, View == "d")#creates subset of lateral views only
NeoCombOut.d.efour <- efourier(NeoCombOut.d, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoCombOut.d.l <- combine (NeoCombOut.d, NeoCombOut.l)#dataset of just dorsal and lateral views
#?LDA
NeoCombOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoCombOut.d.l.lda <- NeoCombOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(NeoCombOut.d.l.lda, col=c("#f018cf","#fca605","#24947A", "#0000FF", "#c44918", "#AA29E1"), points= TRUE , labelspoints= FALSE, zoom= 1.35, cex.labelsgroups= 0.95, rect.labelsgroups= TRUE, cex = 0.1)

###
#Stage 2.2 comparison with archaeobotanical material - Land races
#Catal
NeoCombCoords_LRCatal <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_LR_Catal", full.names = TRUE)
NeoCombFrame_LRCatal <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Combined_Neo_LR_Catal.csv", header = TRUE)
NeoCombTxt_LRCatal <- import_txt(NeoCombCoords_LRCatal, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoCombOut_LRCatal <- Out(NeoCombTxt_LRCatal, fac=NeoCombFrame_LRCatal) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoCombOut2_LRCatal <-coo_scale (NeoCombOut_LRCatal)
#lateral
NeoCombOut_LRCatal.l <- filter(NeoCombOut2_LRCatal, View == "l")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoCombOut_LRCatal.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoCombOut_LRCatal.l.efour <- efourier(NeoCombOut_LRCatal.l, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoCombOut_LRCatal.d <- filter(NeoCombOut2_LRCatal, View == "d")#creates subset of lateral views only
NeoCombOut_LRCatal.d.efour <- efourier(NeoCombOut_LRCatal.d, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsalLDA
NeoCombOut_LRCatal.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoCombOut_LRCatal.d.lda <- NeoRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#combining just dorsal and lateral
NeoCombOut_LRCatal.d.l <- combine (NeoCombOut_LRCatal.d, NeoCombOut_LRCatal.l)#dataset of just dorsal and lateral views
#?LDA
NeoCombOut_LRCatal.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('Accession', scale=FALSE, center= TRUE) 
NeoCombOut_LRCatal.d.l.lda <- NeoCombOut_LRCatal.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(NeoCombOut_LRCatal.d.l.lda, col=c("#f018cf","#24947A", "#0000FF", "#c44918", "#AA29E1"), points= TRUE , center= TRUE, labelspoints= FALSE, zoom= 1.8, cex.labelsgroups= 0.90, rect.labelsgroups= TRUE, cex = 0.1, labelsgroups= TRUE)
?plot.LDA

#same but with Kouph
NeoCombCoords_LRKouph <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_LR_Kouph", full.names = TRUE)
NeoCombFrame_LRKouph <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Combined_Neo_LR_Kouph.csv", header = TRUE)
NeoCombTxt_LRKouph <- import_txt(NeoCombCoords_LRKouph, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoCombOut_LRKouph <- Out(NeoCombTxt_LRKouph, fac=NeoCombFrame_LRKouph) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoCombOut2_LRKouph <-coo_scale (NeoCombOut_LRKouph)
#lateral
NeoCombOut_LRKouph.l <- filter(NeoCombOut2_LRKouph, View == "l")#creates subset of lateral views only
calibrate_harmonicpower_efourier(NeoCombOut_LRKouph.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
NeoCombOut_LRKouph.l.efour <- efourier(NeoCombOut_LRKouph.l, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoCombOut_LRKouph.d <- filter(NeoCombOut2_LRKouph, View == "d")#creates subset of lateral views only
NeoCombOut_LRKouph.d.efour <- efourier(NeoCombOut_LRKouph.d, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsalLDA
NeoCombOut_LRKouph.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoCombOut_LRKouph.d.lda <- NeoRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#combining just dorsal and lateral
NeoCombOut_LRKouph.d.l <- combine (NeoCombOut_LRKouph.d, NeoCombOut_LRKouph.l)#dataset of just dorsal and lateral views
#?LDA
NeoCombOut_LRKouph.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoCombOut_LRKouph.d.l.lda <- NeoCombOut_LRKouph.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(NeoCombOut_LRKouph.d.l.lda, col=c("#f018cf","#24947A", "#0000FF", "#c44918", "#AA29E1" ), points= TRUE , labelspoints= FALSE, zoom= 1.5, cex.labelsgroups= 0.90, rect.labelsgroups= TRUE, cex = 0.1, labelsgroups= TRUE)

#Reclassification
##unclass Catal
UnclassCoords <- list.files("D:\\Coordinates\\Free Threshing\\Unclass_Catal", full.names = TRUE)
UnclassFrame <- read.csv("D:\\Coordinate_Frames\\Unclass_Catal.csv", header = TRUE)
UnclassTxt <- import_txt(UnclassCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
UnclassOut <- Out(UnclassTxt, UnclassFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
UnclassOut1 <-coo_scale(UnclassOut)
#reclass
UnclassOut.d <- filter(UnclassOut1, View == "d")
UnclassOut.l <- filter(UnclassOut1, View == "l")
UnclassOut.d.l <- combine (UnclassOut.d, UnclassOut.l)
UnclassOut.d.l.efour <- UnclassOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
reLDA(UnclassOut.d.l.efour, NeoRefOut_LR.d.l.lda)
NeoRefOut_LR.d.l.lda
#no priors
mod<-NeoRefOut_LR.d.l.lda$mod
newdata<-UnclassOut.d.l.efour
predicted<-predict(mod, newdata$coe, prior= c(1,1,1,1)/4)
predicted

#Unclass Kouph
UnclassCoordsK <- list.files("D:\\Coordinates\\Free Threshing\\Unclass_Kouph", full.names = TRUE)
UnclassFrameK <- read.csv("D:\\Coordinate_Frames\\Unclass_Kouph.csv", header = TRUE)
UnclassTxtK <- import_txt(UnclassCoordsK, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
UnclassOutK <- Out(UnclassTxtK, UnclassFrameK) #creates an Out object from a specified list- you can specify landmarks in "ldk"
UnclassOut1K <-coo_scale(UnclassOutK)
#reclass
UnclassOut.dK <- filter(UnclassOut1K, View == "d")
UnclassOut.lK <- filter(UnclassOut1K, View == "l")
UnclassOut.d.lK <- combine (UnclassOut.dK, UnclassOut.lK)
UnclassOut.d.l.efourK <- UnclassOut.d.lK %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
reLDA(UnclassOut.d.l.efourK, NeoRefOut_LR.d.l.lda)
#no priors
mod<-NeoRefOut_LR.d.l.lda$mod
newdata<-UnclassOut.d.l.efourK
predicted<-predict(mod, newdata$coe, prior= c(1,1,1,1)/4)
predicted


#pairwise comparisons 
NeoRefCoords_PW <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_LR", full.names = TRUE)
NeoRefFrame_PW <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Combined_Neo_LR.csv", header = TRUE)
NeoRefTxt_PW <- import_txt(NeoRefCoords_PW, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut_PW <- Out(NeoRefTxt_PW, fac=NeoRefFrame_PW) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2_PW <-coo_scale (NeoRefOut_PW)
#lateral
NeoRefOut.l_PW <- filter(NeoRefOut2_PW, View == "l")#creates subset of lateral views only
NeoRefOut.l.efour_PW <- efourier(NeoRefOut.l_PW, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut.d_PW <- filter(NeoRefOut2_PW, View == "d")#creates subset of lateral views only
NeoRefOut.d.efour_PW <- efourier(NeoRefOut.d_PW, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut.d.l_PW <- combine (NeoRefOut.d_PW, NeoRefOut.l_PW)#dataset of just dorsal and lateral views

#MANOVA/Mean shapes
MeanDorsal<-MSHAPES(NeoRefOut.l.efour_PW, fac = 'taxon.code', FUN = mean, nb.pts = 120)
MeanDorsal2<-MSHAPES(NeoRefOut.d.efour_PW, fac = 'taxon.code', FUN = mean, nb.pts = 120)
Neo.ms <- MSHAPES(NeoRefOut.d.efour_PW, 'taxon.code')
Neo.ms2<- MSHAPES(NeoRefOut.l.efour_PW, 'taxon.code')
#pairwise comparison
plot_mshapes (Neo.ms, size= 3/4)
plot_mshapes (Neo.ms2, size= 3/4)

#Stage 3: Varietal comparisons 
#Catal


###
NeoRefCoords_var_C <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_Hex", full.names = TRUE)
NeoRefFrame_var_C <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Combined_Neo_Hex.csv", header = TRUE)
NeoRefTxt_var_C <- import_txt(NeoRefCoords_var_C, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut_var_C <- Out(NeoRefTxt_var_C, fac=NeoRefFrame_var_C) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2_var_C <-coo_scale (NeoRefOut_var_C)
#lateral
NeoRefOut.l_var_C <- filter(NeoRefOut2_var_C, View == "l")#creates subset of lateral views only
NeoRefOut.l.efour_var_C <- efourier(NeoRefOut.l_var_C, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut.d_var_C <- filter(NeoRefOut2_var_C, View == "d")#creates subset of lateral views only
NeoRefOut.d.efour_var_C <- efourier(NeoRefOut.d_var_C, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut.d.l_var_C <- combine (NeoRefOut.d_var_C, NeoRefOut.l_var_C)#dataset of just dorsal and lateral views
#?LDA
NeoRefOut.d.l_var_C %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoRefOut.d.l.lda_var_C <- NeoRefOut.d.l_var_C %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(NeoRefOut.d.l.lda_var_C, col=c("#24947A", "#0000FF", "#629F53", "#AA29E1", "#584848", "#CE8DA3", "#8A8989", "#F56505", "#3D0AA2"), points= TRUE , labelspoints= FALSE, zoom= 1.25, cex.labelsgroups= 0.95, rect.labelsgroups= TRUE, cex = 0.1)

#Reclass (same as above but without arch material)
NeoRefCoords_var_C <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_Hex2", full.names = TRUE)
NeoRefFrame_var_C <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Combined_Neo_Hex2.csv", header = TRUE)
NeoRefTxt_var_C <- import_txt(NeoRefCoords_var_C, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut_var_C <- Out(NeoRefTxt_var_C, fac=NeoRefFrame_var_C) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2_var_C <-coo_scale (NeoRefOut_var_C)
#lateral
NeoRefOut.l_var_C <- filter(NeoRefOut2_var_C, View == "l")#creates subset of lateral views only
NeoRefOut.l.efour_var_C <- efourier(NeoRefOut.l_var_C, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut.d_var_C <- filter(NeoRefOut2_var_C, View == "d")#creates subset of lateral views only
NeoRefOut.d.efour_var_C <- efourier(NeoRefOut.d_var_C, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut.d.l_var_C <- combine (NeoRefOut.d_var_C, NeoRefOut.l_var_C)#dataset of just dorsal and lateral views
#?LDA
NeoRefOut.d.l_var_C %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoRefOut.d.l.lda_var_C <- NeoRefOut.d.l_var_C %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#reclass
reLDA(UnclassOut.d.l.efour, NeoRefOut.d.l.lda_var_C)

#Varietal comparisons: Kouphovouno
NeoRefCoords_var_K <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_Kouph_Dur", full.names = TRUE)
NeoRefFrame_var_K <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Kouph_Dur.csv", header = TRUE)
NeoRefTxt_var_K <- import_txt(NeoRefCoords_var_K, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut_var_K <- Out(NeoRefTxt_var_K, fac=NeoRefFrame_var_K) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2_var_K <-coo_scale (NeoRefOut_var_K)
#lateral
NeoRefOut.l_var_K <- filter(NeoRefOut2_var_K, View == "l")#creates subset of lateral views only
NeoRefOut.l.efour_var_K <- efourier(NeoRefOut.l_var_K, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut.d_var_K <- filter(NeoRefOut2_var_K, View == "d")#creates subset of lateral views only
NeoRefOut.d.efour_var_K <- efourier(NeoRefOut.d_var_K, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut.d.l_var_K <- combine (NeoRefOut.d_var_K, NeoRefOut.l_var_K)#dataset of just dorsal and lateral views
#?LDA
NeoRefOut.d.l_var_K %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoRefOut.d.l.lda_var_K <- NeoRefOut.d.l_var_K %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(NeoRefOut.d.l.lda_var_K, col=c("#24947A", "#0000FF", "#629F53", "#AA29E1", "#584848", "#CE8DA3", "#8A8989"), points= TRUE , labelspoints= FALSE, zoom= 1.25, cex.labelsgroups= 0.95, rect.labelsgroups= TRUE, cex = 0.1)

#Reclass (same as above but without arch material)
NeoRefCoords_var_K <- list.files("D:\\Coordinates\\Free Threshing\\Combined_Neo_Kouph_Dur2", full.names = TRUE)
NeoRefFrame_var_K <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT_Kouph_Dur2.csv", header = TRUE)
NeoRefTxt_var_K <- import_txt(NeoRefCoords_var_K, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
NeoRefOut_var_K <- Out(NeoRefTxt_var_K, fac=NeoRefFrame_var_K) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NeoRefOut2_var_K <-coo_scale (NeoRefOut_var_K)
#lateral
NeoRefOut.l_var_K <- filter(NeoRefOut2_var_K, View == "l")#creates subset of lateral views only
NeoRefOut.l.efour_var_K <- efourier(NeoRefOut.l_var_K, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#dorsal
NeoRefOut.d_var_K <- filter(NeoRefOut2_var_K, View == "d")#creates subset of lateral views only
NeoRefOut.d.efour_var_K <- efourier(NeoRefOut.d_var_K, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
#combining just dorsal and lateral
NeoRefOut.d.l_var_K <- combine (NeoRefOut.d_var_K, NeoRefOut.l_var_K)#dataset of just dorsal and lateral views
#?LDA
NeoRefOut.d.l_var_K %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
NeoRefOut.d.l.lda_var_K <- NeoRefOut.d.l_var_K %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
reLDA(UnclassOut.d.l.efourK, NeoRefOut.d.l.lda_var_K)
mod<-NeoRefOut.d.l.lda_var_K$mod
newdata<-UnclassOut.d.l.efourK
predicted<-predict(mod, newdata$coe, prior= c(1,1,1,1,1)/5)
predicted


#3d LDA

#3d plot for Part 1
library(rgl)
gcolours<-c("#24947A", "#0000FF", "#c44918", "#AA29E1")
open3d()
par3d(windowRect = c(100, 100, 612, 612))
NeoRefOut.d.l.lda.colour<-gcolours[as.factor(NeoRefOut.d.l.lda$fac)]
LD<-as.data.frame(NeoRefOut.d.l.lda[["mod.pred"]][["x"]])
LD
plot3d(LD$LD1, LD$LD2, LD$LD3,col=NeoRefOut.d.l.lda.colour, xlab= "LD1", ylab="LD2", zlab="LD3")
rgl.snapshot("modernaccessions.png","png")

#3d plot comparing landrace data and Catal 
gcolours<-c("#f018cf","#24947A", "#0000FF", "#c44918", "#AA29E1")
open3d()
par3d(windowRect = c(100, 100, 612, 612))
NeoCombOut_LRCatal.d.l.lda.colour<-gcolours[as.factor(NeoCombOut_LRCatal.d.l.lda$fac)]
LD<-as.data.frame(NeoCombOut_LRCatal.d.l.lda[["mod.pred"]][["x"]])
LD
plot3d(LD$LD1, LD$LD2, LD$LD3,col=NeoCombOut_LRCatal.d.l.lda.colour, xlab= "LD1", ylab="LD2", zlab="LD3")

NeoCombOut_LRCatal.d.l.lda
