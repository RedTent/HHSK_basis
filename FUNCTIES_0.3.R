#library(lubridate)
#library(trend)
#library(dplyr)
#library(sp)
#library(ggplot2)
#library(gridExtra)
#library(sp)
#library(maptools)
#library(gstat)
#library(raster)
#library(microbenchmark)

# GEGEVENS IMPORTEREN ----------------------------------------------------------------------

import_data <- function(datacsv="fys_chem.csv"){
  require(readr)
  require(dplyr)
  require(lubridate)
  #kolomnamen <- c("mp","datum","parnr","par","eenheid","detectiegrens","waarde")
  df <- read_csv2(file=datacsv,col_types=cols(datum=col_date(format="%d-%m-%Y %H:%M:%S")))
  df <- dplyr::filter(df,!is.na(waarde)) # alle metingen moeten een meetwaarde hebben
  #Toevoegen jaren en maanden
  df$jaar <-as.integer(year(df$datum))
  df$maand <- as.integer(month(df$datum))

  #info zodat je weet wat je importeert
  print(paste("Laatste meetdatum is",max(df$datum)))
  df
}

import_meetpunten <- function(meetpuntencsv="meetpunten.csv"){
  require(dplyr)
  require(readr)
  meetpuntendf <- read_csv2(meetpuntencsv,col_types = cols())
  names(meetpuntendf) <- tolower(names(meetpuntendf))
  meetpuntendf <- meetpuntendf %>% rename(X=x,Y=y)
  meetpuntendf
}

import_meetpunten_latlong <- function(meetpuntencsv="meetpunten.csv", X = "X", Y = "Y"){
  require(rgdal)
  require(dplyr)
  meetpuntendf <- import_meetpunten(meetpuntencsv)
  longlat <- meetpuntendf %>% filter(X != 0, Y != 0) %>%  mutate(long = X, lat = Y)
  coordinates(longlat) = ~long+lat
  proj4string(longlat) <- CRS("+init=EPSG:28992")
  longlat <- spTransform(longlat,"+init=EPSG:4326")
  meetpuntendf <- left_join(meetpuntendf, select(as_data_frame(longlat), mp, long, lat), by = "mp")
}

import_parameters <- function(parametercsv='parameters.csv'){
  require(readr)
  parameterdf <-read_csv2(parametercsv, col_types = cols())
  parameterdf
}

import_bio_stadia <- function(biologiecsv = "biologie.csv"){
  require(readr)
  require(lubridate)
  require(dplyr)
  biodf <- read_csv2("biologie.csv",col_types=cols(datum=col_date(format="%d-%m-%Y %H:%M:%S")))
  biodf$monsterid <- group_indices(biodf, mp,datum,methode)
  biodf$jaar <- as.integer(year(biodf$datum))
  biodf
}

import_bio <- function(biologiecsv = "biologie.csv"){
  require(readr)
  require(lubridate)
  require(dplyr)
  biodf <- read_csv2(biologiecsv,col_types=cols(datum=col_date(format="%d-%m-%Y %H:%M:%S")))
  biodf$monsterid <- group_indices(biodf, mp,datum,methode)
  biodf$jaar <- as.integer(year(biodf$datum))
  biodf <- biodf %>% select(-stadium,-stadiumwaarde) %>% unique()
  biodf
}

#biodf <- read_csv2("biologie.csv",col_types=cols(datum=col_date(format="%d-%m-%Y %H:%M:%S")))

#biodf$monsterid <- group_indices(biodf, mp,datum,methode)
#biodf$jaar <- as.integer(year(biodf$datum))

#biostadia <- biodf
#biodf <- biodf %>% select(-stadium,-stadiumwaarde) %>% unique() %>% filter(jaar<2017)


import_bio_kenmerken <- function(bio_kenmerkencsv="biologie_kenmerken.csv"){
  require(readr)
  bio_km_df <- read_csv2(bio_kenmerkencsv,col_types = cols())
  bio_km_df
}

# INFORMATIE ----------------------------------------------------------------------

meetpuntinformatie <- function(mpcode="00016", meetpuntendf){
    mpomsch <- meetpuntendf[meetpuntendf[["mp"]]==mpcode,"mpomsch"]
    X <- meetpuntendf[meetpuntendf[["mp"]]==mpcode,"X"]
    Y <- meetpuntendf[meetpuntendf[["mp"]]==mpcode,"Y"]
    gebied <- meetpuntendf[meetpuntendf[["mp"]]==mpcode,"gebied"]
    meetnet <- meetpuntendf[meetpuntendf[["mp"]]==mpcode,"meetnet"]
    meetpuntsoort <- meetpuntendf[meetpuntendf[["mp"]]==mpcode,"meetpuntsoort"]
    return(list("mpcode"=mpcode,"mpomsch"=mpomsch,"X"=X,"Y"=Y,"gebied"=gebied,"meetnet"=meetnet,"meetpuntsoort"=meetpuntsoort))
}#end of function

parameterinformatie <- function(parnr="3",parameterdf){
    parameterkort <- as.character(parameterdf[parameterdf[["parnr"]]==parnr,"par"])
    parameterlang <-  as.character(parameterdf[parameterdf[["parnr"]]==parnr,"parnaamlang"])
    eenheid <-  as.character(parameterdf[parameterdf[["parnr"]]==parnr,"eenheid"])
    return(list("parnr"=parnr,"parameterkort"=parameterkort,"parameterlang"=parameterlang,"eenheid"=eenheid))
}#end of function

# ANALYSE ----------------------------------------------------------------------

# trendtests MK en Theil-Sen, zowel gewoon als seasonal, tegelijk!
trendanalysemk <- function(df,mpkolom="mp",parkolom="par",season="maand",waarde="waarde",datum="datum",jaarkolom="jaar", grenswaarde=0.1){
#functie berekend zowel de seasonal variant als de totale variant
require(zyp)
require(trend)

outputdf <-data.frame()

mpvec <- unique(df[[mpkolom]])
#mpvec <-c("00049")
for(meetpunt in mpvec){
    subset1 <- df[df[[mpkolom]]==meetpunt,]
    parvec <-unique(subset1[[parkolom]])
    for(parameter in parvec){
        subset2 <- subset1[subset1[[parkolom]]==parameter,]
        subset2 <- subset2[order(subset2$datum),]
        aantal_jaren <- length(subset2[[jaarkolom]])
        if(nrow(subset2)<13|aantal_jaren<2){next}#ga naar de volgende loop als er te weinig data is, mogelijk nog optimaliseren

        seasonvec <- unique(subset2[[season]])

        #Mann Kendalltest zonder seizoenen
        mkrestot <- mk.test(ts(subset2[[waarde]]))
        Stot <- mkrestot[["Sg"]]
        if(Stot==0){Ztot<-0} else{Ztot <- (Stot-sign(Stot))/sqrt(mkrestot[["varSg"]])}
        probtrendtot = min(pnorm(Ztot),1-pnorm(Ztot))
        if(Stot<0&probtrendtot<grenswaarde){mktekst<-"Dalende trend"}else if(Stot>0&probtrendtot<grenswaarde){mktekst<-"Stijgende trend"}else{mktekst<-"Geen trend"}

        #Theil Sen Slope estimator
        ndatums <- as.numeric(subset2[[datum]])-min(as.numeric(subset2[[datum]]))#dag 0 is de eerste dag van de meetreeks het is de vraag of dat handig blijkt
        waarden <- subset2[[waarde]]
        aantaldagen <- max(ndatums, na.rm=TRUE)
        TStot <-zyp.sen(waarden~ndatums) #per dag, intercept is dag 1 van meetreeks
        dayslope <- TStot[[1]][2]
        dayslopeup95 <- quantile(TStot[[2]],0.95, na.rm=TRUE)
        dayslopedown05 <- quantile(TStot[[2]],0.05, na.rm=TRUE)
        intercept <- TStot[[1]][1]
        totalchange <- dayslope*aantaldagen

        sTStot<- data.frame()

        #Seasonal Mann Kendalltest en Seasonal TheilSen Slope estimator
        S = 0
        varS = 0
        for(x in seasonvec){
            subsetseason <- subset2[subset2[[season]]==x,]
            if(nrow(subsetseason)<2){next}#ga naar de volgende loop als er te weinig data is

            #Seasonal Mann Kendall part1
            mkres <- mk.test(ts(subsetseason[[waarde]]))
            S<- S+mkres[["Sg"]]
            varS <- varS+mkres[["varSg"]]

            #Seasonal Theil Sen part1
            sndatums <- as.numeric(subsetseason[[datum]])-min(as.numeric(subsetseason[[datum]]))#dag 0 is de eerste dag van de meetreeks het is de vraag of dat handig blijkt
            swaarden <- subsetseason[[waarde]]
            TStemp <-zyp.sen(swaarden~sndatums)
            tempresTS <- data.frame(seizoen=x,slope = TStemp[[1]][2],slopeup = quantile(TStemp[[2]],0.95, na.rm=TRUE),slopedown = quantile(TStemp[[2]],0.05, na.rm=TRUE),intercept=TStemp[[1]][1])
            sTStot<-rbind(sTStot,tempresTS)

        }#end of seasonal loop

        #Seasonal Mann Kendall part2
        if(S==0){Z<-0} else{Z <- (S-sign(S))/sqrt(varS)}
        probtrend = min(pnorm(Z),1-pnorm(Z))
        if(S<0&probtrend<grenswaarde){smktekst<-"Dalende trend"}else if(S>0&probtrend<grenswaarde){smktekst<-"Stijgende trend"}else{smktekst<-"Geen trend"}

        #Seasonal Theil Sen part2
        sdayslope <- median(sTStot$slope)
        #print(sdayslope)
        sdayslopeup95 <- median(sTStot$slopeup)
        sdayslopedown05 <- median(sTStot$slopedown)
        sintercept <- median(sTStot$intercept)
        stotalchange <- sdayslope*aantaldagen

        #Verwerking resultaten
        beginjaar<-min(subset2[[jaarkolom]])
        eindjaar <-max(subset2[[jaarkolom]])
        uitkomst<-data.frame(meetpunt=meetpunt, parameter=parameter,beginjaar=beginjaar,eindjaar=eindjaar,aantal_jaren=aantal_jaren,aantaldagen=aantaldagen,aantal_waarnemingen=nrow(subset2),
                            mk_Z=Ztot, mk_kans_trend=probtrendtot, mk_trendrichting=mktekst ,smk_Z=Z,smk_kans_trend=probtrend,smk_trendrichting=smktekst,
                            daghelling=dayslope, daghellingmax95=dayslopeup95, daghellingmin05=dayslopedown05, intercept=intercept, totale_verandering=totalchange,
                            s_daghelling=sdayslope, s_daghellingmax95=sdayslopeup95, s_daghellingmin05=sdayslopedown05, s_intercept=sintercept, s_totale_verandering=stotalchange)
        outputdf <- rbind(outputdf,uitkomst)

    } #end par loop
} #end mp loop

return(outputdf)

}#end function


# creert een df met genormaliseerde waarde en het verschil van de genormaliseerde waarde met de vorige en volgende
# kan sterk versimpeld worden door gebruik te maken van group_by i.c.m. mutate
controleuitschieters <- function(df, mpkolom="mp",parkolom="parnr",waarde="waarde",datum="datum",parameterdf,meetpuntendf ){
  require(dplyr)
  groupdf <- group_by(df, mp, parnr)
  stats <- summarise(groupdf,gemiddelde = mean(waarde),sd = sd(waarde))
  dfstats <- merge(df,stats, by=c("mp","parnr"),sort=FALSE)
  dfstats$genorm_waarde <-  (dfstats[["waarde"]] - dfstats[["gemiddelde"]]) / dfstats[["sd"]]
  dfstats <- arrange(dfstats,mp,parnr,datum)

  nrijen <- nrow(dfstats)
  dfstats$vorigenorm <- lag(dfstats$genorm_waarde)
  dfstats$vorigemp <- lag(dfstats$mp)
  dfstats$vorigeparnr <- lag(dfstats$parnr)
  dfstats$diffvorige <- ifelse(dfstats$mp==dfstats$vorigemp&dfstats$parnr==dfstats$vorigeparnr,abs(dfstats$genorm_waarde-dfstats$vorigenorm),0)

  dfstats$volgendenorm <- lead(dfstats$genorm_waarde)
  dfstats$volgendemp <- lead(dfstats$mp)
  dfstats$volgendeparnr <- lead(dfstats$parnr)
  dfstats$diffvolgende <- ifelse(dfstats$mp==dfstats$volgendemp&dfstats$parnr==dfstats$volgendeparnr,abs(dfstats$genorm_waarde-dfstats$volgendenorm),0)
  dfstats <- select(dfstats,-starts_with("volgende"),-starts_with("vorige"))
  dfstats
}#end of function

controle_sprongen <- function(data){
  require(dplyr)
  data %>%
    arrange(mp,parnr,datum) %>%
    group_by(mp,parnr) %>%
    mutate(verschil = waarde-lag(waarde), abs_norm_verschil = (abs(verschil)-mean(abs(verschil), na.rm=TRUE))/sd(verschil, na.rm=TRUE)) %>%
    ungroup()

}

# GRAFIEKEN PLOTTEN ----------------------------------------------------------------------

HHSKthema <- function(){
  require(ggplot2)
  hhskgroen <<- "#8dc63f"
  hhskblauw <<- "#0079c2"
  hhskthema <<- theme_light() + theme(plot.title = element_text(color = hhskgroen,face="bold",hjust=0.5),
                                     axis.title = element_text(color=hhskblauw,face="bold"),
                                     axis.text = element_text(color=hhskblauw),
                                     axis.ticks = element_line(color=hhskblauw),
                                     axis.line.x = element_line(color=hhskblauw, size=0.5),
                                     panel.border = element_rect(color = hhskblauw, size=1),
                                     panel.grid.major = element_line(color=hhskgroen,linetype="dotted", size = 0.5),
                                     panel.grid.minor = element_line(color=hhskgroen,linetype="dotted", size = 0.5)

  )
}


#één grafiek plotten
tijdreeksgrafiek <- function(df, meetpunt, parameternr, parameterdf,meetpuntendf){
  require(ggplot2)
  require(dplyr)
  require(lubridate)

  HHSKthema()

    mpinfo <- meetpuntinformatie(meetpunt,meetpuntendf=meetpuntendf)
    paraminfo <- parameterinformatie(parnr = parameternr, parameterdf = parameterdf)
      min <- min(df$waarde)
      max <- max(df$waarde)
      if(min/(max-min)>1){ylimieten <- c(min*0.9,max*1.1)}else {ylimieten <- c(0,max*1.1)}

      #GGPLOT opjecten opbouwen
      lijn <- geom_line(col=hhskblauw)
      punten <- geom_point(col=hhskblauw)
      grafiektitel <- ggtitle(paste("Meetpunt:",meetpunt,"-",mpinfo$mpomsch,"\n","Parameter:",paraminfo$parameterlang))
      y_label <- ylab(paraminfo$eenheid)
      x_label <- xlab("")
      x_axis <- scale_x_date(date_breaks = "years",labels=year)
      y_axis <- scale_y_continuous(limits=ylimieten,expand=c(0,0),oob=scales::rescale_none)#oob makes sure that the CI for LOESS is always plotted
      loess_lijn <- geom_smooth(se=TRUE,col=hhskgroen,linetype="dashed",fill=hhskblauw,alpha=0.08,fullrange=TRUE)

      plot <- ggplot(df, aes(x=datum,y=waarde))+lijn+punten+grafiektitel+x_label+y_label+x_axis+y_axis+loess_lijn+hhskthema
      plot

}#end of function

grafieken <- function(df,parameterdf,meetpuntendf, export_pad = "export/grafieken"){
  require(ggplot2)
  require(dplyr)
  require(lubridate)

  HHSKthema()

  for(meetpunt in sort(unique(df$mp))){
    print(meetpunt)
    subset1 <- dplyr::filter(df,mp==meetpunt)
    mpinfo <- meetpuntinformatie(meetpunt,meetpuntendf=meetpuntendf)

    filename <- paste0(export_pad,"/",meetpunt,".pdf")
    pdf(file=filename,width=16,height=8)
    aantalplots <- 0 #om lege plots later te verwijderen

    for (parameternr in sort(unique(subset1$parnr))){
      if(!(parameternr %in% c(1:99,107,200:401,403:505,507:899,1000:2999))){next}
      subset2 <- dplyr::filter(subset1,parnr==parameternr)
      if(nrow(subset2)<13){next}#geen grafiek bij minder dan 13 waarden
      if(min(subset2$waarde)==max(subset2$waarde)){next}#geen grafiek als alle waarden gelijk zijn
      paraminfo <- parameterinformatie(parnr = parameternr, parameterdf = parameterdf)
      min <- min(subset2$waarde)
      max <- max(subset2$waarde)
      if(min/(max-min)>1){ylimieten <- c(min*0.9,max*1.1)}else {ylimieten <- c(0,max*1.1)}

      #GGPLOT opjecten opbouwen
      lijn <- geom_line(col=hhskblauw)
      punten <- geom_point(col=hhskblauw)
      grafiektitel <- ggtitle(paste("Meetpunt:",meetpunt,"-",mpinfo$mpomsch,"\n","Parameter:",paraminfo$parameterlang))
      y_label <- ylab(paraminfo$eenheid)
      x_label <- xlab("")
      x_axis <- scale_x_date(date_breaks = "years",labels=year)
      y_axis <- scale_y_continuous(limits=ylimieten,expand=c(0,0),oob=scales::rescale_none)#oob makes sure that the CI for LOESS is always plotted
      loess_lijn <- geom_smooth(se=TRUE,col=hhskgroen,linetype="dashed",fill=hhskblauw,alpha=0.08,fullrange=TRUE)


      plot <- ggplot(subset2, aes(x=datum,y=waarde))+lijn+punten+grafiektitel+x_label+y_label+x_axis+y_axis+loess_lijn+hhskthema
      print(plot)
      aantalplots <- aantalplots + 1
    }#end parameterloop
    dev.off()
    if(aantalplots==0){file.remove(filename)}#verwijdert lege plots
  }#end meetpuntloop


}#end of function



grafiekeninternet <- function(df,mpkolom="mp",parkolom="parnr",parnaam="par",season="maand",waarde="waarde",datum="datum",jaarkolom="jaar",parameterdf,meetpuntendf){
  mpvec <- unique(df[[mpkolom]])
  #mpvec <- c("00016")

  for(meetpunt in mpvec){
    subset1 <- df[df[[mpkolom]]==meetpunt,]
    parnrvec <-sort(unique(subset1[[parkolom]]))

    #creeren bestand voor plots
    filename <- paste0("export/internet/grafieken 3-10-2017/",meetpunt,".pdf")
    pdf(file=filename,width=16,height=8)

    #ophalen meetpunt info
    print(meetpunt)
    meetpuntinfo <- meetpuntinformatie(meetpunt,meetpuntendf)
    meetpuntomsch <- meetpuntinfo[["mpomsch"]]

    aantalplots <- 0 #om lege plots later te verwijderen

    for(parameternr in parnrvec){
      if(parameternr>2999|(parameternr>99&parameternr<107)|(parameternr>107&parameternr<200)|parameternr==402|parameternr==506|(parameternr>899&parameternr<1000)){next}
      subset2 <- subset1[subset1[[parkolom]]==parameternr,]
      aantaljaar <- max(subset2[[jaarkolom]]) - min(subset2[[jaarkolom]])

      #bepalen of er een grafiek wordt gemaakt
      minwaarde <- min(subset2[[waarde]])
      maxwaarde <- max(subset2[[waarde]])

      if(nrow(subset2)<13){next}
      if(minwaarde==maxwaarde){next}

      #ophalen parameter info
      parameterinfo <- parameterinformatie(parameternr,parameterdf)
      parkort <- parameterinfo[["parameterkort"]]
      parameterlang <- parameterinfo[["parameterlang"]]
      eenheid <- parameterinfo[["eenheid"]]

      #plotinfo
      oldpar<-par()

      if(minwaarde==maxwaarde){ymin<-0} else if(minwaarde/(maxwaarde-minwaarde)>1){ymin <- minwaarde*0.9}else{ymin<-0} #ondergrens y is 90% van min als de onderste helft leeg zou zijn, anders is de ondergrens 0
      ymax <- maxwaarde*1.1

      titel <- paste(meetpunt,meetpuntomsch,"\n",parameterlang)
      par(mar=c(8, 4, 4, 2))

      #plot grafiek
      plot(subset2$datum,subset2$waarde,type="o",pch=16,lwd=1.5,col="blue",bg="blue", cex=1, ylim=c(ymin,ymax),yaxs="i",xaxs="r",
           xlab="",ylab=eenheid, main=titel,sub="", lab=c(5,5,12),
           panel.first = abline(v=pretty(subset2$datum,n=aantaljaar),lty=2,col="lightgray", h=pretty(c(ymin,ymax),n=6)))

      aantalplots <- aantalplots + 1

    } #end par loop

    # afronden bestanden
    dev.off()
    if(aantalplots==0){file.remove(filename)}#verwijdert lege plots

  } #end mp loop

}#end of function

# per meetpunt een pdf met grafieken en trendtests per meetpunt
plotting <- function(df,mpkolom="mp",parkolom="parnr",parnaam="par",season="maand",waarde="waarde",datum="datum",jaarkolom="jaar",parameterdf,meetpuntendf){
  require(ggplot2)
  require(dplyr)
  mpvec <- unique(df[[mpkolom]])
  #mpvec <- c("00016")

  for(meetpunt in mpvec){
    subset1 <- df[df[[mpkolom]]==meetpunt,]
    parnrvec <-sort(unique(subset1[[parkolom]]))

    #creeren bestand voor plots
    filename <- paste0("export/Grafieken/",meetpunt,".pdf")
    pdf(file=filename,width=16,height=8)

    #ophalen meetpunt info
    print(meetpunt)
    meetpuntinfo <- meetpuntinformatie(meetpunt,meetpuntendf)
    meetpuntomsch <- meetpuntinfo[["mpomsch"]]

    aantalplots <- 0 #om lege plots later te verwijderen

    for(parameternr in parnrvec){
      if(parameternr>2999|(parameternr>99&parameternr<107)|(parameternr>107&parameternr<200)|parameternr==402|parameternr==506|(parameternr>899&parameternr<1000)){next}
      subset2 <- subset1[subset1[[parkolom]]==parameternr,]
      parameternr <- as.integer(parameternr)
      aantaljaar <- max(subset2[[jaarkolom]]) - min(subset2[[jaarkolom]])

      #bepalen of er een grafiek wordt gemaakt
      minwaarde <- min(subset2[[waarde]])
      maxwaarde <- max(subset2[[waarde]])

      if(nrow(subset2)<13){next}
      if(minwaarde==maxwaarde){next}

      #ophalen parameter info
      parameterinfo <- parameterinformatie(parameternr,parameterdf)
      parkort <- parameterinfo[["parameterkort"]]
      parameterlang <- parameterinfo[["parameterlang"]]
      eenheid <- parameterinfo[["eenheid"]]

      #trendinfo
      trendinfo <- as.list(trendanalysemk(subset2))
      if(length(trendinfo)==0){
      typetest <- "Geen test"
      trendkans <- ""
      trendtekst <- ""
      TSintercept <- ""
      TSverandering <- ""
      }
      else if(parameternr<100|(parameternr>399&parameternr<600)){
      typetest <- "Mann-Kendalltest en Theil-Sen hellingschatter - seizoenstesten"
      trendkans <- signif(trendinfo$smk_kans_trend, digits=3)
      trendtekst <- trendinfo$smk_trendrichting
      TSintercept <- signif(trendinfo$s_intercept, digits=3)
      TSverandering <- signif(trendinfo$s_totale_verandering, digits=3)
      #aantaljaar <- trendinfo$aantal_jaren #geeft het aantal meetjaren maar niet de periode
      }else{
      typetest <- "Mann-Kendalltest en Theil-Sen hellingschatter - totaal"
      trendkans <- signif(trendinfo$mk_kans_trend, digits=3)
      trendtekst <- trendinfo$mk_trendrichting
      TSintercept <- signif(trendinfo$intercept, digits=3)
      TSverandering <- signif(trendinfo$totale_verandering, digits=3)
      #aantaljaar <- trendinfo$aantal_jaren #geeft het aantal meetjaren, maar niet de periode
      }

      #plotinfo
      oldpar<-par()

      if(minwaarde==maxwaarde){ymin<-0} else if(minwaarde/(maxwaarde-minwaarde)>1){ymin <- minwaarde*0.9}else{ymin<-0} #ondergrens y is 90% van min als de onderste helft leeg zou zijn, anders is de ondergrens 0
      ymax <- maxwaarde*1.1

      titel <- paste(meetpunt,meetpuntomsch,"\n",parameterlang)
      par(mar=c(8, 4, 4, 2))
      trendstring <- paste(typetest,"\n","\n","Mann-Kendall:",trendtekst,"met kans",trendkans,"\n","Theil-Sen: startpunt =",TSintercept,eenheid,", totale verandering =",TSverandering,eenheid,"in",aantaljaar,"jaar")
      if(typetest=="Geen test"){trendstring <- "Geen trendtest beschikbaar"}

      #plot grafiek
      plot(subset2$datum,subset2$waarde,type="l",pch=16,lwd=1.5,col="blue",bg="blue", cex=1, ylim=c(ymin,ymax),yaxs="i",xaxs="r",
           xlab="",ylab=eenheid, main=titel,sub="", lab=c(5,5,12),
           panel.first = abline(v=pretty(subset2$datum,n=aantaljaar),lty=2,col="lightgray", h=pretty(c(ymin,ymax),n=6)))

      mtext(trendstring,side=1,line=6)
      #LOESS-model
      loessspan <- min(1/(aantaljaar**0.3),1, na.rm = TRUE)
      loessmodel<-loess(subset2$waarde~as.numeric(subset2$datum),span=loessspan)
      lines(subset2[[datum]],fitted(loessmodel),col="red")
      aantalplots <- aantalplots + 1

    } #end par loop

    # afronden bestanden
    dev.off()
    if(aantalplots==0){file.remove(filename)}#verwijdert lege plots

  } #end mp loop

}#end of function

#histogrammen per parameter, maak vooraf een df met de group die interessant is, ook mogelijk is een loop met subsets te creeren waarin deze functie steeds wordt aangeroepen
histbasis <- function(df,naam="Schieland",parkolom="parnr",waarde="waarde",jaarkolom="jaar",parameterdf,meetpuntendf){
    graphics.off()
    parnrvec <-sort(unique(df[[parkolom]]))

    #creeren bestand voor plots
    filename <- paste0("export/Histogrammen/",naam,".pdf")
    pdf(file=filename,width=10,height=6)

    #ophalen meetpunt info
    for(parameternr in parnrvec){
        parsubset <- df[df[[parkolom]]==parameternr,]
        aantaljaar <- max(parsubset[[jaarkolom]]) - min(parsubset[[jaarkolom]]) #nodig?

        #bepalen of er een grafiek wordt gemaakt
        minwaarde <- min(parsubset[[waarde]])
        maxwaarde <- max(parsubset[[waarde]])
        p99waarde <- quantile(parsubset[[waarde]],0.99)
        if(nrow(parsubset)<13){next}
        if(minwaarde==maxwaarde){next}
        if(parameternr>999|(parameternr>99&parameternr<200)|parameternr==402|(parameternr>899&parameternr<1000)){next}

        print(parameternr)
        #ophalen parameter info
        parameterinfo <- parameterinformatie(parameternr,parameterdf)
        parkort <- parameterinfo[["parameterkort"]]
        parameterlang <- parameterinfo[["parameterlang"]]
        eenheid <- parameterinfo[["eenheid"]]

        nrbreaks = (maxwaarde-minwaarde)/(p99waarde-minwaarde)*20

        titel <- paste(naam,"-",parameterlang)
        breuken <- pretty(c(minwaarde,maxwaarde),n=nrbreaks)
        colrange <- colorRampPalette(c("lightskyblue1","navy")) #creeert een FUNCTIE
        if(minwaarde==maxwaarde){xmin<-0} else if(minwaarde/(maxwaarde-minwaarde)>1){xmin <- minwaarde*0.9}else{xmin<-0} #ondergrens y is 90% van min als de onderste helft leeg zou zijn, anders is de ondergrens 0

        hist(parsubset[["waarde"]],main=titel ,xlim = c(xmin,p99waarde),breaks=breuken, xlab=eenheid, ylab="Aantal metingen", col=colrange(28))
        #abline(v=quantile(parsubset[["waarde"]],0.5),col="navy",lwd=1.5,lty=3)
        #mtext("50%",col="navy",side=1,at=quantile(parsubset[["waarde"]],0.5),line=2.5)

    } #end par loop
    # afronden bestanden
    dev.off()

    graphics.off()
}#end of function

# histogrammen per parameter waarin alle groepen weergegeven worden
histgroups <- function(df,mpkolom="mp",parkolom="parnr",parnaam="par",season="maand",waarde="waarde",datum="datum",jaarkolom="jaar",parameterdf,meetpuntendf,groupkolom="gebied"){

    df<- df[complete.cases(df[[groupkolom]]),]

    groupvec <- unique(df[[groupkolom]])
    nrgroups <- length(groupvec)

    parnrvec <-sort(unique(df[[parkolom]]))
    #parnrvec <- c(1)
    #creeren bestand voor plots
    filename <- paste0("export/Grafieken/",groupkolom,".pdf")
    pdf(file=filename,width=7,height=4*nrgroups)

    #ophalen meetpunt info
    for(parameternr in parnrvec){
      subset1 <- df[df[[parkolom]]==parameternr,]
      parameternr <- as.integer(parameternr)
      aantaljaar <- max(subset1[[jaarkolom]]) - min(subset1[[jaarkolom]])

      #bepalen of er een grafiek wordt gemaakt
      minwaarde <- min(subset1[[waarde]])
      maxwaarde <- max(subset1[[waarde]])
      p99waarde <- quantile(subset1[[waarde]],0.99)
      if(nrow(subset1)<13){next}
      if(minwaarde==maxwaarde){next}
      if(parameternr>999|(parameternr>99&parameternr<200)|parameternr==402|(parameternr>899&parameternr<1000)){next}
      print(parameternr)
      #ophalen parameter info
      parameterinfo <- parameterinformatie(parameternr,parameterdf)
      parkort <- parameterinfo[["parameterkort"]]
      parameterlang <- parameterinfo[["parameterlang"]]
      eenheid <- parameterinfo[["eenheid"]]
      oldpar <- par()
      par(mfrow=c(nrgroups,1))
      nrbreaks = (maxwaarde-minwaarde)/(p99waarde-minwaarde)*20


      for(x in groupvec){
        ifelse(x=="SCHIEL",titelnaam<-"Schieland",ifelse(x=="KRIWA",titelnaam<-"Krimpenerwaard",titelnaam<-x))
        subset2 <- subset1[subset1[[groupkolom]]==x,]
        titel <- paste(titelnaam,"\n",parameterlang)
        breuken <- pretty(c(minwaarde,maxwaarde),n=nrbreaks)
        colrange <- colorRampPalette(c("lightskyblue1","navy"))
        if(minwaarde==maxwaarde){xmin<-0} else if(minwaarde/(maxwaarde-minwaarde)>1){xmin <- minwaarde*0.9}else{xmin<-0} #ondergrens y is 90% van min als de onderste helft leeg zou zijn, anders is de ondergrens 0

        hist(subset2[["waarde"]],main=titel ,xlim = c(xmin,p99waarde),breaks=breuken, xlab=eenheid, ylab="Aantal metingen", col=colrange(28))

      }#end group loop

    } #end par loop

    # afronden bestanden
    dev.off()

}#end of function

#histogrammen waarbij 1 group tegen de rest wordt uitgezet.
hist2groups <- function(df,mpkolom="mp",parkolom="parnr",parnaam="par",season="maand",waarde="waarde",datum="datum",jaarkolom="jaar",parameterdf,meetpuntendf,groupkolom="gebied",group="Schieland"){
    # 1 tegen de rest
    df<- df[complete.cases(df[[groupkolom]]),]

    parnrvec <-sort(unique(df[[parkolom]]))
    #parnrvec <- c(7)
    #creeren bestand voor plots
    filename <- paste0("export/Grafieken/",group,".pdf")
    pdf(file=filename,width=7,height=8)
    oldpar <- par()
    par(mfrow=c(2,1))
    #ophalen meetpunt info

    for(parameternr in parnrvec){
      subset1 <- df[df[[parkolom]]==parameternr,]
      parameternr <- as.integer(parameternr)
      aantaljaar <- max(subset1[[jaarkolom]]) - min(subset1[[jaarkolom]])

      #bepalen of er een grafiek wordt gemaakt
      minwaarde <- min(subset1[[waarde]])
      maxwaarde <- max(subset1[[waarde]])
      p99waarde <- quantile(subset1[[waarde]],0.99)
      if(nrow(subset1)<13){next}
      if(minwaarde==maxwaarde){next}
      if(parameternr>999|(parameternr>99&parameternr<200)|parameternr==402|(parameternr>899&parameternr<1000)){next}
      #print(parameternr)
      #ophalen parameter info
      parameterinfo <- parameterinformatie(parameternr,parameterdf)
      parkort <- parameterinfo[["parameterkort"]]
      parameterlang <- parameterinfo[["parameterlang"]]
      eenheid <- parameterinfo[["eenheid"]]

      nrbreaks = (maxwaarde-minwaarde)/(p99waarde-minwaarde)*20

      subset2 <- subset1[subset1[[groupkolom]]==group,]
      inversesubset2 <- subset1[subset1[[groupkolom]]!=group,]
      titel <- paste(parameterlang,"\n\n",group)
      titelinverse <- paste0("Niet-",group)
      #titelinverse <- paste0("De rest van Schieland")
      breuken <- pretty(c(minwaarde,maxwaarde),n=nrbreaks)
      colrange <- colorRampPalette(c("lightskyblue1","navy"))
      if(minwaarde==maxwaarde){xmin<-0} else if(minwaarde/(maxwaarde-minwaarde)>1){xmin <- minwaarde*0.9}else{xmin<-0} #ondergrens y is 90% van min als de onderste helft leeg zou zijn, anders is de ondergrens 0
      hist(subset2[["waarde"]],main=titel ,xlim = c(xmin,p99waarde),breaks=breuken, xlab=eenheid, ylab="Aantal metingen",col=colrange(28))
      hist(inversesubset2[["waarde"]],main=titelinverse ,xlim = c(xmin,p99waarde),breaks=breuken, xlab=eenheid,ylab="Aantal metingen",col=colrange(28))


    } #end par loop

    # afronden bestanden
    dev.off()

}#end of function


# GIS ----------------------------------------------------------------------

#geef een object een projectie
project <- function(df){
  require(sp)
  myCRS <- "+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +towgs84=565.2369,50.0087,465.658, -0.406857330322398,0.350732676542563,-1.8703473836068, 4.0812"
  proj4string(df) <- CRS(myCRS)

  df
}

#maak van een df meet coordinaten een spdf - wist ongeldige waarden
set.coords = function(df,X="X",Y="Y"){
  # Function converts dataframe to spatialpointsdataframe
  require(sp)

  #Remove empty coordinates
  df <- df[complete.cases(df$X),]
  df <- df[complete.cases(df$Y),]
  df <- df[df$X!=0&df$Y!=0,]

  # Convert to spatial points data frame
  coordinates(df) = ~X+Y

  # Assign the datum and compute the lat-lon of the x-y coordinates
  myCRS <- "+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +towgs84=565.2369,50.0087,465.658, -0.406857330322398,0.350732676542563,-1.8703473836068, 4.0812"
  proj4string(df) <- CRS(myCRS)
  df
}

#voegt coordinaten aan een df met meetpunten en zet om in een spdf
add.coords = function(df,meetpuntendf){
  # Function adds coordinates and converts dataframe to spatialpointsdataframe
  require(sp)
  require(dplyr)
  meetpuntendfkort <- select(meetpuntendf, mp, X,Y)
  df <- merge(df, meetpuntendfkort, by="mp")

  #Remove empty coordinates
  df <- df[complete.cases(df$X),]
  df <- df[complete.cases(df$Y),]
  df <- df[df$X!=0&df$Y!=0,]


  # Convert to spatial points data frame
  coordinates(df) = ~X+Y
  myCRS <- "+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +towgs84=565.2369,50.0087,465.658, -0.406857330322398,0.350732676542563,-1.8703473836068, 4.0812"
  # Assign the datum and compute the lat-lon of the x-y coordinates
  proj4string(df) <- CRS(myCRS)
  df
}

# wijzigt de CRS - borgt dat overal dezelfde CRS wordt gebruikt
reset.CRS <- function(spatialobject){ #functie borgt gelijke coordinaatsystemen
  require(rgdal)
  myCRS <- "+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +towgs84=565.2369,50.0087,465.658, -0.406857330322398,0.350732676542563,-1.8703473836068, 4.0812"
  spatialobject <- spTransform(spatialobject,myCRS)
  return(spatialobject)
}#end of function

convert_to_wgs84 <- function(spatialobject){
  require(rgdal)
  spatialobject <- spTransform(spatialobject,"+init=EPSG:4326")
  return(spatialobject)
}

#voegt aan een sp meetpuntendf de info van een (polygonen)laag toe
spatialjoinmp <- function(mp,layer){
  require(sp)
  require(dplyr)

  mp <- reset.CRS(mp)#zou overbodig moeten zijn
  layer <- reset.CRS(layer)
  joined <- over(mp,layer)
  mpjoined <- cbind(data.frame(mp),data.frame(joined)) %>%select(-optional)
  set.coords(mpjoined)

}#end of function

# UTILITIES ----------------------------------------------------------------------

df_to_named_list <- function(df, waarden=1, namen=2){
  require(dplyr)
  values <- c(dplyr::select(df, waarden), use.names = FALSE, recursive = TRUE)
  names_values <- c(dplyr::select(df, namen), use.names = FALSE, recursive = TRUE)
  names(values) <- names_values
  values
}




# EINDE FUNCTIES ----------------------------------------------------------------------


graphics.off()
#plotting(df,parameterdf=parameterdf,meetpuntendf=meetpuntendf)
#histbasis(df,naam="Kralingse Plas",parameterdf=parameterdf,meetpuntendf=meetpuntendf)
#write.table(x, file="export/trendanalyse.csv",quote=FALSE,sep=";",dec=",",row.names=FALSE)


