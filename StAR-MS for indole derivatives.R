#####
##### 0. Libraries and descriptions
### 0.1 Libraries
library(tidyverse)

### 0.2 Descriptions
rm(list = ls()) 
#####
##### 1. Import data
a <- 'C:/ZHD/UNC/UNC_Research/AHR_Project/IndolePRDataAnalysis/IndolePRZHD/89_4.1.csv' %>% 
  read_csv()

### 1.1 Write a function to tidy imported peak lists (in total 64 peak lists)
plrt.tidy <- function(csv_path = 'C:/ZHD/UNC/UNC_Research/AHR_Project/IndolePRDataAnalysis/IndolePRZHD/89_4.23.csv')
{
  #csv_path <- '/Users/aaronthecloud/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC_Research/LuLabResearch/OtherLabMemberStudy/HDZ/IndolePrecursorScan/IndolePrecursorScanList/89_4.25-4.5.csv'
  d <- csv_path %>% 
    read_csv(file = csv_path)
  rfn <- d[2 , 1] %>% unlist() %>% as.character()
  pr <- d[3 , 1] %>% unlist() %>% as.character()
  sr <- d[4 , 1] %>% unlist() %>% as.character()
  rtr <- d[5 , 1] %>% unlist() %>% as.character()
  
  d1 <- d[9:nrow(d) , ]
  colnames(d1) <- c('precursor m/z','intensity')
  d1 <- d1 %>% 
    mutate(. , 
           `rawfile name` = rfn , 
           `precursor scan` = pr , 
           `scan range` = sr , 
           `retention time range` = rtr , 
           `retention time` = substr(rtr, 5, 8),
           `retention time` = as.numeric(`retention time`),
           `retention time range left boundary` = as.integer(floor(`retention time`/0.1))/10,
           `retention time range right boundary` = `retention time range left boundary`+0.1,
           `retention time range left boundary` = as.numeric(`retention time range left boundary`) , 
           `retention time range right boundary` = as.numeric(`retention time range right boundary`))
  return(d1)
}

### 1.2 rbind all imported peak lists
csv.vector <- list.files(path = 'C:/ZHD/UNC/UNC_Research/AHR_Project/IndolePRDataAnalysis/IndolePRZHD/')
plrt.tidy.total <- tibble()
for(i in 1:length(csv.vector))
{
  print(i)
  p <- csv.vector[i]
  plrt.tidy.i <- p %>% 
    paste0('C:/ZHD/UNC/UNC_Research/AHR_Project/IndolePRDataAnalysis/IndolePRZHD/', .) %>% 
    plrt.tidy(csv_path = .)
  
  if(i == 1)
  {
    plrt.tidy.total <- plrt.tidy.i
  }
  else
  {
    plrt.tidy.total <- rbind(plrt.tidy.total , 
                             plrt.tidy.i)
  }
}
write.csv(plrt.tidy.total, "plrt.tidy.total.csv")
##### clean the bug data, then import
##### 2. Filter combined peak list
### 2.1 Filter by the unit mass and keep the row with the highest intensity
plrt.tidy.total<- read.csv("plrt.tidy.total.csv")
plrt.tidy.total.1 <- plrt.tidy.total %>% # nrow=6400
  mutate(. , 
         rownum = 1:nrow(.) , 
         intensity = as.numeric(intensity) , 
         `precursor scan target` = ifelse(str_detect(pattern = 'pr 89' , 
                                                     string = `precursor.scan`) , 
                                          89 , 
                                          ifelse(str_detect(pattern = 'pr 103' , 
                                                            string = `precursor.scan`) , 
                                                 103 , 
                                                 ifelse(str_detect(pattern = 'pr 116' , 
                                                                   string = `precursor.scan`) , 
                                                        116 , 
                                                        ifelse(str_detect(pattern = 'pr 130' , 
                                                                          string = `precursor.scan`) , 
                                                               130 , NA)))) , 
         `precursor.m.z` = round(as.numeric(`precursor.m.z`) , 0)) %>% 
  arrange(desc(intensity) , rownum) %>% 
  group_by(`precursor.m.z` , `rawfile.name` , `precursor.scan` , 
           `scan.range` , `retention.time.range`) %>% 
  slice(1) %>% 
  arrange(rownum) # nrow=2417

### 2.2 Build summary table of the PRM target candidates
prm.il.summary.0 <- plrt.tidy.total.1 %>% 
  as_tibble() %>% 
  select(. , 
         `precursor.m.z` , `retention.time.range` , `retention.time.range.left.boundary` , `retention.time.range.right.boundary`) %>% 
  unique()

pr89.detected.total <- c()
pr103.detected.total <- c()
pr116.detected.total <- c()
pr130.detected.total <- c()
pr89.intensity.total <- c()
pr103.intensity.total <- c()
pr116.intensity.total <- c()
pr130.intensity.total <- c()
for(i in 1:nrow(prm.il.summary.0))
{
  print(i)
  #i <- 1
  pr <- prm.il.summary.0[i , ]$`precursor.m.z`
  rtr <- prm.il.summary.0[i , ]$`retention.time.range.left.boundary`
  
  prlist.89.i <- plrt.tidy.total.1 %>% 
    filter(. , 
           `precursor.m.z` == pr & 
             `retention.time.range.left.boundary` == rtr & 
             `precursor scan target` == 89)
  prlist.103.i <- plrt.tidy.total.1 %>% 
    filter(. , 
           `precursor.m.z` == pr & 
             `retention.time.range.left.boundary` == rtr & 
             `precursor scan target` == 103)
  prlist.116.i <- plrt.tidy.total.1 %>% 
    filter(. , 
           `precursor.m.z` == pr & 
             `retention.time.range.left.boundary` == rtr & 
             `precursor scan target` == 116)
  prlist.130.i <- plrt.tidy.total.1 %>% 
    filter(. , 
           `precursor.m.z` == pr & 
             `retention.time.range.left.boundary` == rtr & 
             `precursor scan target` == 130)
  
  pr89.detected.i <- ifelse(nrow(prlist.89.i) == 0 , FALSE , TRUE)
  pr103.detected.i <- ifelse(nrow(prlist.103.i) == 0 , FALSE , TRUE)
  pr116.detected.i <- ifelse(nrow(prlist.116.i) == 0 , FALSE , TRUE)
  pr130.detected.i <- ifelse(nrow(prlist.130.i) == 0 , FALSE , TRUE)
  
  pr89.intensity.i <- ifelse(nrow(prlist.89.i) == 0 , NA , 
                             sort(prlist.89.i$intensity , decreasing = TRUE)[1])
  pr103.intensity.i <- ifelse(nrow(prlist.103.i) == 0 , NA , 
                         sort(prlist.103.i$intensity , decreasing = TRUE)[1])
  pr116.intensity.i <- ifelse(nrow(prlist.116.i) == 0 , NA , 
                              sort(prlist.116.i$intensity , decreasing = TRUE)[1])
  pr130.intensity.i <- ifelse(nrow(prlist.130.i) == 0 , NA , 
                              sort(prlist.130.i$intensity , decreasing = TRUE)[1])
  
  pr89.detected.total <- c(pr89.detected.total , pr89.detected.i)
  pr103.detected.total <- c(pr103.detected.total , pr103.detected.i)
  pr116.detected.total <- c(pr116.detected.total , pr116.detected.i)
  pr130.detected.total <- c(pr130.detected.total , pr130.detected.i)
  pr89.intensity.total <- c(pr89.intensity.total , pr89.intensity.i)
  pr103.intensity.total <- c(pr103.intensity.total , pr103.intensity.i)
  pr116.intensity.total <- c(pr116.intensity.total , pr116.intensity.i)
  pr130.intensity.total <- c(pr130.intensity.total , pr130.intensity.i)
}

prm.il.summary.1 <- prm.il.summary.0 %>% 
  mutate(. , 
         `product ion detection: 89` = pr89.detected.total , 
         `product ion detection: 103` = pr103.detected.total , 
         `product ion detection: 116` = pr116.detected.total , 
         `product ion detection: 130` = pr130.detected.total , 
         `product ion detection: total detected product ions` = pr89.detected.total + pr103.detected.total + pr116.detected.total + pr130.detected.total , 
         `product ion intensity: 89` = pr89.intensity.total , 
         `product ion intensity: 103` = pr103.intensity.total , 
         `product ion intensity: 116` = pr116.intensity.total , 
         `product ion intensity: 130` = pr130.intensity.total , 
         `product ion intensity: intensity sum of total intensity product ions` = rowSums(tibble(a = pr89.intensity.total , 
                                                                                                 b = pr103.intensity.total , 
                                                                                                 c = pr116.intensity.total , 
                                                                                                 d = pr130.intensity.total) , 
                                                                                          na.rm = TRUE)) %>% 
  arrange(.data = . , 
          desc(`product ion detection: total detected product ions`) , 
          desc(`product ion intensity: intensity sum of total intensity product ions`))

table(prm.il.summary.1$`product ion detection: total detected product ions`)

write_csv(x = prm.il.summary.1 , 
          file ='C:/ZHD/UNC/UNC_Research/AHR_Project/IndolePRDataAnalysis/Newzhd.csv')
rm(list = ls()) 
a<-read.csv("C:/ZHD/UNC/UNC_Research/AHR_Project/IndolePRDataAnalysis/IndolePR3.0/IndolePR330.csv")
Fragment3<-filter(a, a$"product.ion.detection..total.detected.product.ions" ==3)
#Select MS1 with four fragments
Fragment4<-filter(a, a$"product.ion.detection..total.detected.product.ions" ==4)
Frag130_any<-filter(a, a$"product.ion.detection..130" ==TRUE)
#Select MS1 with two fragments with specific pairs
a1<-filter(a, a$"product.ion.detection..total.detected.product.ions" ==2)
Frag89_116<-filter(a1, a1$"product.ion.detection..89" ==TRUE & a1$"product.ion.detection..116" ==TRUE)
Frag103_130<-filter(a1, a1$"product.ion.detection..103" ==TRUE & a1$"product.ion.detection..130" ==TRUE)
a2<-filter(a, a$"product.ion.detection..total.detected.product.ions" ==1)
Frag130<-filter(a2, a2$"product.ion.detection..130" ==TRUE)
Frag89<-filter(a2, a2$"product.ion.detection..89" ==TRUE)
Frag116<-filter(a2, a2$"product.ion.detection..116" ==TRUE)
Frag103<-filter(a2, a2$"product.ion.detection..103" ==TRUE)
#Combine with Orbitrap MS1 data, generate more accurate inclusion list
b<-read.csv("IndoleInFecesMS1.csv")
Frag103_130New<-merge(Frag103_130,b,by=c("m.z_Floor","RT_Floor"))
Frag89_116New<-merge(Frag89_116,b,by=c("m.z_Floor","RT_Floor"))
Frag3New<-merge(Fragment3,b,by=c("m.z_Floor","RT_Floor"))
Frag4New<-merge(Fragment4,b,by=c("m.z_Floor","RT_Floor"))
Frag130_anyNew<-merge(Frag130_any,b,by=c("m.z_Floor","RT_Floor"))
write.csv(Frag103_130New,"Frag103_130NewPrecise.csv")
write.csv(Frag89_116New,"Frag89_116NewPrecise.csv")
write.csv(Frag3New,"Fragment3NewPrecise.csv")
write.csv(Frag4New,"Fragment4NewPrecise.csv")
write.csv(Frag130New,"Frag130NewPrecise.csv")
