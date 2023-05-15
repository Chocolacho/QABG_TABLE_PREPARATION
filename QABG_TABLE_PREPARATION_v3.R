library("tidyverse")

#### Compilation and selection ATB

#### Choice your specie and year
SPEC = "Klebsiella pneumoniae"
YEAR = "2019"
DIR = paste(".",SPEC,YEAR,sep="/")
if(SPEC=="Klebsiella pneumoniae"){AB="KP"}

#### Select the number of ATB to extract
N_ATB =6

setwd(DIR)
getwd()
#### Extract your raw data (qABG) and group (Enz, ST, Ag...)
List_ATB = read.table("List_ATB.txt", sep="\t")
List_ATB <- as.vector(List_ATB)
List_ATB <- unlist(List_ATB)


#### Import your data (verify one name column commun) and ajust data
RAW_DATA <- read_tsv("ALL_BLSE_2019.txt")
RAW_DATA$Germe[RAW_DATA$Germe %in% "Klebsiella pneumoniae pneumoniae"] <- "Klebsiella pneumoniae"
RAW_DATA$Germe[RAW_DATA$Germe %in% c("Enterobacter cloacae",
                                     "Klebsiella aerogenes (Enterobacter aero)")] <- "Enterobacter cloacae complex"
RAW_DATA <- RAW_DATA %>% select(-"Date de prél.",-"Date entrée",
                                -"Topographies", -"Protocole",-List_ATB) %>%
  filter(Germe==SPEC)
N= 13
for(i in List_ATB){
  NEW_NAME = i
  NEW_NAME
  OLD_NAME = paste("Valeur...",N,sep="")
  OLD_NAME
  names(RAW_DATA)[names(RAW_DATA) == OLD_NAME] <- NEW_NAME
  TYPE = paste("TYPE",NEW_NAME,sep="_")
  RAW_DATA <- separate(RAW_DATA, col=NEW_NAME, into = c(NEW_NAME,TYPE), sep = " ")
  N=N+2
}
rm(i, N, NEW_NAME, OLD_NAME, TYPE)
DATA_SELECT <- read_tsv(paste("DATA_UF0358_",AB,"_",YEAR,".txt",sep=""))
DATA_SELECT <- DATA_SELECT %>% select(-"Numero.sejour",-"Origine",-"Rank", -"Support")
DATA_SELECT <- DATA_SELECT %>% filter(Taxon==SPEC)
DATA <- merge(x = RAW_DATA, y = DATA_SELECT, by = "N.demande")

#### Clean data by remove null column and convert min and max of MIC in numeric value
DATA_CLEAN <- DATA %>% select_if(~!all(is.na(.)))
DATA_CLEAN <- DATA_CLEAN %>% mutate_if(is.character, ~ifelse(grepl("[<>=]",.),
                                                                    as.numeric(gsub("[><=].*","",.)),.))


#### Highlight ATB (ABG data)
List_ATB = names(DATA_CLEAN)[grep("^9",names(DATA_CLEAN))]
TABLE_ABG <- data.frame(matrix(ncol=0, nrow=7))
TABLE_MIC <- data.frame(matrix(ncol=0, nrow=7))
for(i in List_ATB){
  ATB = i
  TYPE = paste("TYPE",ATB,sep="_")
  PROV <- DATA_CLEAN %>% select(ATB, TYPE) %>% drop_na()
  PROV <- PROV[(PROV[,2]=="mm"),]
  PROV <- as.numeric(PROV[,1])
  if(length(PROV)>0){
    Col_PROV <- as.data.frame(unclass(summary(PROV)))
    names(Col_PROV)[names(Col_PROV) == "unclass(summary(PROV))"] <- ATB
    Col_PROV <- rbind(Col_PROV, length(PROV))
    row.names(Col_PROV)[7] <-"N"
    TABLE_ABG <- cbind(TABLE_ABG,Col_PROV)
  }
  PROV <- DATA_CLEAN %>% select(ATB, TYPE) %>% drop_na()
  PROV <- PROV[(PROV[,2]=="mg/L"),]
  PROV <- as.numeric(PROV[,1])
  PROV <- log2(PROV)
  if(length(PROV)>0){
    Col_PROV <- as.data.frame(unclass(summary(PROV)))
    names(Col_PROV)[names(Col_PROV) == "unclass(summary(PROV))"] <- ATB
    Col_PROV <- rbind(Col_PROV, length(PROV))
    row.names(Col_PROV)[7] <-"N"
    TABLE_MIC <- cbind(TABLE_MIC,Col_PROV)
  }
}
rm(PROV, i, TYPE, ATB, Col_PROV)
ECART_ABG <- TABLE_ABG[5,]-TABLE_ABG[2,]
TABLE_ABG <- rbind(TABLE_ABG, ECART_ABG)
row.names(TABLE_ABG)[8] <-"ECARTQ1_Q3"
ECART_MIC <- TABLE_MIC[5,]-TABLE_MIC[2,]
TABLE_MIC <- rbind(TABLE_MIC, ECART_MIC)
row.names(TABLE_MIC)[8] <-"ECARTQ1_Q3"

#### BEST ATB
BEST_ABG <- TABLE_ABG[8,] %>% unlist() %>% 
  sort(decreasing = TRUE) %>% head(6) %>% names() %>% as.character()
BEST_MIC <- TABLE_MIC[8,] %>% unlist() %>% 
  sort(decreasing = TRUE) %>% head(6) %>% names() %>% as.character()
#### CHOICE ATB by ABG or MIC
ATB = c("9CRO1","9CAZ1","9GM1","9CIP1","9SXT1")
TYPE = "ABG"
if(TYPE == "ABG"){
  UNIT="mm"
  ATB <- BEST_ABG
  TABLE_ABG_SELECT <- TABLE_ABG %>% select(ATB)
}
if(TYPE=="MIC"){
  UNIT="mg/L"
  ATB <- BEST_MIC
  TABLE_MIC_SELECT <- TABLE_MIC %>% select(ATB)
  }

MATRIX_ATB <- DATA_CLEAN %>% select(Nom.etude, ATB, paste("TYPE_",ATB, sep=""))
MATRIX <- MATRIX_ATB %>% rowwise() %>%
  filter(sum(c_across(-1:-N_ATB) == UNIT, na.rm = TRUE)==N_ATB) %>%
  select(0:N_ATB+1) %>% ungroup()
write_tsv(x=MATRIX, file = paste("MATRIX_", TYPE, "_", SPEC, "_", YEAR, ".txt", sep=""))
Strain_to_comp <- MATRIX %>% select(Nom.etude) %>% unlist()
ST_group <- DATA_SELECT %>% select(Nom.etude, ST) %>% filter(Nom.etude %in% Strain_to_comp)
write_tsv(x=ST_group, file = paste("ST_group_", TYPE, "_", SPEC, "_", YEAR, ".txt", sep=""))
BL_group <- DATA_SELECT %>% select(Nom.etude, Beta.lactam) %>% filter(Nom.etude %in% Strain_to_comp)
write_tsv(x=BL_group, file = paste("betalactam_group_", TYPE, "_", SPEC, "_", YEAR, ".txt", sep=""))
