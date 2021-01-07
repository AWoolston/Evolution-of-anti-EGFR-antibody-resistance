##############################################################################
##############################################################################

## CODE TO REPRODUCE FIGURE S7

##############################################################################
##############################################################################

setwd("/Users/awoolston/Desktop/Reanalysis_Versions/Latest/GitHub/")

## Synapse ID - syn12026190
cosmic_filename<-'COSMIC_Mutational_Signatures_v3.1_SBS_GRCh37_WXS.xlsx'

## Exome regions file
target_filename<-'SSV5_regions.bed'

## signature weights as shown in supplementary table S2
weights_filename<-"Supplementary_Table_S2.xlsx"

########################### Set Patients/Signatures ###########################

## Select BL sample from tumours with prolonged benefit
selected.cases<-c('C1005BL', 'C1007BL', 'C1014BL', 'C1018BL', 'C1020BL', 'C1024BL', 'C1025BL', 'C1026BL', 'C1027BL', 'C1030BL', 'C1037BL', 'C1044BL')

## Select COSMIC signatures
selected.signatures<- c('SBS1','SBS3','SBS5','SBS6','SBS15','SBS17a','SBS17b','SBS35','SBS40')

################## Calculate Signature Exposures ###############################

library('readxl')

## Load the canonical signature profiles
Canonical_Signatures <- as.matrix(read_excel(cosmic_filename)); 

## set rownames to common trinucleotide mutation context notation
rownames(Canonical_Signatures)<-paste0(substr(Canonical_Signatures[,'SubType'],1,1),'[',Canonical_Signatures[,'Type'],']',substr(Canonical_Signatures[,'SubType'],3,3));

## extract signatures of interest
Canonical_Signatures_selected<-t(Canonical_Signatures)[selected.signatures,]
class(Canonical_Signatures_selected)<-"numeric"

## Load signature data and assign signature labels as rownames (see Supplementary Table S2)
Signature.Matrix <- as.matrix(read_excel(weights_filename))
rownames(Signature.Matrix)<-Signature.Matrix[,'Signatures']

## Extract assigned signature weights and convert to numeric
Signature.Weights<-Signature.Matrix[selected.signatures,][,selected.cases]
class(Signature.Weights)<-"numeric"

## Extract the explained variance
R2<-as.numeric(Signature.Matrix['R2',][selected.cases])

## Adjust signature weights by the R2 value
Signature.Weights.mat.adjusted<-t(t(Signature.Weights) * R2)

## Extract the mutation load of each sample
Mutation.Load<-as.numeric(Signature.Matrix['MutationLoad',][selected.cases])

## Calculate the mutation load attributed to each signature
Signature.AssignedMutations<-rbind(t(t(Signature.Weights.mat.adjusted) * Mutation.Load),"unexplained proportion"=(1-R2)*Mutation.Load);

## Prolonged Benefit BL
BL.ProlongedBenefit.Summary<-rowSums(Signature.AssignedMutations[,selected.cases])

## Calculate the overall proportional signature contributions
prop.sig.contribution<-(rowSums(Signature.AssignedMutations)/sum(BL.ProlongedBenefit.Summary))[selected.signatures]

## Scale the signature constributions to sum to the median mutation load
prop.sig.contribution.nomiss<-prop.sig.contribution*(1/sum(prop.sig.contribution))

######################## Normalize Signatures to Exome ########################

library('GenomicRanges')
library('SigsPack')

## determine the trinucleotide across the exome (we are using the signatures derived exclusively from exome samples and so the conversion is exome context to even). 
regions <- read.delim(target_filename, header=FALSE, stringsAsFactors=FALSE)

## convert to Granges object
gr<-GenomicRanges::GRanges(seqnames=regions[,1],ranges=IRanges::IRanges(start=regions[,2],end=regions[,3]),strand=c("+"))

## gather all trinucleotide contexts across the whole exome
exome_context <- SigsPack::get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)

## set the conversion proportions
source_context<-target_context<-exome_context
target_context[,1]<-sum(source_context)/length(source_context)

Canonical_Signatures_selected_ExomeScaled<-SigsPack::normalize(t(Canonical_Signatures_selected), source_context, target_context = target_context)

################################ Figure S7A ################################

## weight the signature contributions by observed proportions
select.default.weighted<-as.matrix(prop.sig.contribution.nomiss*t(Canonical_Signatures_selected_ExomeScaled))
class(select.default.weighted)<-"numeric";

colours<-rbind(c("SBS1","lightblue"),
               c("SBS3","darkblue"),
               c("SBS5","lightgreen"),
               c("SBS6","yellow"),
               c("SBS15","darkgreen"),
               c("SBS17a","orangered"),
               c("SBS17b","red"),
               c("SBS35","purple"),
               c("SBS40","orange"))

colour_code<-colours[,2][match(rownames(select.default.weighted),colours[,1])]

## plot modelled mutation profile
barplot(select.default.weighted,col=colour_code)

########################### Figure S7B #####################################

######################### reported hotspot mutations ################################

## Prospect C reported mutations (Khan et al., Cancer Discovery, 2018)
Prospect1_KnownMuts<-rbind(c("chr12","25398285","25398285","C","A","C1005","KRAS","G12C"), ## KRAS G12C
                        c("chr1","115256530","115256530","G","T","C1005","NRAS","Q61K"), ## NRAS Q61K
                        c("chr12","25398285","25398285","C","A","C1044","KRAS","G12C"),  ## KRAS G12C
                        c("chr12","25380275","25380275","T","G","C1044","KRAS","Q61H"),  ## KRAS Q61H
                        c("chr1","115256530","115256530","G","T","C1044","NRAS","Q61K") ## NRAS Q61K
)

## Prospect C mutations (Woolston et al, Cancer Cell, 2019)
Prospect2_KnownMuts<-rbind(c("chr1","115256530","115256530","G","T","C1041","NRAS","Q61K"), ## NRAS Q61K
                          c("chr12","25398284","25398284","C","A","C1030","KRAS","G12V"),  ## KRAS G12V
                          c("chr12","25398284","25398284","C","T","C1030","KRAS","G12D"),  ## KRAS G12D
                          c("chr12","25398281","25398281","C","T","C1030","KRAS","G13D"),  ## KRAS G13D
                          c("chr12","25380275","25380275","T","G","C1025","KRAS","Q61H"),  ## KRAS Q61H
                          c("chr12","25380275","25380275","T","G","C1030","KRAS","Q61H"),  ## KRAS Q61H
                          c("chr12","25380275","25380275","T","G","C1030","KRAS","Q61H"),  ## KRAS Q61H
                          c("chr12","25380275","25380275","T","A","C1030","KRAS","Q61H"),  ## KRAS Q61H
                          c("chr12","25380275","25380275","T","G","C1037","KRAS","Q61H"))  ## KRAS Q61H
                          
## Detection of circulating tumor ... (Bettegowda et al., Science Translational Medicine, 2014)
Bettegowda_KnownMuts<-rbind(c("chr12","25398284","25398284","C","T","Patient5","KRAS","G12D"),  ## KRAS G12D
                            c("chr12","25380275","25380275","T","A","Patient16","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25398284","25398284","C","A","Patient18","KRAS","G12V"), ## KRAS G12V
                            c("chr1","115258748","115258748","C","T","Patient18","NRAS","G12S"), ## NRAS G12S
                            c("chr1","115256528","115256528","T","A","Patient19","NRAS","Q61H"), ## NRAS Q61H
                            c("chr12","25398284","25398284","C","A","Patient21","KRAS","G12V"), ## KRAS G12V
                            c("chr1","115256528","115256528","T","A","Patient21","NRAS","Q61H"), ## NRAS Q61H
                            c("chr1","115256529","115256529","T","C","Patient21","NRAS","Q61R"), ## NRAS Q61R
                            c("chr12","25398284","25398284","C","G","Patient22","KRAS","G12A"), ## KRAS G12A
                            c("chr12","25398285","25398285","C","A","Patient22","KRAS","G12C"), ## KRAS G12C
                            c("chr12","25398285","25398285","C","T","Patient22","KRAS","G12S"), ## KRAS G12S
                            c("chr12","25398284","25398284","C","T","Patient22","KRAS","G12D"), ## KRAS G12D
                            c("chr12","25398284","25398284","C","A","Patient22","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25398284","25398284","C","G","Patient24","KRAS","G12A"), ## KRAS G12A
                            c("chr1","115256530","115256530","G","T","Patient24","NRAS","Q61K"), ## NRAS Q61K
                            c("chr12","25398284","25398284","C","A","Patient26","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25398284","25398284","C","A","Patient26","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25398284","25398284","C","G","Patient1","KRAS","G12A"), ## KRAS G12A
                            c("chr12","25398285","25398285","C","A","Patient1","KRAS","G12C"), ## KRAS G12C
                            c("chr12","25398285","25398285","C","G","Patient1","KRAS","G12R"), ## KRAS G12R
                            c("chr12","25398284","25398284","C","T","Patient1","KRAS","G12D"), ## KRAS G12D
                            c("chr12","25398284","25398284","C","A","Patient1","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25380275","25380275","T","A","Patient1","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25380276","25380276","T","A","Patient1","KRAS","Q61L"), ## KRAS Q61L
                            c("chr12","25380276","25380276","T","C","Patient1","KRAS","Q61R"), ## KRAS Q61R
                            c("chr1","115256530","115256530","G","T","Patient1","NRAS","Q61K"), ## NRAS Q61K
                            c("chr1","115256528","115256528","T","G","Patient1","NRAS","Q61H"), ## NRAS Q61H
                            c("chr1","115256529","115256529","T","A","Patient1","NRAS","Q61L"), ## NRAS Q61L
                            c("chr1","115256529","115256529","T","C","Patient1","NRAS","Q61R"), ## NRAS Q61R
                            c("chr12","25398284","25398284","C","A","Patient2","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25398285","25398285","C","G","Patient4","KRAS","G12R"), ## KRAS G12R
                            c("chr12","25380275","25380275","T","A","Patient4","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25398284","25398284","C","G","Patient7","KRAS","G12A"), ## KRAS G12A
                            c("chr12","25398284","25398284","C","A","Patient7","KRAS","G12V"), ## KRAS G12V
                            c("chr1","115256528","115256528","T","A","Patient7","NRAS","Q61H"), ## NRAS Q61H
                            c("chr1","115256530","115256530","G","T","Patient7","NRAS","Q61K"), ## NRAS Q61K
                            c("chr1","115256529","115256529","T","A","Patient7","NRAS","Q61L"), ## NRAS Q61L
                            c("chr12","25398284","25398284","C","T","Patient9","KRAS","G12D"), ## KRAS G12D
                            c("chr1","115256530","115256530","G","T","Patient9","NRAS","Q61K"), ## NRAS Q61K
                            c("chr12","25398284","25398284","C","A","Patient10","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25380275","25380275","T","A","Patient10","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25380275","25380275","T","G","Patient10","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25398284","25398284","C","G","Patient12","KRAS","G12A"), ## KRAS G12A
                            c("chr12","25398285","25398285","C","A","Patient12","KRAS","G12C"), ## KRAS G12C
                            c("chr12","25380275","25380275","T","G","Patient12","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25398285","25398285","C","G","PatientBARD101","KRAS","G12R"), ## KRAS G12R
                            c("chr12","25380275","25380275","T","A","PatientBARD101","KRAS","Q61H"), ## KRAS Q61H
                            c("chr1","115256529","115256529","T","C","PatientBARD101","NRAS","Q61R"), ## NRAS Q61R
                            c("chr12","25398285","25398285","C","G","PatientBARD102","KRAS","G12R"), ## KRAS G12R
                            c("chr12","25398285","25398285","C","G","PatientBARD102","KRAS","G12R"), ## KRAS G12R
                            c("chr12","25380277","25380277","G","C","PatientBARD102","KRAS","Q61E"), ## KRAS Q61E
                            c("chr12","25398284","25398284","C","T","PatientBARD103","KRAS","G12D"), ## KRAS G12D
                            c("chr12","25398284","25398284","C","A","PatientBARD103","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25398285","25398285","C","A","CRC188","KRAS","G12C"), ## KRAS G12C
                            c("chr12","25380275","25380275","T","A","CRC189","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25380275","25380275","T","G","CRC189","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25398284","25398284","C","G","CRC190","KRAS","G12A"), ## KRAS G12A
                            c("chr12","25380275","25380275","T","G","CRC190","KRAS","Q61H"), ## KRAS Q61H
                            c("chr12","25380276","25380276","T","C","CRC190","KRAS","Q61R"), ## KRAS Q61R
                            c("chr1","115256529","115256529","T","A","CRC190","NRAS","Q61L"), ## NRAS Q61L
                            c("chr1","115256529","115256529","T","C","CRC190","NRAS","Q61R"), ## NRAS Q61R
                            c("chr12","25398284","25398284","C","T","CRC191","KRAS","G12D"), ## KRAS G12D
                            c("chr12","25398285","25398285","C","G","CRC191","KRAS","G12R"), ## KRAS G12R
                            c("chr12","25398284","25398284","C","A","CRC191","KRAS","G12V"), ## KRAS G12V
                            c("chr12","25380275","25380275","T","G","CRC191","KRAS","Q61H")) ## KRAS Q61H

##############################################################################

KRASNRASEGFR_KnownMuts<-rbind(Prospect1_KnownMuts,Prospect2_KnownMuts,Bettegowda_KnownMuts); KRASNRASEGFR_KnownMuts=KRASNRASEGFR_KnownMuts[which(KRASNRASEGFR_KnownMuts[,7]%in%c("KRAS","NRAS")),];

All.Mut.Prop<-dim(KRASNRASEGFR_KnownMuts)[1]
Q61.Prop<-dim(KRASNRASEGFR_KnownMuts[which(KRASNRASEGFR_KnownMuts[,8]=="Q61H"),])[1]/All.Mut.Prop
notQ61.Prop<-dim(KRASNRASEGFR_KnownMuts[which(KRASNRASEGFR_KnownMuts[,8]!="Q61H"),])[1]/All.Mut.Prop

## make barplot
mids<-barplot(c(notQ61.Prop,Q61.Prop),col=c("red","blue"),ylim=c(0,1),ylab="Percentage of KRAS NRAS Mutations Generated among G12/G13/Q61")
axis(1,c("G12 G13 Q61 (Not Q61H)","Q61H"),at=mids)

## fold change
text(x=mids[1]+(mids[2]-mids[1])/2,y=0.9,labels=paste0(round(notQ61.Prop/Q61.Prop,2),"-fold fewer Q61H generated"))

########################### Figure S7C-F #####################################

########################## Clinical Data #####################################

## age (years) at diagnosis of mets
AGE_AT_DIAGNOSIS=68.4
## median time (years) since diagnosis of mets
DIAGNOSIS_TO_BIOPSY=2.7
## median age (years) at BL biopsy
AGE_AT_BIOPSY=71.1
## 2x the time from diagnosis as tumours are likely there for 6 months to a few years before diagnosed
AGE_AT_TUMOUR_INITIATION=5.4

##################### Define hotspot mutations ################################

## Every KRAS/NRAS mutation considered in Figure 3A
hotspot.contexts=rbind(
  c("A[C>G]C","KRAS G12A"),
  c("C[C>A]A","KRAS G12C"),
  c("A[C>T]C","KRAS G12D"),
  c("C[C>G]A","KRAS G12R"),
  c("C[C>T]A","KRAS G12S"),
  c("A[C>A]C","KRAS G12V"),
  c("G[C>T]C","KRAS G13D"),
  c("C[C>A]A","KRAS G13C"),
  c("T[C>G]A","KRAS Q61E"),
  c("C[T>A]T","KRAS Q61H (T>A)"),
  c("C[T>G]T","KRAS Q61H (T>G)"),
  c("T[C>A]A","KRAS Q61K"),
  c("T[T>A]G","KRAS Q61L"),
  c("T[T>G]G","KRAS Q61P"),
  c("T[T>C]G","KRAS Q61R"),
  c("A[C>G]C","NRAS G12A"),
  c("C[C>A]T","NRAS G12C"),
  c("A[C>T]C","NRAS G12D"),
  c("C[C>T]T","NRAS G12S"),
  c("A[C>A]C","NRAS G12V"),
  c("A[C>T]C","NRAS G13D"),
  c("C[C>G]A","NRAS G13R"),
  c("C[T>A]T","NRAS Q61H (T>A)"),
  c("C[T>G]T","NRAS Q61H (T>G)"),
  c("T[T>A]G","NRAS Q61L"),
  c("A[C>A]A","NRAS Q61K"),
  c("T[T>C]G","NRAS Q61R"))

############### define constant and temporary signatures ####################

## signatures with HR/MSI signatures treated as clock like
signatures.constant.HRMSIclock<-c('SBS1','SBS3','SBS5','SBS6','SBS15','SBS17a','SBS40')
signatures.temp.HRMSIclock<-c('SBS17b','SBS35')

## signatures with HR/MSI signatures treated as temporary
signatures.constant.HRMSItemp<-c('SBS1','SBS5','SBS17a','SBS40')
signatures.temp.HRMSItemp<-c('SBS3','SBS6','SBS15','SBS17b','SBS35')

rates=c(1,5,10)

############### define constant and temporary signatures ####################

MutationRateAdjustment<-function(rate.acceleration,constant.signatures,acceleration.start,acceleration.end){
  
  acceleration.diff=acceleration.end-acceleration.start
  
  rate.matrix<-Signature.AssignedMutations[constant.signatures,]/(acceleration.start+acceleration.diff*rate.acceleration)
  
  return(rate.matrix)
  
}

##############################################################################
