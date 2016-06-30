#To have this script running, copy and paste the following line at R console:
#source("../masking.R")

step <- 0
steps <- 0
masking.R <- scan("../masking_package/masking.R",what="character",quiet=TRUE)
for (i in 1:length(masking.R)) {
	ifelse(masking.R[i] == "task",steps <- steps+1,NA)
	}
steps <- steps-1

#Read configuration file ("masking.cfg").

task <- "\nReading configuration file (masking.cfg)..."
message(task)

library(seqinr)
masking.cfg <- scan("../masking_package/masking.cfg",what="character",quiet=TRUE)
charset <- masking.cfg[1:(length(masking.cfg)-3)]
num.charset <- length(charset)
agreement <- as.numeric(masking.cfg[length(masking.cfg)-2])
seqtype <- masking.cfg[length(masking.cfg)-1]
ifaa <- paste("_",seqtype,sep="")
ifelse(seqtype == "DNA",ifaa <- "",ifelse(seqtype == "AA",ifaa <- "_aa",warning(paste("Unexpected 'seqtype'! ",seqtype," was found in .cfg file, while 'AA' or 'DNA' is expected...",sep=""))))
nbchar <- as.numeric(masking.cfg[length(masking.cfg)])

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Set up names of files to work on.

	#Originals

	task <- "\nOriginal files..."
	message(task)

	files_original <- character()
	for (i in 1:num.charset) {
	files_original <- c(files_original,paste("../T-Coffee/",charset[i],ifaa,".fas",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")


	#Aliscore

	task <- "\nAliscore files..."
	message(task)

	files_Aliscore <- character()
	for (i in 1:num.charset) {
	files_Aliscore <- c(files_Aliscore,paste("../Aliscore/",charset[i],ifaa,".fasta_aln_Profile_random.txt",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#BMGE

	task <- "\nBMGE files..."
	message(task)

	files_BMGE <- character()
	for (i in 1:num.charset) {
	files_BMGE <- c(files_BMGE,paste("../BMGE/",charset[i],ifaa,".fasta_aln_BMGE.htm",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#GBlocks

	task <- "\nGBlocks files..."
	message(task)

	files_GBlocks <- character()
	for (i in 1:num.charset) {
	files_GBlocks <- c(files_GBlocks,paste("../GBlocks/",charset[i],ifaa,".fasta_aln-gb.htm",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#Noisy

	task <- "\nNoisy files..."
	message(task)

	files_Noisy <- character()
	for (i in 1:num.charset) {
	files_Noisy <- c(files_Noisy,paste("../Noisy/",charset[i],ifaa,"_sta.gr",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#T-Coffee

	task <- "\nT-Coffee files..."
	message(task)

	files_TCoffee <- character()
	for (i in 1:num.charset) {
	files_TCoffee <- c(files_TCoffee,paste("../T-Coffee/",charset[i],ifaa,".fasta_aln",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

#Prepare the final table with gene lengths extracting data from Aliscore output - contextually, read and record different softwares' opinions.

task <- "\nPreparing the final table with gene lengths extracting data from Aliscore output..."
message(task)

sites <- 0
gene <- character()
site <- numeric()
Aliscore <- numeric()
BMGE <- numeric()
GBlocks <- numeric()
Noisy <- numeric()
consensus <- numeric()

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Read and record different softwares' opinions.

task <- "\nReading and recording different softwares' opinions..."
message(task)

lengths <- numeric()
previous_sites <- 0

for (c in 1:num.charset) {
	current_Aliscore <- read.table(files_Aliscore[c],header=TRUE)
	current_sites <- length(current_Aliscore$position)
	lengths <- c(lengths,current_sites)
	max_sites <- max(c(previous_sites,current_sites),na.rm=TRUE)
	max_y <- trunc(max_sites/100)*100+100
	sites <- sites + current_sites
	current_gene <- character()
	current_gene[1:current_sites] <- charset[c]
	gene <- c(gene,current_gene)
	site <- c(site,1:current_sites)


	#Read and record Aliscore opinion.

	for (i in 1:current_sites) {
		ifelse(current_Aliscore$negative[i]<0,Aliscore <- c(Aliscore,0),Aliscore <- c(Aliscore,1))
		}

	#Read and record BMGE opinion.

	current_BMGE_selected <- read.table(file=files_BMGE[c],sep=c("\n"," ","-"),stringsAsFactors=TRUE)
	current_BMGE_selected_vector <- as.vector(current_BMGE_selected$V1[length(scan(file=files_BMGE[c],what="character",sep="\n",quiet=TRUE))-6-current_sites])
	current_BMGE_split <- strsplit(current_BMGE_selected_vector,split=" ")
	current_BMGE_starts <- numeric()
	current_BMGE_ends <- numeric()
	for (i in 1:(length(current_BMGE_split[[1]])-4)) {
		current_BMGE_starts <- c(current_BMGE_starts,as.numeric(strsplit(current_BMGE_split[[1]][i+4],split="-")[[1]][1]))
		}
	for (i in 1:(length(current_BMGE_split[[1]])-4)) {
		ifelse(is.na(as.numeric(strsplit(current_BMGE_split[[1]][i+4],split="-")[[1]][2])) == TRUE,current_BMGE_ends <- c(current_BMGE_ends,as.numeric(strsplit(current_BMGE_split[[1]][i+4],split="-")[[1]][1])),current_BMGE_ends <- c(current_BMGE_ends,as.numeric(strsplit(current_BMGE_split[[1]][i+4],split="-")[[1]][2])))
		}
	for (i in 1:current_sites) {
		wsIw <- 0
		for (j in 1:length(current_BMGE_starts)) {
			ifelse(wsIw == 1,NA,ifelse(i>current_BMGE_starts[j]-1,ifelse(i<current_BMGE_ends[j]+1,wsIw <- 1,wsIw <- 0),wsIw <- 0))
			}
		BMGE <- c(BMGE,wsIw)
		}

	#Read and record GBlocks opinion.

	current_GBlocks_selected <- read.table(file=files_GBlocks[c],sep=c("\n"," ","-"),stringsAsFactors=TRUE)
	current_GBlocks_selected_vector <- as.vector(current_GBlocks_selected$V1[length(scan(file=files_GBlocks[c],what="character",sep="\n",quiet=TRUE))-3])
	current_GBlocks_split <- strsplit(current_GBlocks_selected_vector,split="")
	current_GBlocks_starts <- numeric()
	current_GBlocks_ends_raw <- numeric()
	current_GBlocks_ends <- numeric()
	wsId <- 0
	current_number <- 0
	for (i in 9:length(current_GBlocks_split[[1]])) {
		current_GBlocks_starts_length <- length(current_GBlocks_starts)
		current_GBlocks_ends_length <- length(current_GBlocks_ends_raw)
		ifelse(current_GBlocks_split[[1]][i] == "[",wsId <- 1,ifelse(any(current_GBlocks_split[[1]][i] == " ",current_GBlocks_split[[1]][i] == "]"),ifelse(wsId == 1,current_GBlocks_starts <- c(current_GBlocks_starts,current_number),current_GBlocks_ends_raw <- c(current_GBlocks_ends_raw,current_number)),current_number <- current_number*10+as.numeric(current_GBlocks_split[[1]][i])))
		ifelse(length(current_GBlocks_starts) == current_GBlocks_starts_length,NA,wsId <- 0)
		ifelse(length(current_GBlocks_starts) == current_GBlocks_starts_length,NA,current_number <- 0)
		ifelse(length(current_GBlocks_ends_raw) == current_GBlocks_ends_length,NA,current_number <- 0)
		}
	for (i in 1:length(current_GBlocks_ends_raw)) {
		ifelse(current_GBlocks_ends_raw[i] == 0,NA,current_GBlocks_ends <- c(current_GBlocks_ends,current_GBlocks_ends_raw[i]))
		}
	for (i in 1:current_sites) {
		wsIw <- 0
		for (j in 1:length(current_GBlocks_starts)) {
			ifelse(wsIw == 1,wsIw <- 1,ifelse(i>current_GBlocks_starts[j]-1,ifelse(i<current_GBlocks_ends[j]+1,wsIw <- 1,wsIw <- 0),wsIw <- 0))
			}
		GBlocks <- c(GBlocks,wsIw)
		}

	#Read and record Noisy opinion.

	current_Noisy <- read.table(files_Noisy[c],header=FALSE,skip=14+current_sites,nrows=current_sites)
	current_Noisy_cutoff <- read.table(file="../Noisy/noisy_cmd.sh",header=FALSE,skip=2,nrows=1)$V4
	for (i in 1:current_sites) {
		ifelse(current_Noisy$V2[i]<current_Noisy_cutoff,Noisy <- c(Noisy,0),Noisy <- c(Noisy,1))
		}

	previous_sites <- current_sites

	}

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Construction of the final file.

task <- "\nConstructing the final file..."
message(task)

length(consensus) <- sites
for (i in 1:sites) {
	consensus[i] <- Aliscore[i]+BMGE[i]+GBlocks[i]+Noisy[i]
	}
masking.table <- data.frame(gene=gene,site=site,Aliscore=Aliscore,BMGE=BMGE,Gblocks=GBlocks,Noisy=Noisy,consensus=consensus)
write.table(masking.table,file="../masking.msk",quote=FALSE,sep=",",row.names=FALSE)
remove(gene,site,Aliscore,BMGE,GBlocks,Noisy,consensus)

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Compute some basic statistics about selected sites.

task <- "\nComputing some basic statistics about selected sites..."
message(task)

total_original_sites <- numeric()
for (i in 1:num.charset) {
	current_original <- read.fasta(files_original[i],seqtype=seqtype,strip.desc=TRUE)
	current_num.taxa <- length(getName(current_original))
	total_original_sites[i] <- 0
	for (j in 1:current_num.taxa) {
		total_original_sites[i] <- total_original_sites[i]+getLength(current_original[j])
		}
	}

masking.table$consensusR <- ifelse(masking.table$consensus>1,1,0)
attach(masking.table)
at_least_2 <- tapply(consensusR,gene,sum)
sum_at_least_2 <- sum(at_least_2)
percentages_over_original_2 <- numeric()
for (i in 1:num.charset) {
	current_original <- read.fasta(files_original[i],seqtype=seqtype,strip.desc=TRUE)
	current_num.taxa <- length(getName(current_original))
	percentages_over_original_2 <- c(percentages_over_original_2,((at_least_2[i]*current_num.taxa*100)/total_original_sites[i]))
	}
percentages_over_aligned_2 <- (at_least_2*100)/lengths
detach(masking.table)

masking.table$consensusR <- ifelse(masking.table$consensus>2,1,0)
attach(masking.table)
at_least_3 <- tapply(consensusR,gene,sum)
sum_at_least_3 <- sum(at_least_3)
percentages_over_original_3 <- numeric()
for (i in 1:num.charset) {
	current_original <- read.fasta(files_original[i],seqtype=seqtype,strip.desc=TRUE)
	current_num.taxa <- length(getName(current_original))
	percentages_over_original_3 <- c(percentages_over_original_3,((at_least_3[i]*current_num.taxa*100)/total_original_sites[i]))
	}
percentages_over_aligned_3 <- (at_least_3*100)/lengths
detach(masking.table)

masking.table$consensusR <- ifelse(masking.table$consensus>3,1,0)
attach(masking.table)
at_least_4 <- tapply(consensusR,gene,sum)
sum_at_least_4 <- sum(at_least_4)
percentages_over_original_4 <- numeric()
for (i in 1:num.charset) {
	current_original <- read.fasta(files_original[i],seqtype=seqtype,strip.desc=TRUE)
	current_num.taxa <- length(getName(current_original))
	percentages_over_original_4 <- c(percentages_over_original_4,((at_least_4[i]*current_num.taxa*100)/total_original_sites[i]))
	}
percentages_over_aligned_4 <- (at_least_4*100)/lengths
detach(masking.table)

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Draw basic barplots about character keeping.

task <- "\nDrawing basic barplots about character keeping..."
message(task)

vector_at_least_2 <- as.vector(at_least_2)
vector_at_least_3 <- as.vector(at_least_3)
vector_at_least_4 <- as.vector(at_least_4)
list_at_least <- list(sites_2=vector_at_least_2,sites_3=vector_at_least_3,sites_4=vector_at_least_4)
data_frame_at_least <- as.data.frame(list_at_least,row.names=charset)
vector_percentages_over_original_2 <- as.vector(percentages_over_original_2)
vector_percentages_over_original_3 <- as.vector(percentages_over_original_3)
vector_percentages_over_original_4 <- as.vector(percentages_over_original_4)
list_percentages_over_original <- list(percentages_over_original_2=vector_percentages_over_original_2,percentages_over_original_3=vector_percentages_over_original_3,percentages_over_original_4=vector_percentages_over_original_4)
data_frame_percentages_over_original <- as.data.frame(list_percentages_over_original,row.names=charset)
vector_percentages_over_aligned_2 <- as.vector(percentages_over_aligned_2)
vector_percentages_over_aligned_3 <- as.vector(percentages_over_aligned_3)
vector_percentages_over_aligned_4 <- as.vector(percentages_over_aligned_4)
list_percentages_over_aligned <- list(percentages_over_aligned_2=vector_percentages_over_aligned_2,percentages_over_aligned_3=vector_percentages_over_aligned_3,percentages_over_aligned_4=vector_percentages_over_aligned_4)
data_frame_percentages_over_aligned <- as.data.frame(list_percentages_over_aligned,row.names=charset)

par(mfrow=c(3,1))
barplot(height=t(data_frame_at_least),beside=TRUE,space=c(0,1),ylim=c(0,max_y),ylab="selected",border=NA,main="selected",col=c("seagreen3","steelblue4","khaki"),cex.names=0.5,cex.axis=0.4)
barplot(height=t(data_frame_percentages_over_original),beside=TRUE,space=c(0,1),ylim=c(0,100),ylab="% original",border=NA,main="% original",col=c("seagreen3","steelblue4","khaki"),cex.names=0.5,cex.axis=0.4)
barplot(height=t(data_frame_percentages_over_aligned),beside=TRUE,space=c(0,1),ylim=c(0,100),ylab="% aligned",border=NA,main="% aligned",col=c("seagreen3","steelblue4","khaki"),cex.names=0.5,cex.axis=0.4)
#dev.copy(pdf,"../masking.pdf")
dev.off()

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Write a boundary table of all charsets - useful for PartitionFinderProtein.

task <- "\nWriting a boundary table of all charsets - useful for PartitionFinderProtein..."
message(task)

ifelse(agreement == 2,desired <- vector_at_least_2,ifelse(agreement == 3,desired <- vector_at_least_3,ifelse(agreement == 4,desired <- vector_at_least_4,warning("Non-sense agreement!"))))
start_position <- numeric()
end_position <- numeric()
length(start_position) <- num.charset
length(end_position) <- num.charset
start_position[1] <- 1
end_position[1] <- desired[1]
for (i in 2:num.charset) {
	start_position[i] <- end_position[i-1]+1
	end_position[i] <- start_position[i]+desired[i]-1
	}
equal_to <- character()
hyphen <- character()
semicolon <- character()
length(equal_to) <- num.charset
length(hyphen) <- num.charset
length(semicolon) <- num.charset
equal_to[1:num.charset] <- " = "
hyphen[1:num.charset] <- "-"
semicolon[1:num.charset] <- ";"
boundaries.cfg <- data.frame(genes=charset,equal_to=equal_to,start=start_position,hyphen=hyphen,end=end_position,semicolon=semicolon)
write.table(boundaries.cfg,file="../boundaries.cfg",quote=FALSE,sep="",row.names=FALSE,col.names=FALSE)

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Prepare and write meta-masked alignments.

task <- "\nPreparing and writing meta-masked alignments..."
message(task)

previous_sites <- 0
kept <- character()

attach(masking.table)
for (i in 1:num.charset) {
	current_kept <- numeric()
	kept[i] <- paste("Gene ",charset[i],ifaa,":",sep="")
	current_Aliscore <- read.table(files_Aliscore[i],header=TRUE)
	current_sites <- length(current_Aliscore$position)
	for (j in 1:current_sites) {
		ifelse(consensus[j+previous_sites]>(agreement-1),current_kept <- c(current_kept,site[j+previous_sites]),NA)
		ifelse(consensus[j+previous_sites]>(agreement-1),kept[i] <- paste(kept[i]," ",j,sep=""),NA)
		}
	current_TCoffee <- read.fasta(files_TCoffee[i],seqtype=seqtype,strip.desc=TRUE)
	num.taxa <- length(getName(current_TCoffee))
	for (j in 1:num.taxa) {
		current_sequence <- character()
		current_name <- getName(current_TCoffee)[j]
		for (k in 1:length(current_kept)) {
			current_sequence <- c(current_sequence,as.character(getFrag(current_TCoffee[j],current_kept[k],current_kept[k])))
			}
		write.fasta(current_sequence,names=current_name,nbchar=nbchar,file.out=paste("./",charset[i],ifaa,"_T-Coffee_masked.fas",sep=""),open="a")
		}
	previous_sites <- previous_sites + current_sites
	}
detach(masking.table)

writeLines(kept,"../kept.msk",sep="\n")

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Concatenate meta-masked alignments.

task <- "\nConcatenating meta-masked alignments..."
message(task)

files_metamasked <- character()
previous_metamasked_num.taxa <- 0
for (i in 1:num.charset) {
	current_metamasked <- read.fasta(paste("./",charset[i],ifaa,"_T-Coffee_masked.fas",sep=""),seqtype=seqtype,strip.desc=TRUE)
	current_metamasked_num.taxa <- length(getName(current_metamasked))
	ifelse(i == 1,NA,ifelse(current_metamasked_num.taxa == previous_metamasked_num.taxa,NA,warning("Charset #",i," has ",current_metamasked_num.taxa," sites, while charset #",i-1," had ",previous_metamasked_num.taxa,"! This will lead to errors in concatenating!...")))
	previous_metamasked_num.taxa <- current_metamasked_num.taxa
	files_metamasked <- c(files_metamasked,paste("./",charset[i],ifaa,"_T-Coffee_masked.fas",sep=""))
	}

for (i in 1:current_metamasked_num.taxa) {
	current_name <- character()
	previous_name <- character()
	current_concatenated_sequence <- character()
	for (j in 1:num.charset) {
		current_metamasked <- read.fasta(files_metamasked[j],seqtype=seqtype,strip.desc=TRUE)
		current_name <- getName(current_metamasked[i])
		ifelse(i == 1,NA,ifelse(current_name == previous_name,NA,warning("Sequence #",i," in charset #",j," has a different name with respect to sequence #",i," in charset #",j-1,"! Please double-check it!")))
		previous_name <- current_name
		current_concatenated_sequence <- c(current_concatenated_sequence,getSequence(current_metamasked[i][[1]]))
		}
	write.fasta(current_concatenated_sequence,names=current_name,nbchar=nbchar,file.out="../dataset_masked.fas",open="a")
	}

step <- step+1
message("Completed step ",step,"/",steps,".\n")

message("\nEnd.\n")
