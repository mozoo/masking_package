#To have this script running, copy and paste the following line at R console:
#source("../masking.R")

library(seqinr)
num.softwares <- 5
step <- 0
steps <- 0
masking.R <- scan("../masking_package/masking.R",what="character",quiet=TRUE)
for (i in 1:length(masking.R)) {
	ifelse(masking.R[i] == "task",steps <- steps+1,NA)
	}
steps <- steps-1
agreement <- ceiling(num.softwares/2)
seqtype <- "AA"
nbchar= 60
suffix <- ""
Zorro_cutoff <- 4

#Read configuration file ("masking.cfg").

task <- "\nReading configuration file (masking.cfg)..."
message(task)

masking.cfg <- scan("../masking_package/masking.cfg",what="character",quiet=TRUE,blank.lines.skip=FALSE)
for (l in 1:length(masking.cfg)) {
	if (masking.cfg[l] == "###CHARSETS###") {
		charset.start <- l+1
		m <- l+1
		charset.end.found <- FALSE
		while (charset.end.found == FALSE) {
			m <- m+1
			if (length(strsplit(masking.cfg[m],split="")[[1]]) > 2) {
				if (paste(strsplit(masking.cfg[m],split="")[[1]][1],strsplit(masking.cfg[m],split="")[[1]][2],strsplit(masking.cfg[m],split="")[[1]][3],sep="") == "###") {
					charset.end.found <- TRUE
					charset.end <- m-1
					}
				}
			if (m == length(masking.cfg)) {
				charset.end.found <- TRUE
				charset.end <- length(masking.cfg)
				}
			}
		charset <- masking.cfg[charset.start:charset.end]
		charset <- charset[charset != ""]
		num.charset <- length(charset)
		if (num.charset == 0) stop("No charsets were detected. Please check the configuration file!")
		}
	else if (masking.cfg[l] == "###AGREEMENT###") agreement <- as.numeric(masking.cfg[l+1])
	if (agreement == 0 || agreement > num.softwares) stop("Wrong agreement setting. Please check the configuration file!")
	else if (masking.cfg[l] == "###DATATYPE###") seqtype <- masking.cfg[l+1]
	else if (masking.cfg[l] == "###LINELENGTH###") nbchar <- as.numeric(masking.cfg[l+1])
	else if (masking.cfg[l] == "###SUFFIX###") suffix <- masking.cfg[l+1]
	else if (masking.cfg[l] == "###ZORROCUTOFF###") zorro_cutoff <- as.numeric(masking.cfg[l+1])
	}

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Set up names of files to work on.

	#Originals

	task <- "\nOriginal files..."
	message(task)

	files_original <- character()
	for (i in 1:num.charset) {
	files_original <- c(files_original,paste("../T-Coffee/",charset[i],suffix,".fas",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")


	#Aliscore

	task <- "\nAliscore files..."
	message(task)

	files_Aliscore <- character()
	for (i in 1:num.charset) {
	files_Aliscore <- c(files_Aliscore,paste("../Aliscore/",charset[i],suffix,".fasta_aln_Profile_random.txt",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#BMGE

	task <- "\nBMGE files..."
	message(task)

	files_BMGE <- character()
	for (i in 1:num.charset) {
	files_BMGE <- c(files_BMGE,paste("../BMGE/",charset[i],suffix,".fasta_aln_BMGE.htm",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#GBlocks

	task <- "\nGBlocks files..."
	message(task)

	files_GBlocks <- character()
	for (i in 1:num.charset) {
	files_GBlocks <- c(files_GBlocks,paste("../GBlocks/",charset[i],suffix,".fasta_aln-gb.htm",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#Noisy

	task <- "\nNoisy files..."
	message(task)

	files_Noisy <- character()
	for (i in 1:num.charset) {
	files_Noisy <- c(files_Noisy,paste("../Noisy/",charset[i],suffix,"_sta.gr",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#T-Coffee

	task <- "\nT-Coffee files..."
	message(task)

	files_TCoffee <- character()
	for (i in 1:num.charset) {
	files_TCoffee <- c(files_TCoffee,paste("../T-Coffee/",charset[i],suffix,".fasta_aln",sep=""))
		}

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	#Zorro

	task <- "\nZorro files..."
	message(task)

	files_Zorro <- character()
	for (i in 1:num.charset) {
	files_Zorro <- c(files_Zorro,paste("../Zorro/",charset[i],suffix,".fasta_aln.zorro",sep=""))
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
Zorro <- numeric()
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

	#Read and record Zorro opinion.

	current_Zorro <- as.numeric(scan(files_Zorro[c],what="character",quiet=TRUE))
	for (i in 1:length(current_Zorro)) ifelse(current_Zorro[i] < Zorro_cutoff,Zorro <- c(Zorro,0),Zorro <- c(Zorro,1))

	previous_sites <- current_sites

	}

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Construction of the final file.

task <- "\nConstructing the final file..."
message(task)

length(consensus) <- sites
for (i in 1:sites) {
	consensus[i] <- Aliscore[i]+BMGE[i]+GBlocks[i]+Noisy[i]+Zorro[i]
	}
masking.table <- data.frame(gene=gene,site=site,Aliscore=Aliscore,BMGE=BMGE,Gblocks=GBlocks,Noisy=Noisy,Zorro=Zorro,consensus=consensus)
write.table(masking.table,file="../masking.msk",quote=FALSE,sep=",",row.names=FALSE)
remove(gene,site,Aliscore,BMGE,GBlocks,Noisy,Zorro,consensus)

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

at_least <- numeric()
sums <- numeric()
percentages_over_original <- numeric()
percentages_over_aligned <- numeric()
for (s in 1:num.softwares) {
	masking.table$consensusR <- ifelse(masking.table$consensus>s-1,1,0)
	attach(masking.table)
	at_least <- as.vector(c(at_least,tapply(consensusR,gene,sum)))
	sums <- c(sums,sum(at_least[((s-1)*num.charset+1):(s*num.charset)]))
	for (i in 1:num.charset) {
		current_original <- read.fasta(files_original[i],seqtype=seqtype,strip.desc=TRUE)
		current_num.taxa <- length(getName(current_original))
		percentages_over_original <- c(percentages_over_original,((at_least[(s-1)*num.charset+i]*current_num.taxa*100)/total_original_sites[i]))
		}
	percentages_over_aligned <- c(percentages_over_aligned,(at_least[((s-1)*num.charset+1):(s*num.charset)]*100)/lengths)
	detach(masking.table)
	}

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Draw basic barplots about character keeping.

task <- "\nDrawing basic barplots about character keeping..."
message(task)

data_frame_at_least <- as.data.frame(matrix(at_least,nrow=num.charset,ncol=num.softwares))
data_frame_percentages_over_original <- as.data.frame(matrix(percentages_over_original,nrow=num.charset,ncol=num.softwares))
data_frame_percentages_over_aligned <- as.data.frame(matrix(percentages_over_aligned,nrow=num.charset,ncol=num.softwares))
par(mfrow=c(3,1))
barplot(height=t(data_frame_at_least),beside=TRUE,space=c(0,1),ylim=c(0,max_y),names.arg=charset,ylab="selected",border=NA,main="selected",col=c("seagreen3","steelblue4"),cex.names=0.5,cex.axis=0.4)
barplot(height=t(data_frame_percentages_over_original),beside=TRUE,space=c(0,1),ylim=c(0,100),names.arg=charset,ylab="% original",border=NA,main="% original",col=c("seagreen3","steelblue4"),cex.names=0.5,cex.axis=0.4)
barplot(height=t(data_frame_percentages_over_aligned),beside=TRUE,space=c(0,1),ylim=c(0,100),names.arg=charset,ylab="% aligned",border=NA,main="% aligned",col=c("seagreen3","steelblue4"),cex.names=0.5,cex.axis=0.4)
dev.copy(pdf,"../masking.pdf")
dev.off()

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Write a boundary table of all charsets - useful for PartitionFinderProtein.

task <- "\nWriting a PartitionFinder-formatted boundary table of all charsets..."
message(task)

desired <- data_frame_at_least[[agreement]]
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
equal_to <- rep(" = ",num.charset)
hyphen <- rep("-",num.charset)
semicolon <- rep(";",num.charset)
boundaries.cfg <- data.frame(genes=charset,equal_to=equal_to,start=start_position,hyphen=hyphen,end=end_position,semicolon=semicolon)
write.table(boundaries.cfg,file="../boundaries.cfg",quote=FALSE,sep="",row.names=FALSE,col.names=FALSE)

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#Prepare and write meta-masked alignments.

task <- "\nPreparing and writing meta-masked alignments..."
message(task)

previous_sites <- 0
kept <- character()
files_metamasked <- character()

attach(masking.table)
for (i in 1:num.charset) {
	files_metamasked <- c(files_metamasked,paste("./",charset[i],suffix,"_T-Coffee_masked.fas",sep=""))
	current_kept <- numeric()
	kept[i] <- paste("Gene ",charset[i],suffix,":",sep="")
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
		write.fasta(current_sequence,names=current_name,nbchar=nbchar,file.out=files_metamasked[i],open="a")
		}
	previous_sites <- previous_sites + current_sites
	}
detach(masking.table)

writeLines(kept,"../kept.msk",sep="\n")

step <- step+1
message("Completed step ",step,"/",steps,".\n")

#If necessary, concatenate meta-masked alignments.

if (num.charset == 1) {
	write.fasta(read.fasta(files_metamasked,seqtype=seqtype,strip.desc=TRUE),names=getName(read.fasta(files_metamasked,seqtype=seqtype,strip.desc=TRUE)),nbchar=nbchar,file.out="../dataset_masked.fas",open="w")
	} else {

	task <- "\nConcatenating meta-masked alignments..."
	message(task)

	first_metamasked <- read.fasta(paste("./",charset[1],suffix,"_T-Coffee_masked.fas",sep=""),seqtype=seqtype,strip.desc=TRUE)
	found_taxa <- getName(first_metamasked)
	found_sequences <- getSequence(first_metamasked)
	current_metamasked_length <- length(found_sequences[[1]])
	charset_nums <- length(found_taxa)
	for (j in 2:num.charset) {
		current_metamasked <- read.fasta(files_metamasked[j],seqtype=seqtype,strip.desc=TRUE)
		for (k in 1:length(found_taxa)) {
			if (found_taxa[k] %in% getName(current_metamasked)) {
				found_sequences[[k]] <- c(found_sequences[[k]],getSequence(current_metamasked)[[c(1:length(current_metamasked))[getName(current_metamasked) == found_taxa[k]]]])
				}
				else {
				found_sequences[[k]] <- c(found_sequences[[k]],rep("-",length(getSequence(current_metamasked)[[1]])))
				}
			}
		for (k in 1:length(getName(current_metamasked))) {
			if (getName(current_metamasked)[k] %in% found_taxa) {
				}
				else {
				found_taxa <- c(found_taxa,getName(current_metamasked)[k])
				found_sequences[length(found_taxa)] <- c(rep("-",current_metamasked_length),getSequence(current_metamasked)[[k]])
				}
			}
		current_metamasked_length <- current_metamasked_length+length(getSequence(current_metamasked)[[1]])
		charset_nums <- c(charset_nums,length(current_metamasked))
		}
	write.fasta(found_sequences,names=found_taxa,nbchar=nbchar,file.out="../dataset_masked.fas",open="w")
	charset_nums_data_frame <- data.frame(charset=charset,taxa=charset_nums)
	write.table(charset_nums_data_frame,file="../charset_nums.txt",quote=FALSE,sep="\t",row.names=FALSE)

	step <- step+1
	message("Completed step ",step,"/",steps,".\n")

	}

message("\nEnd.\n")
