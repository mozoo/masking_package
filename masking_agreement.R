masking <- read.csv(file="../masking.msk",header=TRUE,sep=",",check.names=TRUE)
length <- length(masking[[1]])

masking$code <- paste(masking[[3]],masking[[4]],masking[[5]],masking[[6]],masking[[7]],sep="")

masking$noone[1:length] <- 0

masking$onlyAliscore[1:length] <- 0
masking$onlyBMGE[1:length] <- 0
masking$onlyGBlocks[1:length] <- 0
masking$onlyNoisy[1:length] <- 0
masking$onlyZorro[1:length] <- 0

masking$AliscoreBMGE[1:length] <- 0
masking$AliscoreGBlocks[1:length] <- 0
masking$AliscoreNoisy[1:length] <- 0
masking$AliscoreZorro[1:length] <- 0
masking$BMGEGBlocks[1:length] <- 0
masking$BMGENoisy[1:length] <- 0
masking$BMGEZorro[1:length] <- 0
masking$GBlocksNoisy[1:length] <- 0
masking$GBlocksZorro[1:length] <- 0
masking$NoisyZorro[1:length] <- 0

masking$GBlocksNoisyZorro[1:length] <- 0
masking$BMGENoisyZorro[1:length] <- 0
masking$BMGEGBlocksZorro[1:length] <- 0
masking$BMGEGBlocksNoisy[1:length] <- 0
masking$AliscoreNoisyZorro[1:length] <- 0
masking$AliscoreGBlocksZorro[1:length] <- 0
masking$AliscoreGBlocksNoisy[1:length] <- 0
masking$AliscoreBMGEZorro[1:length] <- 0
masking$AliscoreBMGENoisy[1:length] <- 0
masking$AliscoreBMGEGBlocks[1:length] <- 0

masking$BMGEGBlocksNoisyZorro[1:length] <- 0
masking$AliscoreGBlocksNoisyZorro[1:length] <- 0
masking$AliscoreBMGENoisyZorro[1:length] <- 0
masking$AliscoreBMGEGBlocksZorro[1:length] <- 0
masking$AliscoreBMGEGBlocksNoisy[1:length] <- 0

masking$all[1:length] <- 0

steps.10 <- 0
for (i in 1:length) {
	if (masking$code[i] == "00000") masking$noone[i] <- 1
	else if (masking$code[i] == "10000") masking$onlyAliscore[i] <- 1
	else if (masking$code[i] == "01000") masking$onlyBMGE[i] <- 1
	else if (masking$code[i] == "00100") masking$onlyGBlocks[i] <- 1
	else if (masking$code[i] == "00010") masking$onlyNoisy[i] <- 1
	else if (masking$code[i] == "00001") masking$onlyZorro[i] <- 1
	else if (masking$code[i] == "11000") masking$AliscoreBMGE[i] <- 1
	else if (masking$code[i] == "10100") masking$AliscoreGBlocks[i] <- 1
	else if (masking$code[i] == "10010") masking$AliscoreNoisy[i] <- 1
	else if (masking$code[i] == "10001") masking$AliscoreZorro[i] <- 1
	else if (masking$code[i] == "01100") masking$BMGEGBlocks[i] <- 1
	else if (masking$code[i] == "01010") masking$BMGENoisy[i] <- 1
	else if (masking$code[i] == "01001") masking$BMGEZorro[i] <- 1
	else if (masking$code[i] == "00110") masking$GBlocksNoisy[i] <- 1
	else if (masking$code[i] == "00101") masking$GBlocksZorro[i] <- 1
	else if (masking$code[i] == "00011") masking$NoisyZorro[i] <- 1
	else if (masking$code[i] == "00111") masking$GBlocksNoisyZorro[i] <- 1
	else if (masking$code[i] == "01011") masking$BMGENoisyZorro[i] <- 1
	else if (masking$code[i] == "01101") masking$BMGEGBlocksZorro[i] <- 1
	else if (masking$code[i] == "01110") masking$BMGEGBlocksNoisy[i] <- 1
	else if (masking$code[i] == "10011") masking$AliscoreNoisyZorro[i] <- 1
	else if (masking$code[i] == "10101") masking$AliscoreGBlocksZorro[i] <- 1
	else if (masking$code[i] == "10110") masking$AliscoreGBlocksNoisy[i] <- 1
	else if (masking$code[i] == "11001") masking$AliscoreBMGEZorro[i] <- 1
	else if (masking$code[i] == "11010") masking$AliscoreBMGENoisy[i] <- 1
	else if (masking$code[i] == "11100") masking$AliscoreBMGEGBlocks[i] <- 1
	else if (masking$code[i] == "01111") masking$BMGEGBlocksNoisyZorro[i] <- 1
	else if (masking$code[i] == "10111") masking$AliscoreGBlocksNoisyZorro[i] <- 1
	else if (masking$code[i] == "11011") masking$AliscoreBMGENoisyZorro[i] <- 1
	else if (masking$code[i] == "11101") masking$AliscoreBMGEGBlocksZorro[i] <- 1
	else if (masking$code[i] == "11110") masking$AliscoreBMGEGBlocksNoisy[i] <- 1
	else if (masking$code[i] == "11111") masking$all[i] <- 1
	if (floor(i/length*10) == steps.10+1) {
		message(paste((steps.10+1)*10,"% successfully coded.",sep=""))
		steps.10 <- steps.10+1
		}
	}

write.table(masking,file="../masking_agreement.msk",quote=FALSE,sep="\t",row.names=FALSE)
noone <- sum(masking$noone)
onlyAliscore <- sum(masking$onlyAliscore)
onlyBMGE <- sum(masking$onlyBMGE)
onlyGBlocks <- sum(masking$onlyGBlocks)
onlyNoisy <- sum(masking$onlyNoisy)
onlyZorro <- sum(masking$onlyZorro)
AliscoreBMGE <- sum(masking$AliscoreBMGE)
AliscoreGBlocks <- sum(masking$AliscoreGBlocks)
AliscoreNoisy <- sum(masking$AliscoreNoisy)
AliscoreZorro <- sum(masking$AliscoreZorro)
BMGEGBlocks <- sum(masking$BMGEGBlocks)
BMGENoisy <- sum(masking$BMGENoisy)
BMGEZorro <- sum(masking$BMGEZorro)
GBlocksNoisy <- sum(masking$GBlocksNoisy)
GBlocksZorro <- sum(masking$GBlocksZorro)
NoisyZorro <- sum(masking$NoisyZorro)
GBlocksNoisyZorro <- sum(masking$GBlocksNoisyZorro)
BMGENoisyZorro <- sum(masking$BMGENoisyZorro)
BMGEGBlocksZorro <- sum(masking$BMGEGBlocksZorro)
BMGEGBlocksNoisy <- sum(masking$BMGEGBlocksNoisy)
AliscoreNoisyZorro <- sum(masking$AliscoreNoisyZorro)
AliscoreGBlocksZorro <- sum(masking$AliscoreGBlocksZorro)
AliscoreGBlocksNoisy <- sum(masking$AliscoreGBlocksNoisy)
AliscoreBMGEZorro <- sum(masking$AliscoreBMGEZorro)
AliscoreBMGENoisy <- sum(masking$AliscoreBMGENoisy)
AliscoreBMGEGBlocks <- sum(masking$AliscoreBMGEGBlocks)
BMGEGBlocksNoisyZorro <- sum(masking$BMGEGBlocksNoisyZorro)
AliscoreGBlocksNoisyZorro <- sum(masking$AliscoreGBlocksNoisyZorro)
AliscoreBMGENoisyZorro <- sum(masking$AliscoreBMGENoisyZorro)
AliscoreBMGEGBlocksZorro <- sum(masking$AliscoreBMGEGBlocksZorro)
AliscoreBMGEGBlocksNoisy <- sum(masking$AliscoreBMGEGBlocksNoisy)
all <- sum(masking$all)
masking_agreement <- c(length,noone,onlyAliscore,onlyBMGE,onlyGBlocks,onlyNoisy,onlyZorro,AliscoreBMGE,AliscoreGBlocks,AliscoreNoisy,AliscoreZorro,BMGEGBlocks,BMGENoisy,BMGEZorro,GBlocksNoisy,GBlocksZorro,NoisyZorro,GBlocksNoisyZorro,BMGENoisyZorro,BMGEGBlocksZorro,BMGEGBlocksNoisy,AliscoreNoisyZorro,AliscoreGBlocksZorro,AliscoreGBlocksNoisy,AliscoreBMGEZorro,AliscoreBMGENoisy,AliscoreBMGEGBlocks,BMGEGBlocksNoisyZorro,AliscoreGBlocksNoisyZorro,AliscoreBMGENoisyZorro,AliscoreBMGEGBlocksZorro,AliscoreBMGEGBlocksNoisy,all)
masking_agreement.names <- c("length","noone","onlyAliscore","onlyBMGE","onlyGBlocks","onlyNoisy","onlyZorro","AliscoreBMGE","AliscoreGBlocks","AliscoreNoisy","AliscoreZorro","BMGEGBlocks","BMGENoisy","BMGEZorro","GBlocksNoisy","GBlocksZorro","NoisyZorro","GBlocksNoisyZorro","BMGENoisyZorro","BMGEGBlocksZorro","BMGEGBlocksNoisy","AliscoreNoisyZorro","AliscoreGBlocksZorro","AliscoreGBlocksNoisy","AliscoreBMGEZorro","AliscoreBMGENoisy","AliscoreBMGEGBlocks","BMGEGBlocksNoisyZorro","AliscoreGBlocksNoisyZorro","AliscoreBMGENoisyZorro","AliscoreBMGEGBlocksZorro","AliscoreBMGEGBlocksNoisy","all")
par(fin=c(5,5))
barplot(height=masking_agreement,space=0,names.arg=masking_agreement.names,ylim=c(0,(trunc(length/1000)+1)*1000),ylab="sites",border="black",main="masking_agreement",col=sample(c(1:657),32),cex.names=0.5,cex.axis=1,las=2)
dev.copy(pdf,"../masking_agreement.pdf")
dev.off()
