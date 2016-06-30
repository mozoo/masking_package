masking <- read.csv(file="../masking.msk",header=TRUE,sep=",",check.names=TRUE)
length <- length(masking[[1]])
masking$code <- paste(masking[[3]],masking[[4]],masking[[5]],masking[[6]],sep="")
masking$noone[1:length] <- 0
masking$onlyAliscore[1:length] <- 0
masking$onlyBMGE[1:length] <- 0
masking$onlyGBlocks[1:length] <- 0
masking$onlyNoisy[1:length] <- 0
masking$AliscoreBMGE[1:length] <- 0
masking$AliscoreGBlocks[1:length] <- 0
masking$AliscoreNoisy[1:length] <- 0
masking$BMGEGBlocks[1:length] <- 0
masking$BMGENoisy[1:length] <- 0
masking$GBlocksNoisy[1:length] <- 0
masking$BMGEGBlocksNoisy[1:length] <- 0
masking$AliscoreGBlocksNoisy[1:length] <- 0
masking$AliscoreBMGENoisy[1:length] <- 0
masking$AliscoreBMGEGBlocks[1:length] <- 0
masking$all[1:length] <- 0
steps.10 <- 0
for (steps.10 in 0:9) {
	for (i in (trunc(length/10*(steps.10))+1):(trunc(length/10*(steps.10+1)))) {
		ifelse(masking$code[i] == "0000",masking$noone[i] <- 1,ifelse(masking$code[i] == "1000",masking$onlyAliscore[i] <- 1,ifelse(masking$code[i] == "0100",masking$onlyBMGE[i] <- 1,ifelse(masking$code[i] == "0010",masking$onlyGBlocks[i] <- 1,ifelse(masking$code[i] == "0001",masking$onlyNoisy[i] <- 1,ifelse(masking$code[i] == "1100",masking$AliscoreBMGE[i] <- 1,ifelse(masking$code[i] == "1010",masking$AliscoreGBlocks[i] <- 1,ifelse(masking$code[i] == "1001",masking$AliscoreNoisy[i] <- 1,ifelse(masking$code[i] == "0110",masking$BMGEGBlocks[i] <- 1,ifelse(masking$code[i] == "0101",masking$BMGENoisy[i] <- 1,ifelse(masking$code[i] == "0011",masking$GBlocksNoisy[i] <- 1,ifelse(masking$code[i] == "0111",masking$BMGEGBlocksNoisy[i] <- 1,ifelse(masking$code[i] == "1011",masking$AliscoreGBlocksNoisy[i] <- 1,ifelse(masking$code[i] == "1101",masking$AliscoreBMGENoisy[i] <- 1,ifelse(masking$code[i] == "1110",masking$AliscoreBMGEGBlocks[i] <- 1,masking$all[i] <- 1)))))))))))))))
		}
	message(paste((steps.10+1)*10,"% successfully coded.",sep=""))
	}
write.table(masking,file="../masking_agreement.msk",quote=FALSE,sep="\t",row.names=FALSE)
noone <- sum(masking$noone)
onlyAliscore <- sum(masking$onlyAliscore)
onlyBMGE <- sum(masking$onlyBMGE)
onlyGBlocks <- sum(masking$onlyGBlocks)
onlyNoisy <- sum(masking$onlyNoisy)
AliscoreBMGE <- sum(masking$AliscoreBMGE)
AliscoreGBlocks <- sum(masking$AliscoreGBlocks)
AliscoreNoisy <- sum(masking$AliscoreNoisy)
BMGEGBlocks <- sum(masking$BMGEGBlocks)
BMGENoisy <- sum(masking$BMGENoisy)
GBlocksNoisy <- sum(masking$GBlocksNoisy)
BMGEGBlocksNoisy <- sum(masking$BMGEGBlocksNoisy)
AliscoreGBlocksNoisy <- sum(masking$AliscoreGBlocksNoisy)
AliscoreBMGENoisy <- sum(masking$AliscoreBMGENoisy)
AliscoreBMGEGBlocks <- sum(masking$AliscoreBMGEGBlocks)
all <- sum(masking$all)
masking_agreement <- c(length,noone,onlyAliscore,onlyBMGE,onlyGBlocks,onlyNoisy,AliscoreBMGE,AliscoreGBlocks,AliscoreNoisy,BMGEGBlocks,BMGENoisy,GBlocksNoisy,BMGEGBlocksNoisy,AliscoreGBlocksNoisy,AliscoreBMGENoisy,AliscoreBMGEGBlocks,all)
masking_agreement.names <- c("length","noone","onlyAliscore","onlyBMGE","onlyGBlocks","onlyNoisy","AliscoreBMGE","AliscoreGBlocks","AliscoreNoisy","BMGEGBlocks","BMGENoisy","GBlocksNoisy","BMGEGBlocksNoisy","AliscoreGBlocksNoisy","AliscoreBMGENoisy","AliscoreBMGEGBlocks","all")
par(fin=c(5,5))
barplot(height=masking_agreement,space=0,names.arg=masking_agreement.names,ylim=c(0,(trunc(length/1000)+1)*1000),ylab="sites",border="black",main="masking_agreement",col=c("white","burlywood4","darkorange4","darkorange2","darkorange4","darkorange2","olivedrab4","olivedrab3","olivedrab4","olivedrab3","olivedrab4","olivedrab3","sandybrown","salmon3","sandybrown","salmon3","yellow"),cex.names=1,cex.axis=1,las=2)
#dev.copy(pdf,"../masking_agreement.pdf")
dev.off()
