##################################
#Function to calculate Scores as matrix
#or single sample
##################################

calcScoreMat <- function(myMat=NULL, mySetsUp=NULL)
{
	getScoreSet <- function(x, myMat=myMat)
	{	
		return(colMeans(myMat[x,]))
	}

	myMatUp <- data.frame(lapply(mySetsUp,FUN=getScoreSet, myMat));
	return(myMatUp)
}


