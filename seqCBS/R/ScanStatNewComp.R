ScanStatNewComp <-
function(combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, max.win, statistic) {
	if(statistic=="rabinowitz") {
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		newRes = .Call(C_ScanStatNewCompRabinC, combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, max.win)
	}
	else if(statistic=="normal") {
		SijFactorN = p*(1-p)
		newRes = .Call(C_ScanStatNewCompNormalC, combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorN, p, nTotal, grid.cur, max.win)
	}
	else if(statistic=="binomial") {
		newRes = .Call(C_ScanStatNewCompBinomC, combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, max.win)
	}
	else {
		print("Name of Statistic is not Recognized, Rabinowitz used")
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		newRes = .Call(C_ScanStatNewCompRabinC, combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, max.win)
	}
	return(newRes)
}

