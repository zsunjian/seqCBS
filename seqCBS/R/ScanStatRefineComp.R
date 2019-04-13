ScanStatRefineComp <-
function(combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, grid.LR, max.win, statistic) {
	if(statistic=="rabinowitz") {
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		#refineRes = ScanStatRefineCompRabin(combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, grid.L, grid.R, max.win)
		refineRes = .Call(C_ScanStatRefineCompRabinC, combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, grid.LR, max.win)
	}
	else if(statistic=="normal") {
		SijFactorN = p*(1-p)
		refineRes = .Call(C_ScanStatRefineCompNormalC, combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorN, p, nTotal, grid.cur, grid.LR, max.win)
	}
	else if(statistic=="binomial") {
		refineRes = .Call(C_ScanStatRefineCompBinomC, combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, grid.LR, max.win)
	}
	else {
		print("Name of Statistic is not Recognized, Rabinowitz used")
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		refineRes = .Call(C_ScanStatRefineCompRabinC, combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, grid.LR, max.win)
	}
	return(refineRes)
}

