


## hopScore = function(q, s, dl, p=0.5) {
##     gd = data.table(q_idx=1:length(q), chr=as.character(seqnames(q)))
##     d = merge(gd, dl[,c("q_idx", "s_idx"),with=FALSE], by="q_idx")
##     ## flanks between q and s
##     sq = start(q[d$q_idx])
##     ss = start(s[d$s_idx])
##     flanks = GRanges(d$chr, IRanges(pmin(sq, ss), pmax(sq, ss)), s_idx=d$s_idx, q_idx=d$q_idx)
##     ## 
##     hits = data.table(data.frame(findOverlaps(q, flanks)))
##     counts = hits[,list(n=.N),by=list(subjectHits)]
##     setnames(counts, c("flank_idx", "n"))
##     counts$s_idx = flanks[counts$flank_idx]$s_idx
##     counts$q_idx = flanks[counts$flank_idx]$q_idx
##     counts$hs = 2*dgeom(counts$n-1, p)
##     setkey(counts, "q_idx", "s_idx")
##     dl[counts,hscore:=hs,]
##     return(dl)
## }
