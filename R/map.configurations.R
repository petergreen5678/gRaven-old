map.configurations<-function (dname, nodes, pmin) 
{
	propagate(dname)
    if (is.null(getOption("mapcon")) || getOption("mapcon") == 
        "sim") {
        n <- 10000
        sg <- simulate.grain(get(dname)$net, n)
        tz <- table(sg[as.vector(nodes)])
        w <- which(as.array(tz) > pmin * n, arr.ind = TRUE)
        o <- order(tz[w], decreasing = TRUE)
        wo <- w[o, , drop = FALSE]
        sw <- wo
        for (v in nodes) sw[, v] <- get(dname)$states[[v]][wo[, 
            v]]
        zz <- cbind(sw, tz[w[o, , drop = FALSE]]/n)
        dimnames(zz) <- list(NULL, c(nodes, "Prob"))
        as.data.frame(zz)
    }
    else {
        z <- querygrain(get(dname)$net, nodes, "joint", exclude = FALSE)
	zz<-as.data.frame.table(z,stringsAsFactors=FALSE)
	d<-zz[c(nodes,'Freq')]
	names(d)<-c(nodes,'Prob')
for(n in nodes) storage.mode(zz[[n]])<-storage.mode(get(dname)$states[[n]])
        o <- order(d$Prob, decreasing = TRUE)
        d <- d[o, ]
        trunc <- match(TRUE, cumsum(d$Prob) > 1 - pmin)
        structure(d[1:trunc, ], row.names = 1:trunc)
    }
}
