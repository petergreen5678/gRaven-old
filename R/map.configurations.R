map.configurations<-function (dom, nodes, pmin) 
{
	propagate(dom)
    if (is.null(getOption("mapcon")) || getOption("mapcon") == 
        "sim") {
        n <- 10000
        sg <- simulate.grain(dom$net, n)
        tz <- table(sg[as.vector(nodes)])
        w <- which(as.array(tz) > pmin * n, arr.ind = TRUE)
        o <- order(tz[w], decreasing = TRUE)
        wo <- w[o, , drop = FALSE]
        sw <- wo
        for (v in nodes) sw[, v] <- dom$states[[v]][wo[, 
            v]]
        zz <- cbind(sw, tz[w[o, , drop = FALSE]]/n)
        dimnames(zz) <- list(NULL, c(nodes, "Prob"))
        as.data.frame(zz)
    }
    else {
        z <- querygrain(dom$net, nodes, "joint", exclude = FALSE)
	zz<-as.data.frame.table(z,stringsAsFactors=FALSE)
	d<-zz[c(nodes,'Freq')]
	names(d)<-c(nodes,'Prob')
for(n in nodes) storage.mode(zz[[n]])<-storage.mode(dom$states[[n]])
        o <- order(d$Prob, decreasing = TRUE)
        d <- d[o, ]
        trunc <- match(TRUE, cumsum(d$Prob) > 1 - pmin)
        structure(d[1:trunc, ], row.names = 1:trunc)
    }
}
