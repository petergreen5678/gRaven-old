print.gRv<-function(x, ...)
{
cat('domain with slots:',names(x),'\n')
cat(length(x$nodes),'nodes:',x$nodes,fill=60)
cat(sum(sapply(x$parents,length)),'edges\n')
for(n in x$nodes) if(length(x$parents[[n]])>0) cat('  ',x$parents[[n]],'->',n,'\n')
invisible(NULL)
}


hugin.domain<-function () 
{
e <- rlang::env(
   nodes=NULL,
   states=NULL,
   parents=NULL,
   cptables=NULL,
 )
class(e)<-c("gRv","environment")
e
}

clone.domain<-function(dom1)
	{
	dom2<-rlang::env_clone(dom1)
	class(dom2)<-c("gRv","environment")
	dom2
	}

initialize.domain<-function(dom)
	{
	if(is.null(dom$net)) return()
	dom$net<-retractEvidence(dom$net)
	}


add.node<-function (dom, n, states, subtype) 
{
    if (missing(states)) {
        if (subtype == "boolean") 
            states <- c(0, 1)
    }
    else {
        if (missing(subtype)) 
            subtype <- switch(mode(states), character = "labeled", 
                numeric = "numbered", logical = "boolean")
    }
    if (n %in% dom$nodes) 
        stop(n, " already in ", dom, "\n")
    dom$nodes <- c(dom$nodes, n)
    dom$states <- c(dom$states, structure(list(states),names=n))
    dom$parents <- c(dom$parents, structure(list(NULL),names=n))
}

add.edge<-function(dom,child,parent)
{
if((!child%in%dom$nodes)||any(!parent%in%dom$nodes)) stop(child,'',parent,' not all already in ',dom,'\n')
dom$parents[[child]]<-c(dom$parents[[child]],parent)
dom$cptables[[child]]<-NULL
}

delete.node<-function (dom, n) 
{
    dom$nodes <- dom$nodes[dom$nodes != n]
    dom$states[[n]] <- NULL
    dom$parents[[n]] <- NULL
    dom$cptables[[n]] <- NULL
	for(m in dom$nodes) if(n%in%dom$parents[[m]]) 
	{
	pars<-dom$parents[[m]]
	pars<-pars[pars!=n]
	if(length(pars)==0) dom$parents[m]<-list(NULL) else dom$parents[[m]]<-pars
	dom$cptables[[m]]<-NULL
	}
}

delete.edge<-function (dom, n, p) 
{
    pars <- dom$parents[[n]]
    pars <- pars[pars != p]
    if (length(pars) == 0) 
        dom$parents[n] <- list(NULL)
    else dom$parents[[n]] <- pars
    dom$cptables[[n]] <- NULL
}

get.table<-function (dom,n) 
{
# delivers CPT as a data.frame, either by extracting it if it already exists in dom$cptables, 
# or initialised with freq=1
z<-dom$cptables
if(is.null(z)||is.null(z[[n]])) 
{
        Freq<-1
        vpa<-c(n,dom$parents[[n]])
} else {
        Freq<-z[[n]]
        vpa<-names(dimnames(Freq))
}
allstates<-dom$states[vpa]
nlev<-unlist(lapply(allstates,length))
k<-cumprod(nlev)
lk<-length(k); length<-k[lk]; k<-c(1,k[-lk])
for(i in 1:lk) 
{
        w<-allstates[[i]][as.numeric(gl(nlev[i],k[i],length))]
        if(i==1) df<-data.frame(w) else df<-cbind(df,w)
}
df<-cbind(df,as.vector(Freq))
names(df)<-c(vpa,'Freq')
df
}

set.table<-function(dom,n,tab=1,type='cpt')
{
if(is.data.frame(tab)) Freq<-tab$Freq else Freq<-as.vector(tab)
if(is.null(dom$net))
{
	if(is.null(dom$cptables)||is.null(dom$cptables[[n]]))
	{
		vpa<-c(n,dom$parents[[n]])
		allstates<-dom$states[vpa]
		nlev<-unlist(lapply(allstates,length))
		leng<-prod(nlev)
		dom$cptables[[n]]<-cpt(vpa,values=rep_len(Freq,leng),levels=allstates)
	} else {
		dom$cptables[[n]][]<-Freq
	}
} else {
	dom$net<-replaceCPT(dom$net,structure(list(Freq),names=n))
}
}

compile.gRv<-function(object, ...)
	{
	if(!is.null(object$net)) warning(object," already compiled")
# if any nodes are missing cptables, provide dummy table
	for(n in setdiff(object$nodes,names(object$cptables))) {
		vpa<-c(n,object$parents[[n]])
		allstates<-object$states[vpa]
		nlev<-unlist(lapply(allstates,length))
		leng<-prod(nlev)
		object$cptables[[n]]<-cpt(vpa,values=rep_len(1,leng),levels=allstates)
	}
	net<-grain(compileCPT(object$cptables))
	class(net)<-c("cpt_grain","grain")
	object$net<-net
	}

check.compiled<-function(dom)
{
	if(!all(dom$nodes%in%names(dom$cptables))) {
		if(is.null(dom$cptables)) dom$cptables<-list()
		for(n in dom$nodes) if(is.null(dom$cptables[[n]])) {set.table(dom,n,1)} #; cat('set table',dom,n,'\n')}
		dom<-get(dom,envir=.GlobalEnv)
		dom$net<-NULL
	}
	if(is.null(dom$net)) {compile.gRv(dom); cat('compiled',dom,'\n')}
}


set.finding<-function(dom, node, finding)
	{
	check.compiled(dom)
	if (is.list(finding)) finding<-unlist(finding)
	if (length(finding) == 1) finding <- as.integer(finding == dom$states[[node]])
	cache<-dom$net$cache
	if(is.null(cache)) cache<-list()
	cache[[node]]<-finding
	dom$net$cache<-cache
	}

get.belief<-function (dom, n) 
{
    unlist(querygrain(dom$net, n,exclude=FALSE))
}

propagate.gRv<-function(object, ...) 
	{
	check.compiled(object)
	if(!is.null(object$net$cache))
		{
		net1<-setEvidence(object$net,evidence=object$net$cache,propagate=TRUE)
		net1$cache<-NULL
		} else {
		net1<-propagate(object$net)
		}
	object$net<-net1
	}

get.normalization.constant<-function(dom,log=FALSE) 
	{
	if(!is.null(dom$net$cache))
		{
		if(!is.null(dom$net$evidence))
			{
			e<-dom$net$evidence$evi_weight
			for(i in 1:length(e))
				{
				n<-names(dimnames(e[[i]]))
				dom$net$cache[[n]]<-as.vector(e[[i]])
				}
			dom$net$evidence<-NULL
			}
		p<-pEvidence(dom$net,evidence=dom$net$cache)		
		} else {
		if(dom$net$isPropagated) 
			p<-pEvidence(dom$net) else
			{
			net1<-propagate(dom$net) 
			p<-pEvidence(net1)
			}
		} 
	if(log) log(p) else p
	}

get.nodes<-function(dom)
{
dom$nodes
}

get.parents<-function(dom,n)
{
dom$parents[[n]]
}

triangulate.gRv<-function(object, ...) {}

compress<-function(dom) {1}


list.domains<-function(print=TRUE)
{
domains<-NULL
for(x in ls(all.names=TRUE,envir=.GlobalEnv))
	if(is(get(x),'gRv')) domains<-c(domains,x)
if(print) cat('domains:',domains,fill=60)
invisible(domains)
}


