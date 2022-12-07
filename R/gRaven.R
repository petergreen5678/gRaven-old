## gRaven interface functions

# domain set up

hugin.domain<-function()
	{
	repeat{
		dname<-paste0('.gRvd',sample(1000,1))
		if(!dname%in%ls(pattern='.gRvd',all.names=TRUE)) break
		}
	dom<-list(nodes=NULL,states=NULL,parents=NULL,cptables=NULL)
	class(dom)<-'gRv'
	attr(dom,'owner')<-NA
	assign(dname,dom,envir=.GlobalEnv)
	class(dname)<-'gRvname'
	dname
	}

clone.domain<-function(dname1)
	{
	repeat{
		dname2<-paste0('.gRvd',sample(1000,1))
		if(!dname2%in%ls(pattern='.gRvd',all.names=TRUE)) break
		}
	dom<-get(dname1,envir=.GlobalEnv)
	assign(dname2,dom,envir=.GlobalEnv)
	class(dname2)<-'gRvname'
	dname2
	}

initialize.domain<-function(dname)
	{
	dom<-get(dname,envir=.GlobalEnv)
	if(is.null(dom$net)) return()
	dom$net<-retractEvidence(dom$net)
	assign(dname,dom,envir=.GlobalEnv)
	}

# printing

print.gRv<-function(x, ...) 
	{
	owner<-attr(x,'owner')
	if(is.null(owner)) cat('orphaned ')
	cat('gRv object')
	if((!is.na(owner))&&!is.null(owner)) cat(' owned by', attr(x,'owner'),'\n')
	cat(' with slots:',names(x),'\n')
	cat(length(x$nodes),'nodes:',x$nodes,fill=60)
	cat(sum(sapply(x$parents,length)),'edges\n')
#	for(n in x$nodes) for(p in x$parents[[n]]) cat('  ',p,'->',n,'\n')
	for(n in x$nodes) if(length(x$parents[[n]])>0) cat('  ',x$parents[[n]],'->',n,'\n')
	invisible(NULL)
	}

print.gRvname<-function(x, ...) 
	{
	dom<-get(x,envir=.GlobalEnv)
	cat('owner of ',as.character(x),': ',sep='')
	cat('gRv object\n')
	cat(' with slots:',names(dom),'\n')
	invisible(NULL)
	}

# creating and modifying BN

add.node<-function (dname, n, states, subtype) 
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
    dom <- get(dname, envir = .GlobalEnv)
    if (n %in% dom$nodes) 
        stop(n, " already in ", dname, "\n")
    dom$nodes <- c(dom$nodes, n)
    dom$states <- c(dom$states, structure(list(states),names=n))
    dom$parents <- c(dom$parents, structure(list(NULL),names=n))
    assign(dname, dom, envir = .GlobalEnv)
}

add.edge<-function(dname,child,parent)
	{
	dom<-get(dname,envir=.GlobalEnv)
	if((!child%in%dom$nodes)||any(!parent%in%dom$nodes)) stop(child,'',parent,' not all already in ',dname,'\n')
	dom$parents[[child]]<-c(dom$parents[[child]],parent)
	dom$cptables[[child]]<-NULL
	assign(dname,dom,envir=.GlobalEnv)
	}

delete.node<-function (dname, n) 
{
    dom <- get(dname, envir = .GlobalEnv)
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
    assign(dname, dom, envir = .GlobalEnv)
}

delete.edge<-function (dname, n, p) 
{
    dom <- get(dname, envir = .GlobalEnv)
    pars <- dom$parents[[n]]
    pars <- pars[pars != p]
    if (length(pars) == 0) 
        dom$parents[n] <- list(NULL)
    else dom$parents[[n]] <- pars
    dom$cptables[[n]] <- NULL
    assign(dname, dom, envir = .GlobalEnv)
}

get.table<-function (dname,n) 
{
# delivers CPT as a data.frame, either by extracting it if it already exists in get(dname)$cptables, 
# or initialised with freq=1
dom<-get(dname,envir=.GlobalEnv)
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

set.table<-function(dname,n,tab=1,type='cpt')
{
if(is.data.frame(tab)) Freq<-tab$Freq else Freq<-as.vector(tab)
dom<-get(dname,envir=.GlobalEnv)
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
assign(dname,dom,envir=.GlobalEnv)
}

# compiling

compile.gRvname<-function(object, ...)
	{
	dom<-get(object,envir=.GlobalEnv)
	if(!is.null(dom$net)) warning(object," already compiled")
# if any nodes are missing cptables, provide dummy table
	for(n in setdiff(dom$nodes,names(dom$cptables))) {
		vpa<-c(n,dom$parents[[n]])
		allstates<-dom$states[vpa]
		nlev<-unlist(lapply(allstates,length))
		leng<-prod(nlev)
		dom$cptables[[n]]<-cpt(vpa,values=rep_len(1,leng),levels=allstates)
	}
	net<-grain(compileCPT(dom$cptables))
	class(net)<-c("cpt_grain","grain")
	dom$net<-net
	assign(object,dom,envir=.GlobalEnv)
	}

# evidence and propagation

set.finding<-function(dname, node, finding)
	{
	check.compiled(dname)
	dom <- get(dname, envir = .GlobalEnv)
	if (is.list(finding)) finding<-unlist(finding)
	if (length(finding) == 1) finding <- as.integer(finding == dom$states[[node]])
	cache<-dom$net$cache
	if(is.null(cache)) cache<-list()
	cache[[node]]<-finding
	dom$net$cache<-cache
	assign(dname, dom, envir = .GlobalEnv)
	}

propagate.gRvname<-function(object, ...) 
	{
	check.compiled(object)
	dom<-get(object,envir=.GlobalEnv)
	if(!is.null(dom$net$cache))
		{
		net1<-setEvidence(dom$net,evidence=dom$net$cache,propagate=TRUE)
		net1$cache<-NULL
		} else {
		net1<-propagate(dom$net)
		}
	dom$net<-net1
	assign(object,dom,envir=.GlobalEnv)
	}

get.normalization.constant<-function(dname,log=FALSE) 
	{
	dom<-get(dname,envir=.GlobalEnv)
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

# basic queries

get.nodes<-function(dname)
{
get(dname,envir=.GlobalEnv)$nodes
}

get.parents<-function(dname,n)
{
get(dname,envir=.GlobalEnv)$parents[[n]]
}

check.compiled<-function(dname)
{
	dom<-get(dname,envir=.GlobalEnv)
	if(!all(dom$nodes%in%names(dom$cptables))) {
		if(is.null(dom$cptables)) dom$cptables<-list()
		for(n in dom$nodes) if(is.null(dom$cptables[[n]])) {set.table(dname,n,1)} #; cat('set table',dname,n,'\n')}
		dom<-get(dname,envir=.GlobalEnv)
		dom$net<-NULL
	}
	if(is.null(dom$net)) {compile.gRvname(dname); cat('compiled',dname,'\n')}
}

get.belief<-function (dname, n) 
{
    dom <- get(dname, envir = .GlobalEnv)
    unlist(querygrain(dom$net, n,exclude=FALSE))
}

triangulate.gRvname<-function(object, ...) {}

compress<-function(dname) {1}

list.domains<-function(print=TRUE,clean=FALSE)
{
domains<-NULL
owners<-NULL
owned<-NULL
for(z in ls(envir = .GlobalEnv, all.names = TRUE)) if(is.character(z)){
x<-get(z)

if(is(x,"gRv"))
{
domains<-c(domains,z)
} else if(is(x,"gRvname")) {
owners<-c(owners,z)
}

if(is.list(x)&&"domains"%in%names(x)) {
xd<-x[['domains']]
for(y in names(xd)) if(is(xd[[y]],"gRvname")) 
{
owners<-c(owners,paste0(z,'$domains$',y))
}
} # end if is.list
} # end for z
if(print&&length(domains)>0) cat('domains:\n',domains,fill=60)
for(r in owners) 
{
dom<-as.character(eval(parse(text=r)))
if(print) cat(r,'owns',dom,'\n')
owned<-c(owned,dom)
}
orphans<-setdiff(domains,owned)
if(length(orphans)>0) 
	{
	if(print) cat('apparent orphans:\n',orphans,fill=60)
	if(clean) 
		{
		rm(list = orphans, envir = .GlobalEnv)
		if(print) cat('.. deleted\n')
		}
	}
invisible(domains)
} # end function
