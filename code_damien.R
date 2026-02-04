##
## D.M. de Vienne
## Sept 2024
## 
## Set of functions used for computing INDUCED TRANSFERS 
## from a species tree and a collection of gene trees
## Fonctions for conversion to Newick and other tree manipulations
## are also available.
##

require(ape)
require(phangorn)


drawtransfers<-function(events, tr, ext) {
	rootedge<-tr$root.edge
	if (is.null(rootedge)) rootedge<-0
	XY<-plotPhyloCoor(tr, direction="up")	
	plot(tr,  cex=0.7, direction="up", edge.color = "grey")
	#
	NG<-setdiff(1:(Ntip(tr)+Nnode(tr)), ext)
	ED<-tr$edge[match(NG, tr$edge[,2]),]
	ED<-ED[!is.na(ED[,1]),]
	apply(ED, 1, function(x,xy) segments(xy[x[2],1], xy[x[1],2], xy[x[2],1], xy[x[2],2], col="red"), xy=XY)
	apply(ED, 1, function(x,xy) segments(xy[x[2],1], xy[x[1],2], xy[x[1],1], xy[x[1],2], col="red"), xy=XY)

#	points(XY)
	nodelabels(cex=0.5)
#	tiplabels(cex=0.5)
#	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	transferdraw<-function(arr, X) {
		ti<-as.numeric(arr$TIME)
		from<-as.numeric(arr$FROM)
		to<-as.numeric(arr$TO)
		yfromto<-ti-rootedge
#		xfromto<-ti-rootedge
		xfrom<-XY[from,1]
		xto<-XY[to,1]
#		yfrom<-lastPP$yy[from]
#		yto<-lastPP$yy[to]
#		arrows(xfromto,yfrom,xfromto,yto, col="blue")
		arrows(xfrom,yfromto,xto,yfromto, col="blue", length=0.1)
		start<-Ntip(tr)+Nnode(tr)+X
		text(xfrom, yfromto, start, cex=0.6, adj=0)
		text(xto, yfromto, -start, cex=0.6, adj=0)
	}
	sapply(1:nrow(events), function(x, ev) transferdraw(ev[x,], x), ev=events) 
#	apply(events, 1, transferdraw)
#	axis(1)
}

### THE DATA PREPARATION

GetExtendedEdgeMatrix<-function(tr, events) {
	EDGE<-tr$edge
	currnode<-Ntip(tr)+Nnode(tr)+1 #start of node numbering
	for (i in 1:nrow(events)) {
		#recipient
		whereedgeaffected<-which(EDGE[,2]==events$TO[i])
		edgeaffected<-EDGE[whereedgeaffected,]
		newedg<-rbind(c(edgeaffected[1],-currnode), c(-currnode,edgeaffected[2]))
		EDGE<-EDGE[-whereedgeaffected,] #remove old edge
		EDGE<-rbind(EDGE, newedg)
		#donor
		whereedgeaffected<-which(EDGE[,2]==events$FROM[i])
		edgeaffected<-EDGE[whereedgeaffected,]
		newedg<-rbind(c(edgeaffected[1],currnode), c(currnode,edgeaffected[2]))
		EDGE<-EDGE[-whereedgeaffected,] #remove old edge
		EDGE<-rbind(EDGE, newedg)
		currnode<-currnode+1
	}
	return(EDGE)
}

GetAllAncestors<-function(EDGE, tr) {
	##GET A LIST OF ALL ANCESTORS FOR EACH TIP
	EDGE2<-EDGE[match(1:Ntip(tr),EDGE[,2]),]
	x<-Ntip(tr)
	while(sum(is.na(EDGE2[,1]))<x) {
		EDGE2<-cbind(EDGE[match(EDGE2[,1],EDGE[,2]),1],EDGE2)
	}
	cbind(EDGE[match(EDGE[,1],EDGE[,2]),1],EDGE)

	AllAnc<-apply(EDGE2,1,function(x) rev(x[!is.na(x)]))  ##ALL ANCESTORS FOR EACH TIP
	return(AllAnc)
}

GetClosestExistingDesc<-function(EDGE, from) {
	subedge<-rbind(EDGE[EDGE[,1]<=(-from),],EDGE[EDGE[,1]>=from,])
	x<-nrow(subedge)
	col<-2
	while(sum(is.na(subedge[,col]))<x) {
		subedge<-cbind(subedge, subedge[match(subedge[,col], subedge[,1]),2])
		col<-col+1  
	}
	res<-t(apply(subedge,1,function(x) x[!is.na(x)][c(1,length(x[!is.na(x)]))]))
	return(res)
}

GetClosestExistingDesc2<-function(EDGE, from) {
	subedge<-rbind(EDGE[EDGE[,1]<=(-from),],EDGE[EDGE[,1]>=from,])
	x<-nrow(subedge)
	col<-2
	while(sum(is.na(subedge[,col]))<x) {
		subedge<-cbind(subedge, subedge[match(subedge[,col], subedge[,1]),2])
		col<-col+1  
	}
	res<-t(apply(subedge,1,function(x) x[!is.na(x)][c(1,length(x[!is.na(x)]))]))
	return(res)
}



GetEdgeFullMatrix<-function(AllAnc) {
	#For a given list of ghosts species, we must be able to provide the list of extinct nodes (without descendants)
	rows<-setdiff(unlist(AllAnc), 1:length(AllAnc))
	DECOMPOSED<-matrix(0, ncol=length(AllAnc), nrow=length(rows))
	rownames(DECOMPOSED)<-rows
	for (i in 1:length(AllAnc)) {
		DECOMPOSED[as.character(AllAnc[[i]][-1]),i]<-1
	}
	return(DECOMPOSED)
}

GetCorrespForAllNodes<-function(tr, ghosts) {
	if (length(ghosts)==0) { #there are no ghosts
		res<-cbind(1:(Ntip(tr)+Nnode(tr)), 1:(Ntip(tr)+Nnode(tr)))
		return(res)
	}
	desc<-Descendants(tr, type="tip")
	descdead<-lapply(desc, function(x,y) setdiff(x,y), y=ghosts)
	deadnodesTF<-unlist(lapply(descdead, function(x) length(x)>0))
	edge<-tr$edge
	edgesmall<-edge[edge[,2] %in% which(deadnodesTF),] 
	singletons<-as.numeric(names(which(table(edgesmall[,1])==1)))
	edgesmall2<-edgesmall[match(singletons, edgesmall[,1]),]
	n<-3
	nr<-nrow(edgesmall2)
	repeat {
		edgesmall2<-cbind(edgesmall2, edgesmall2[match(edgesmall2[,n-1], edgesmall2[,1]),2])
		if (sum(is.na(edgesmall2[,n]))==nr) break
		n<-n+1
	}
	res<-t(apply(edgesmall2, 1, function(x) x[!is.na(x)][c(1, length(x[!is.na(x)]))]))
	#add unaffected tips and nodes
	unaff<-setdiff(1:(Ntip(tr)+Nnode(tr)), res[,1])
	res<-rbind(res, cbind(unaff,unaff))
	return(res)
}

GetSmallEdgeMatrix<-function(tr, ghosts) { #from a tree and a list of absent species, returns the small edge matrix that removes all nodes that have become absent, but keeping the same node names	
	if (length(ghosts)==0) return(list(edge=tr$edge, M=1:(Ntip(tr)+Nnode(tr))))
	ntips<-Ntip(tr)
	nnodes<-Nnode(tr)
	M<-1:(ntips+nnodes) #will keep correspondences. ith element of M indicates the node corresponding to i after removing ghosts.
	edge<-tr$edge	
	edge<-edge[-match(ghosts, edge[,2]),]
	repeat {
		k0<-nrow(edge)
		#remove nodes without descendants
		nodescnodes<-which((edge[,2]>ntips)&!(is.element(edge[,2],edge[,1])))
		M[edge[nodescnodes,2]]<-0
		M[is.element(M, edge[nodescnodes,2])]<-0
#		print(paste(edge[nodescnodes,2],"del"))		
		if (length(nodescnodes)>0) edge<-edge[-nodescnodes,]
		#get list of singletons
		solo<-as.numeric(names(which(table(edge[,1])==1)>0))
		if (length(solo)>0) {
			n<-sort(solo)[1]
			M[edge[match(n, edge[,2]),2]]<-edge[match(n, edge[,1]),2]
			M[which(M==edge[match(n, edge[,2]),2])]<-edge[match(n, edge[,1]),2]
#			print(edge[match(n, edge[,2]),2])
#			print(edge[match(n, edge[,1]),2])
			edge[match(n, edge[,2]),2]<-edge[match(n, edge[,1]),2]
			edge<-edge[-match(n, edge[,1]),]
		}
		k1<-nrow(edge)
		if (k1==k0) break
	}
	return(list(edge=edge, M=M)) 
}

validateorinvalidatetransfers<-function(edgesmall, transfarr) {
	# print(edgesmall)
	# print(transfarr)
	parentsofsourceandrecipient<-edgesmall[match(transfarr[1:2], edgesmall[,2]),1]
	ret<-TRUE
	if (length(unique(parentsofsourceandrecipient))==1) ret<-FALSE #recipient and donor are sister lineages
	if (sum(apply(edgesmall,1,function(x,y) c(sum(abs(x-y))==0, sum(abs(x-rev(y)))==0),y=transfarr[1:2]))) ret<-FALSE
	# if (ret==FALSE) {
	# 	print (paste("a transfer was removed: ",transfarr[1]," --> ", transfarr[2], sep=""))
	# }
	ret
}


PredictTransfers<-function(tr, events, ghosts=c(3,4,5), plot=TRUE, removeUndetectable=TRUE) {
	EDGE<-GetExtendedEdgeMatrix(tr, events)
	newnodestimes<-cbind(-((Ntip(tr)+Nnode(tr)+1):(Ntip(tr)+Nnode(tr)+nrow(events))), events$TIME) #useful to get times at the end
	ALLORDEREDANCESTORS<-GetAllAncestors(EDGE, tr)
	FULLMATRIX<-GetEdgeFullMatrix(ALLORDEREDANCESTORS)
	smalledge<- GetSmallEdgeMatrix(tr, ghosts) #edge matrix of the tree after ghaving removed ghosts
	nonghost<-setdiff(1:Ntip(tr),ghosts)
	if (length(ghosts)==0) extinctnodes<-NULL #there are no extinct nodes
	else extinctnodes<-c(ghosts, as.numeric(names(which(apply(FULLMATRIX[,-ghosts],1,sum)==0))))

	if (plot) { #we draw the transfers
		drawtransfers(events, tr, extinctnodes) #draw tree and transfers
	}
	ntnn<-Ntip(tr)+Nnode(tr)+nrow(events)

	MAT<-matrix(0, nrow=ntnn, ncol=ntnn)


	alreadydeltwith<-Ntip(tr)+1 #liste des noeuds déjà vus et pas la peine de les revisiter. On met la racine au début
	#for (s in nonghost) { ###C'EST PARTI !!
#	TR<-NULL
	TRLIST<-list()
	THEPATH<-array(dim=ntnn)
	nb<-0
	extinctpath<-NULL
	for (s in nonghost) {
		# print(s)
		TR<-NULL
		# print(s)
		transfermode<-0
		extinctmode<-0
		repeat {
			# print(paste(" ",s,sep=""))
			if (s<0) { #transfer is received
				up<-(-s)
				if (transfermode==0) { #we are not already in a transfer mode where we search for the donor. 
					tr.recipient<-s
					# print(paste("    recipient: ",tr.recipient, sep=""))
					transfermode=1
				}
			}
			else {
				up<-EDGE[EDGE[,2]==s,1]
				# print(up)
			}
			if (!is.element(up, extinctnodes)) {
				if (transfermode==1) {
					tr.donor<-up
					# print(paste("    donor: ",tr.donor, sep=""))
					transfermode<-0
					TR<-rbind(TR, c(tr.donor, tr.recipient))
				}
				if (extinctmode==1) { #we were in an extinctmode
					# print("extinctpath ------------>")
					extinctpath<-c(extinctpath, up)
					# print("here is the extinct path before fixing it definitely")
					THEPATH[abs(extinctpath)]<-up
					# print(extinctpath)
					# print(extinctpath)
					# print("< ------------")
					extinctmode<-0
				}

			}
			else { #we are following an extinct path. We want to store all nodes met in this path up to the non-extinct node and s
				if (extinctmode==0) {
					# print("enterextinctmode")
					extinctpath<-up
					extinctmode<-1
				}
				else {
					# print("followextinctmode")
					extinctpath<-c(extinctpath, up)
					# print(extinctpath)
				}
			}
			if (is.element(up, alreadydeltwith)) {
				# print(up)
				# print("I saw it in the past, stop")
				if (transfermode==1) { #This is the missing case! we need to get the donor despite the fact that it has been seen 
					# print(paste("upANDstop", up, sep = "----"))
					tr.donor<-THEPATH[up]
				# 	#and it is not visited again
				 	TR<-rbind(TR, c(tr.donor, tr.recipient))
				 	transfermode<-0
				 	# print(paste("    donor: ",tr.donor, sep=""))				
				}
				##because we exit here in this case, we must store also the THEPATH at this moment. 
				if (extinctmode==1) {
					THEPATH[abs(extinctpath)]<-THEPATH[up]
					extinctmode<-0
				}
				# print("alreadydeltwith")
				# print(subnodes)
				break
			}
			alreadydeltwith<-c(alreadydeltwith, up)
			s<-up
		}
		if (!is.null(TR)) {
			nb<-nb+1
			TRLIST[[nb]]<-TR
		}
#		print(TR)
#		print(THEPATH)
	}
	if (length(TRLIST)==0) TRfFINAL<-NULL
	###ADD TIMES
	else {
		TRLIST<-lapply(TRLIST, function(x,nnt) cbind(x, nnt[match(x[,2], nnt[,1]),2]), nnt=newnodestimes)
		## At the end we must make a recoding of the nodes 
		## so that we return node names that do exist in the tree. 
		from<-Ntip(tr)+Nnode(tr)+1
		correspnodes<-GetClosestExistingDesc(EDGE, from=from) #correspondance between transfer nodes and nodes
		correspnodes.2<- GetCorrespForAllNodes(tr, ghosts) #correspondance between nodes and nodes that exist in small tree
		correspnodes<-cbind(correspnodes, correspnodes.2[match(correspnodes[,2], correspnodes.2[,1]),2])
		correspnodes<-t(apply(correspnodes, 1, function(x) x[!is.na(x)][c(1,length(x[!is.na(x)]))]))
		CORRESPNODES<-rbind(correspnodes, correspnodes.2)

		TRf<-lapply(TRLIST, function(x) t(apply(x,1,function(y,corr) c(corr[match(y[1:2],corr[,1]),2], y[3]),corr=CORRESPNODES)))

		simplifyit<-function(mat) {
			if (nrow(mat)>1) {
				whereiscross1<-mat[1:(nrow(mat)-1),1]-mat[2:(nrow(mat)),2]
				whereiscross2<-mat[1:(nrow(mat)-1),2]-mat[2:(nrow(mat)),1]
				rows2remove<-intersect(which(whereiscross1==0), which(whereiscross2==0))
				if (length(rows2remove)>0) mat<-mat[-rows2remove, ,drop=FALSE]
			}
			return(mat)
		}
		if (removeUndetectable) {
			#REMOVE CASES OF TRANSFERS A->B followed by B->A that are replaced by simply A->B
			TRf<-lapply(TRf, simplifyit) #toremove ? 
		}
		TRf<-do.call(rbind, TRf)
		if (removeUndetectable) {
			TRfFINAL<-TRf[apply(TRf,1, validateorinvalidatetransfers, edgesmall=smalledge$edge),,drop=FALSE]
		}
		else TRfFINAL<-TRf
		if (nrow(TRfFINAL)==0) TRfFINAL<-NULL
		else TRfFINAL<-TRfFINAL[order(TRfFINAL[,3]),, drop=FALSE]
	}
	return(TRfFINAL)
}


readandformatevents<-function(file) {
	# events<-read.table("../data-test-theo/sim/G/Gene_families/1_events.tsv", header=TRUE)
	events<-read.table(file, header=TRUE)
	events<-events[events$EVENT=="T",]
	events<-events[order(events$TIME),]
	# #transform events to add FROM an TO columns (useful for later)
	FROMTO<-do.call(rbind,lapply(strsplit(events$NODES, ";"), function(x) as.numeric(x[c(1,5)])))
	# events$FROM<-as.character(FROMTO[,1])
	# events$TO<-as.character(FROMTO[,2])
	tipandnodes<-c(tr$tip.label, tr$node.label)
	events$FROM<-match(as.character(FROMTO[,1]), tipandnodes)
	events$TO<-match(as.character(FROMTO[,2]), tipandnodes)
	return(events)
}

readandformatsmallgenetree<-function(file) {
	gtr<-read.tree(file)
	gtr$tip.label<-gsub("_.+","",gtr$tip.label)
	gtr$node.label<-gsub("_.+","",gtr$node.label)
	ghosts<-setdiff(tr$tip.label, gtr$tip.label)
	ghosts<-match(ghosts, tr$tip.label)
	return(list(gtr=gtr, ghosts=ghosts))
}

renamenodes<-function(tr, gtr, trmat, ghosts) {
	if (is.null(trmat)) return(NULL)
	else  {
		tipandnode<-c(tr$tip.label, tr$node.label)
		res<-data.frame(FROM=tipandnode[trmat[,1]], TO=tipandnode[trmat[,2]])
		# #for each node we get its descendant species
		# step1<-Descendants(tr, 1:(Ntip(tr)+Nnode(tr)), type="tips")
		# step2<-lapply(step1, function(x,y) setdiff(x,y), y=ghosts)
		# step3<-lapply(step2, function(x,y) y[x],y=tr$tip.label)
		# step4<-lapply(step3, function(x,y) getMRCA(y, x), y=gtr)
		# ntipgtr<-Ntip(gtr)
		# step5<-lapply(step4, function(x,y,w) y[x-w], y=gtr$node.label, w=ntipgtr)
		# step5[1:Ntip(tr)]<-step3[1:Ntip(tr)]
		# ##step5 contains the perfect mapping between big tree and small tree.
		# step6<-t(apply(trmat, 1, function(x,y) unlist(y[x]),y=step5)) 
		# res<-as.data.frame(step6)
		# colnames(res)<-c("FROM","TO")
		return(res)
	}
}

Newick <- function(edge, rmNodeNames=TRUE) {
	withbl<-ifelse (ncol(edge)>2, TRUE, FALSE) #check if branch lengtha re available
	build_newick_tree <- function(edge_matrix, current_node) {
	  # Get children of the current node
	  children <- edge_matrix[edge_matrix[,1] == current_node, ]
	  
	  # If there are no children, return the current node
	  if (nrow(children) == 0) {
	    return(current_node)
	  }
	  
	  child_strings <- character(0)
	  
	  # Loop through children
	  for (i in 1:nrow(children)) {
	    child_node <- children[i, 2]
	    child_length <- ifelse(withbl, children[i, 3],NA)
	    child_string <- build_newick_tree(edge_matrix, child_node)
	    
	    # Add branch length if available
	    if (!is.na(child_length)) {
	      child_string <- paste0(child_string, ":", child_length)
	    }
	    
	    child_strings <- c(child_strings, child_string)
	  }
	  
	  return(paste("(", paste(child_strings, collapse = ","), ")", as.character(current_node), sep = ""))
	}
	root<-unique(edge[!is.element(edge[,1], edge[,2]),1])

	newick_tree <- build_newick_tree(edge, root)
	newick_tree<-paste(newick_tree,";",sep="")
	#remove node names
	if (rmNodeNames) {
		if (withbl) newick_tree<-gsub("\\)\\d+:", "):",newick_tree)
		else newick_tree<-gsub("\\)\\d+", ")",newick_tree)
		newick_tree<-gsub("\\)\\d+;", ");",newick_tree)
	}
	newick_tree
}

#################################################
#################################################
#### NEW FUNCTIONS DEVELOPPED FOR MINI-PAPER ####
#################################################
#################################################

#temp 
# edgemat<-smalledge
# transf<-Induced
# edgsave<-edg
# trsave<-tr



SIMUTRANSF<-function(edgemat, transf) {
	if (is.null(transf)) return(edgemat) #nothing to do if no transfers
	transf<-as.data.frame(transf)
	colnames(transf)<-c("FROM", "TO", "TIME")
	## transf are tghe transfers, from most ancient toi most recent 
	#we do a loop for each transfer. 
	#Each time we do the same : destroy an edge (the one of the recipient, cut an edge in two 
	# (teh one of the donor) and branch the recipient into the new edge. 
	dotheoperation<-function(edg, tr) {
		donbranch<-match(tr$FROM, edg[,2])
		if (is.na(donbranch)) { #the donor is the root
			#we create a new edge to go to the root
			edgeroot<-c(max(edg)+1, tr$FROM)
			edg<-rbind(edg, edgeroot)
			donbranch<-match(tr$FROM, edg[,2])
		}
		recbranch<-match(tr$TO, edg[,2])
		newnod<-max(edg)+1
		newedges<-rbind(c(edg[donbranch,1], newnod), c(newnod, edg[donbranch,2]),c(newnod, edg[recbranch,2]))
		#we must also remove the singleton node that appeared: the dad of the recipient node.
		dadofrecipient<-edg[recbranch,1] #keep before line that goes just after !!
		#remove first set of edges
		edg<-edg[-(c(donbranch, recbranch)),]
		edg<-rbind(edg,newedges)
		edg1dad<-match(dadofrecipient,edg[,2])
		edg2dad<-match(dadofrecipient,edg[,1])
		# #return assignations that change: 
		# #when removing a singleton, we have to make sure that next transfer that involve it
		# #can still be done. To do so we must return the new name of this node
		# #and change the transfers accordingly after
		returnnewassignation<-edg[edg2dad,1:2]

		if (is.na(edg1dad)) {
			#we can remove the root!
			edg<-edg[-edg2dad,]
			#and there a no new edges to create
		}
		else {
			newedges2<-c(edg[edg1dad,1],edg[edg2dad,2])
			#remove second set of edges
			edg<-edg[-(c(edg1dad,edg2dad)),]			
			edg<-rbind(edg, newedges2)			
		}
		rownames(edg)<-NULL
		return(list(edg=edg,newassoc=returnnewassignation))
	}
	for (i in 1:nrow(transf)) {
		if (transf$FROM[i]!=transf$TO[i]) { #this may occur sometimes after some cut/paste operations
		computenewedgeandassign<-dotheoperation(edgemat, transf[i,])
		edgemat<-computenewedgeandassign$edg
		transf$FROM[which(transf$FROM==computenewedgeandassign$newassoc[1])]<-computenewedgeandassign$newassoc[2]
		transf$TO[which(transf$TO==computenewedgeandassign$newassoc[1])]<-computenewedgeandassign$newassoc[2]
		}
	}
	#we must remove a root node that may have emerged : 
	singletonroot<-setdiff(unique(edgemat[,1]), edgemat[duplicated(edgemat[,1]),1])
	if (length(singletonroot)>0) {
		edgemat<-edgemat[-match(singletonroot, edgemat[,1]),]

	}
	return(edgemat)
}

EDGE2NEWICK<-function(edg) {
	ntips<-sum(is.na(match(edg[,2], edg[,1])))
	nnodes<-ntips-1
	thetree<-list(edge=edg, Nnode=nnodes, tip.label=1:ntips)	
}

   




plotAllTransf<-function(tre, tps, EDGg, allhgts) {
		plot(tre, show.tip.label=F, edge.color = "grey", edge.width=1.3)
		xy<-as.data.frame(plotPhyloCoor(tre))
		onlytrue<-EDGg[EDGg[,3]==TRUE,]
		segments(xy$xx[onlytrue[,1]], xy$yy[onlytrue[,2]], xy$xx[onlytrue[,2]], xy$yy[onlytrue[,2]], lwd=1.5)
		segments(xy$xx[onlytrue[,1]], xy$yy[onlytrue[,1]], xy$xx[onlytrue[,1]], xy$yy[onlytrue[,2]], lwd=1.5)
		arrows(allhgts[,3], xy$yy[EDGg[allhgts[,1],2]], allhgts[,3], xy$yy[EDGg[allhgts[,2],2]], length=0.15, angle=20, col="red")
			##PLOT
}

KEEPTRANSFER<-function(tree, EdgeWithGhost, hgt) {
	# HERE WE CHECK IF A TRANSFER IS KEPT OR NOT. IT IS KEPT ONLY IF THE RECIPIENT
	## IS NOT A GHOST
	return(EdgeWithGhost[hgt[2],3])
}

PREPARETREE<-function(tre, tps) {
	#prepare tre for fast lookup of ghost branches
	EDGg<-as.data.frame(tre$edge)
	Desc<-Descendants(tre)
	Descg<-lapply(Desc, function(x, tip) length(intersect(x, tip))>0, tip=tps) #TRUE signifie que cette branche (celle qui mène à ce noeud) existe.
	EDGg<-cbind(EDGg, unlist(Descg[EDGg[,2]]))

	return(EDGg)
}

TRANSF<-function(tree, plot=FALSE) {
	xy<-as.data.frame(plotPhyloCoor(tree))
	EDG<-tree$edge
	EDG<-cbind(EDG, xy$xx[EDG[,1]], xy$xx[EDG[,2]])
	# TRANSFER SOURCE
	cumsumedgelength<-cumsum(tree$edge.length)
	#source du transfert
	sample<-runif(1,0,max(cumsumedgelength))
	if (sample<cumsumedgelength[1]) {
		sourcebranch<-1
		oudanslabranche<-sample
		nodeinvolved<-tree$edge[sourcebranch,1]
	}
	else {
		sourcebranch<-max(which(sample>cumsumedgelength))+1
		oudanslabranche<-sample-cumsumedgelength[sourcebranch-1] #disatnce au début de la branche.
		nodeinvolved<-tree$edge[sourcebranch,1]
	}
	realtimeoftheevent<-xy$xx[nodeinvolved]+oudanslabranche

	# TRANSFER RECIPIENT
	#ensuite pour trouver le receveur, on pioche dans les branches une qui existe au moment du transfert mais qui ne soit pas la source.  
	possiblerecipient<-which((EDG[,3]<=realtimeoftheevent&EDG[,4]>=realtimeoftheevent))
	possiblerecipient<-possiblerecipient[possiblerecipient!=sourcebranch]
	recipientbranch<-ifelse(length(possiblerecipient)==1, possiblerecipient, sample(possiblerecipient,1))
	return(c(sourcebranch, recipientbranch, realtimeoftheevent))
}



COMPUTERF<-function(trs) {
	t1<-trs[[1]]
	t2<-trs[[2]]	
	#we must add a new tip in the tree (o), at the root, so that 
	#we compute the RF of ROOTED tree wich makes sense here 
	#because the tree compared are dated/rooted initially.	
	newtip<-rtree(1, tip.label="o")
	T1<-bind.tree(t1,newtip)
	T2<-bind.tree(t2,newtip)
	return(dist.topo(T1,T2))
}






