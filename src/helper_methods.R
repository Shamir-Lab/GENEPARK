###### load packages
library(sva)
library(ROCR)
library(CMA)
library(e1071)
library(preprocessCore)

applyQuantileNormalization<-function(gdsExpmat){
	if (max(gdsExpmat)>100 && min(gdsExpmat)>0){
		raw_x = log(normalize.quantiles.robust(gdsExpmat))
	}
	else{
		raw_x = normalize.quantiles.robust(gdsExpmat)
	}
	colnames(raw_x) <- colnames(gdsExpmat)
	rownames(raw_x)<-rownames(gdsExpmat)
	return(raw_x)
}


intensityFilter<-function(x,value,p){
	numBelow = length(which(x<value))
	currP = numBelow/length(x)
	return (currP>=p)
}


getROCScore <- function(predictions,labels, lb = NULL,useROCR=T){
	if (useROCR){
		rocr_obj = ROCR::prediction(predictions,labels, label.ordering = lb)
		roc_score = ROCR::performance(rocr_obj, "tpr","fpr")
		roc_score = ROCR::performance(rocr_obj, "auc")
		return (slot(roc_score,"y.values")[[1]])
	}
	return (pROC::auc(labels,predictions))
}

calcAupr <- function(pred, gs) {
	ord.idx <- order(abs(pred), decreasing = T)
	prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
	rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
	fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate
	prec <- c(prec[1], prec)
	rec <- c(0, rec)
	fpr <- c(0, fpr)
	aupr <- areaUnderCurve(rec, prec)
	return (aupr)
	#auroc <- areaUnderCurve(fpr, rec)
	#return(auroc)
}

areaUnderCurve <- function(x, y) {
	dx <- diff(x)
	my <- y[1:(length(y) - 1)] + diff(y) / 2
	return(sum(dx * my))
}



# predictions here is a matrix
getMulticlassROCScore <- function(predictions,labels,all_rocs=T){
		classes = names(summary(as.factor(labels)))
		roc = 0
		rocs = c()
		sizes = c()
		num_samples = length(labels)
		for (cl in classes){
			inds = which(labels == cl)
			if (length(inds)==0){next}
			curr_size = length(inds)
			comp_name = paste("not_",as.character(cl),sep="")
			lb = c(comp_name,as.character(cl))
			curr_vector = rep(comp_name,num_samples)
			curr_vector[inds] = as.character(cl)
			currClassIndex = which(colnames(predictions)==cl)
			curr_roc = getROCScore(predictions[,currClassIndex],curr_vector,lb)
			roc = roc + (curr_roc*curr_size)/num_samples
			rocs = c(rocs,curr_roc)
			sizes = c(sizes,curr_size)
		}
		if (all_rocs){
			results = hash()
			results[["rocs"]] = rocs
			results[["sizes"]]= sizes
			results[["Weighted ROC"]]=roc
			return (results)
		}
		else{
			return(roc)
		}
		
}

# predictions here is a matrix
getMulticlassAccuracy <- function(predictions,labels){
		correct=0
		for (i in 1:length(labels)){
			curr_p = predictions[i,]
			max_p = max(curr_p)
			max_class = colnames(predictions)[which(curr_p==max_p)[1]]
			if (max_class == labels[i]){
				correct = correct+1
			}
		}
		return (correct/length(labels))
}

# data: first row contains the sample names
# path containes name to label mapping
getLablesFromFile<- function(data,path,hasHeader=F,colIndex=1){
	curr_cels = colnames(data)
	cel2label = read.table(path,row.names=1,header=hasHeader,sep="\t")

	# create the mapping of the cel files in the raw data matrix
	curr_cel2labels = c()
	for (i in 1:length(curr_cels)){
		currCel = as.character(curr_cels[i])
		curr_cel2labels[currCel] = as.character(cel2label[currCel,colIndex])
	}
	return (curr_cel2labels)
}

getFSVAObject<- function(training,test,known_variables){
	trainMod = model.matrix(~known_variables$disease,data=known_variables)
	trainMod0 = model.matrix(~1,data=known_variables)
	trainSv = sva(as.matrix(training),trainMod,trainMod0)
	fsvaobj = fsva(dbdat=as.matrix(training),mod=trainMod,sv=trainSv,newdat=as.matrix(test),method="exact")
	return (fsvaobj)
}

getFSVAClassificationAccuracyForEachBatch<-function(training,test,y,test_y,test_batches,known_variables,class1="1",class2="-1",...){

	# run fsva once
	fsvaobj = getFSVAObject (training,test,known_variables)
	tr = fsvaobj$db
	otherSamplesData_after_fsva = fsvaobj$new


	te_batch_set = names(summary(as.factor(test_batches)))
	batch2acc = hash()
	for (b in te_batch_set){
		curr_inds = which(test_batches==b)
		if (length(curr_inds)<2){next}
		print (length(curr_inds))
		curr_test = otherSamplesData_after_fsva[,curr_inds]
		curr_test_y = test_y[curr_inds]
		curr_results = getClassificationPerformance(tr,curr_test,y,curr_test_y,class1,class2,getROC=F,...)
		batch2acc[[b]]=curr_results[["accs"]]
	}
	return (batch2acc)
}

getAccuracy <- function (predictions,y,class1="1",class2="-1"){
        # get accuracy
        good_preds=0
        for (i in 1:length(y)){
                if (y[i]==class2 && predictions[i]<0.5){
                        good_preds = good_preds+1
                }
                if (y[i]==class1 && predictions[i]>0.5){
                        good_preds = good_preds+1
                }
        }
        accuracy = good_preds/length(y)
        return (accuracy)
}


hashToMatrix<-function(h){
	ks = sort(keys(h))
	v1 = h[[ks[1]]]
	m = matrix(nrow=length(ks),ncol=length(v1),0)
	rownames(m)<-ks
	for (k in ks){
		m[k,]=h[[k]]
	}
	return (m)
}

getFeatureLabelAssociation <- function(feature,labels, permute=T,brks = 2){
	h = hist(feature,nclass=brks, plot=F)
	newf = cut(feature,breaks = h$breaks,include.lowest=T)
	return (chisq.test(newf,as.factor(labels),simulate.p.value=permute,correct=T)$p.value)
}


# x = The data matrix. Columns correspond to samples.
#     Colnames = the names of the samples
# y = the labels
# batch_info = assignment of samples to batches
# useFSVA
# classification_function
# useFS = logical, use feature selection?
# fsmethod = relevant only is useFS=T
# sigSizes = the number of features to select if useFS=T
# prediction_args = additional parameters for the predict function
# ... additional parameters to the classification function or the predict function
#	empty for svm
#	c(type="prob") for randomForest
getLeaveBatchOutClassificationPerformance<-function(
	x,y,batch_info,useFSVA=FALSE,classification_function=svm,useFS=TRUE, giveBatchesToClassifier = F,
	fsmethod="limma",sigSizes = seq(from=10,to=200,10),
	prediction_args = c(),...){
	
	predictions = hash()
	batches = unique(batch_info)
	cels_test_order = c()

	known_variables = data.frame(cbind(batch_info,y),row.names=colnames(x))
	colnames(known_variables)<-c("batch","disease")

	# Leave batch out analysis
	for (b in batches){
		curr_inds = which(batch_info == b)
		tr_data = x[,-curr_inds]
		tr_y = y[-curr_inds]
		tr_batches = as.factor(as.character(batch_info[-curr_inds]))
		names(tr_batches)<-names(batch_info[-curr_inds])

		# sanity check		
		numClassesInTraining = length(which(table(as.factor(tr_y)) > 0))
		if (numClassesInTraining < length(table(y))){
			print ("Training set does not contain all classes")
			return (NULL)
		}
		else{
			cels_test_order = c(cels_test_order,curr_inds)
		}

		te_data = x[,curr_inds]
		if (length(curr_inds)==1){te_data = matrix(te_data,ncol=1)}

		trainVars = known_variables[-curr_inds,]
		if (useFSVA){
			trainMod = model.matrix(~trainVars$disease,data=trainVars)
			trainMod0 = model.matrix(~1,data=trainVars)
			trainSv = sva(as.matrix(tr_data),trainMod,trainMod0)

			if (trainSv$n.sv == 0){
				print (" no significant surrogate variables")
				print (" rerun this method with useFSVA=F")
				return (NULL)
			}
			
			fsvaobj = fsva(dbdat=as.matrix(tr_data),mod=trainMod,sv=trainSv,newdat=te_data,method="exact")
			print ("fsva finished")
			train_data = t(fsvaobj$db)
			test_data = t(fsvaobj$new)		
		}
		else{
			train_data= t(tr_data)
			test_data = t(te_data)
		}

		# No Feature selection
		if (giveBatchesToClassifier){
			nofs_model = classification_function(train_data,tr_y,batches=tr_batches,...)
		}
		else{
			nofs_model = classification_function(train_data,tr_y,...)
		}
		
		if (length(prediction_args)>0){
			nofs_preds = predict(nofs_model,newdata=test_data,probability=T,prediction_args)
		}
		else{
			nofs_preds = predict(nofs_model,newdata=test_data,probability=T)
		}
		
		preds = getPredProbabilities(nofs_preds)
		key = "NO_FS"
		if (has.key(key,predictions)){
			predictions[[key]] = rbind(predictions[[key]],preds)
		}
		else{
			predictions[[key]] = preds
		}

		if (!useFS){next}
		
		selection1 = GeneSelection(train_data,tr_y,method=fsmethod)
		for (size in sigSizes){
			selected_inds = toplist(selection1,k=size,show=F)$index
			trainX_reduced = train_data[,selected_inds]
			textX_reduced = test_data[,selected_inds]

			if (giveBatchesToClassifier){
				model = classification_function(trainX_reduced,tr_y,batches=tr_batches,...)
			}
			else{
				model = classification_function(trainX_reduced,tr_y,...)
			}

			if (length(curr_inds)==1){
				textX_reduced = matrix(textX_reduced,nrow=1)
			}

			if (length(prediction_args)>0){
				pred = predict(model,newdata=textX_reduced,probability=T,prediction_args)
			}
			else{
				pred = predict(model,newdata=textX_reduced,probability=T)
			}
			
			curr_preds_mat = getPredProbabilities(pred)
		
			key = as.character(size)
			if (has.key(key,predictions)){
				predictions[[key]] = rbind(predictions[[key]],curr_preds_mat)
			}
			else{
				predictions[[key]] = curr_preds_mat
			}
		}
	}

	sizes = sort(hash::keys(predictions))
	labels = y[cels_test_order]
	results_table= c()
	for (size in sizes){
		curr_preds = predictions[[as.character(size)]]
		accuracy = getMulticlassAccuracy(curr_preds,labels)
		rocs = getMulticlassROCScore(curr_preds,labels)
		curr_roc_scores = c(rocs[["rocs"]],rocs[["Weighted ROC"]])
		results_table = rbind(results_table,c(size,curr_roc_scores,accuracy))
	}

	class_size = rocs[["sizes"]]
	colnames(results_table) = c("Size",class_size,"Weighted ROC","Accuracy")

	return(list(results_table=results_table,predictions=predictions))
}

getPredProbabilities <- function(pred_obj){
	probs = attr(pred_obj,"prob")
	if (!is.null(probs)){
		return(probs)
	}

	probs = attr(pred_obj,"probabilities")
	if (!is.null(probs)){
		return(probs)
	}

	if (is.element("matrix",class(pred_obj))){
		return (pred_obj)
	}

	probs = pred_obj[[2]]
	if (!is.null(probs)){
		return(probs)
	}

	return (NULL)

}

removePrefix<-function(x,pre){
        return (strsplit(x,split=pre)[[1]][2])
}

#################################### Enrichment analysis functions ##############################
getGeneSets<-function(kegg){
	pathways = names(kegg)
	p2g = list()
	for (p in pathways){
		curr_genes = sapply(graphite::nodes(kegg[[p]]),removePrefix,pre="EntrezGene:")
		p2g[[p]]=curr_genes
	}
	return(p2g)
}

getHGPvalue<-function(hitlist,term,bg){
	inters = length(intersect(hitlist,term))
	k = length(term)
	n = length(bg)-length(hitlist)
	m = length(hitlist)
	return (phyper(inters-1,m,n,k,lower.tail=FALSE))
}

getHGEnrichmentAnalysis<-function(userlists,terms,bg){
	results = c()
	for (u in names(userlists)){
		hitlist = intersect(userlists[[u]],bg)
		for (t in names(terms)){
			term = intersect(terms[[t]],bg)
			name = paste(u,t,sep=";")
			results[name] = getHGPvalue(hitlist,term,bg)
		}
	}
	results = p.adjust(results,method="fdr")
	return (results)
}





