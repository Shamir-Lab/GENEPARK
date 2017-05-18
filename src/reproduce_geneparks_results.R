# This is a script for reproducing the results of the GenePark project.
# Below we set the required libraries and input files for the analysis.
# A file with auxiliary functions: 'helper_methods.R'
# The functions below are the main machine learning and batch correction functions.

# The required libraries for running the code
library(sva);library(ROCR);library(CMA);library(e1071);library(hash);library(preprocessCore);library('org.Hs.eg.db')
# Set your working directory to the directory with "helper_methods.R" and all data from 
# http://acgt.cs.tau.ac.il/genepark_data/
# Do not forget to unzip the files
setwd('C:/Users/davidama/Desktop/genepark/analysis/classification_r/reproduce_results')
source("helper_methods.R")

##################################################################################################################
# Text files with the expression profiles
# IMPORTANT(!!!) - if these files were compressed, unzip them first.
expression_files = list()
expression_files[["training_set"]] = "training_set.txt"
expression_files[["training_and_validation_sets"]] = "training_and_validation_sets.txt"
expression_files[["validation_set"]] = "validation_set.txt"
expression_files[["test_set"]] = "test_set.txt"

# Additional text files with labels and batches
cel2date_file = "cel2date.txt"
cel2labels_file = "cel2full_labels.txt"

# Our custom file that mapps probes to Entrez gene
probe2entrez_EXPANDER_file = "HG-U133Plus2_Affy2Entrez.txt"

##############################################################################################
##############################################################################################
# A function for preprocessing of the training set
# The function applies quantile normalization and selects a set of probes based on an intensity filter.
# It also loads the labels and the batches and:
# Exclude samples of undesired classes 
# Exclude samples from small batches.
# minNumCels: remove samples from batches whose size is less than this number 
preprocess_data_set<-function(new_data_raw,labels_path = "../cel2full_labels.txt",l_to_keep = c("IPD","CONTROL"),runQuantileNorm=F,
	batch_path = "../cel2date.txt", minNumCels = 15){
	new_data_raw = data.matrix(new_data_raw)
	if(runQuantileNorm){new_data_raw = applyQuantileNormalization(new_data_raw)}
	toRemove = apply(new_data_raw,1,intensityFilter,value=6,p=0.8)
	new_data_filter1= new_data_raw[-which(toRemove),]
	colnames(new_data_filter1)<-colnames(new_data_raw)
	rownames(new_data_filter1)<-rownames(new_data_raw)[-which(toRemove)]
	# get the labels
	new_data_full_labels = as.factor(getLablesFromFile(new_data_filter1,labels_path))
	# keep the desired labels
	to_keep = which(sapply(new_data_full_labels,is.element,set=l_to_keep))
	data_labels = as.factor(as.character(new_data_full_labels[to_keep]))
	new_data_filter1 = new_data_filter1[,to_keep]
	data_batches = as.factor(getLablesFromFile(new_data_filter1,batch_path))
	# get the batches
	if(minNumCels > 1){
		b_count = table(data_batches)
		small_batches = names(b_count[which(b_count<=minNumCels)])
		to_remove = which(sapply(data_batches,is.element,set=small_batches))
		minus_list = -to_remove
		# filter2: remove batches with a small number of samples
		new_data_filter2 = new_data_filter1[,minus_list]
		curr_cel2date =  as.factor(as.character(data_batches[minus_list]))
		curr_labels = data_labels[minus_list]
		return(list(preprocessed_data = new_data_filter2,cel2batch=curr_cel2date,cel2label=curr_labels))
	}
	return(list(preprocessed_data = new_data_filter1,cel2batch=data_batches,cel2label=data_labels))
}

##############################################################################################
# Preprocessing of the training set
tr_data_raw = read.table(expression_files[["training_set"]],header=T,sep="\t",row.names=1,as.is=T,check.names=F)
tr_data_obj = preprocess_data_set(tr_data_raw,labels_path=cel2labels_file,batch_path = cel2date_file,minNumCels = 10)

# Preprocessing of the union of the training and the validation sets
trv_data_raw = read.table(expression_files[["training_and_validation_sets"]],header=T,sep="\t",row.names=1,as.is=T,check.names=F)
trv_data_obj = preprocess_data_set(trv_data_raw,labels_path=cel2labels_file,batch_path = cel2date_file,runQuantileNorm=T)

##### Data deposited in GEO
genepark_data_rma = cbind(trv_data_raw,read.table(expression_files[["test_set"]],header=T,sep="\t",row.names=1,as.is=T,check.names=F))
pheno_data = cbind(colnames(genepark_data_rma),
	getLablesFromFile(genepark_data_rma,cel2labels_file),
	getLablesFromFile(genepark_data_rma,cel2date_file)
)
v = rep("TEST",ncol(genepark_data_rma))
v[1:ncol(trv_data_raw)] = "VALIDATION"
names(v) = colnames(genepark_data_rma)
v[is.element(names(v),set=colnames(tr_data_raw))] = "TRAINING"
table(v)
pheno_data = cbind(pheno_data,v)
colnames(pheno_data) = c("Title","Label","Batch","Learning set")
write.table(genepark_data_rma,file="GenePark_data_for_GEO_RMA_data.txt",sep="\t",quote=F)
write.table(pheno_data,file="GenePark_data_for_GEO_pheno_data.txt",sep="\t",quote=F)
######

# Leave-batch-out analysis
trv_lbo_results_limma = getLeaveBatchOutClassificationPerformance(x=trv_data_obj[[1]],y=trv_data_obj$cel2label,
				batch_info=trv_data_obj$cel2batch,fsmethod='limma'
				,useFSVA=T,prediction_args = c(probability=T),
				probability=1,kernel="linear",scale=F)
plot(as.numeric(trv_lbo_results_limma[[1]][,4]),type='b',ylim = c(0.5,1))

# A function for the analysis of a new set of samples
# The function uses the training data and:
# Trim the test set to the selected probes
# Partition the test set to samples from the same batches as the training ("large")
# and samples from new batches ("small")
# It then runs the fSVA-Limma-SVM analysis to get predictions
# on these test sets.
get_analysis_on_a_new_sample_set<-function(tr_set_data,tr_set_batch,
		tr_set_labels,test_set_raw,cel2labels_file,cel2date_file,
		sigSizes = c(50,100),runQuantileNorm=F){
	known_variables = data.frame(cbind(tr_set_batch,tr_set_labels),row.names=colnames(tr_set_data))
	colnames(known_variables)<-c("batch","disease")
	tr.y = as.factor(tr_set_labels)
	test_set_raw = data.matrix(test_set_raw)
	if(runQuantileNorm){test_set_raw =  applyQuantileNormalization(test_set_raw)}
	# use the probes in the preprocessed data
	uns_filter_probe_names = rownames(tr_set_data)
	validation_uns_filter_data = test_set_raw[uns_filter_probe_names,]
	colnames(validation_uns_filter_data)<-colnames(test_set_raw)
	# get current labels
	validation_cel2label = getLablesFromFile(test_set_raw,cel2labels_file)
	ipd_c_inds = sort(c(which(validation_cel2label=="IPD"),which(validation_cel2label=="CONTROL")))
	# get the relevant tested labels data
	v_data = validation_uns_filter_data[,ipd_c_inds]
	colnames(v_data)<-colnames(test_set_raw)[ipd_c_inds]
	v_y = as.factor(getLablesFromFile(v_data,cel2labels_file,hasHeader=F))
	#### validation set analysis: get batches that appear in the training
	in_tr_batches =  unique(tr_set_batch)
	te_batches =  getLablesFromFile(v_data,cel2date_file,hasHeader=F)
	selected_inds = c()
	for (batch in in_tr_batches){
		curr_inds = which(te_batches==batch)
		selected_inds = c(selected_inds,curr_inds)
	}
	selected_inds = sort(selected_inds)
	#### Samples in the test from the training's small batches (i.e., not used)
	v_small_batches_data = v_data[,-selected_inds]
	v_small_batches.l = getLablesFromFile(v_small_batches_data,cel2labels_file,hasHeader=F)
	v_small_batches.y = as.factor(v_small_batches.l)
	v_small_batches.b = getLablesFromFile(v_small_batches_data,cel2date_file,hasHeader=T)
	#### Samples in the test from the training's large batches
	v_large_batches_data = v_data[,selected_inds]
	v_large_batches.l = getLablesFromFile(v_large_batches_data,cel2labels_file,hasHeader=F)
	v_large_batches.y = as.factor(v_large_batches.l)
	v_large_batches.b = getLablesFromFile(v_large_batches_data,cel2date_file,hasHeader=T)
	# Run classification with fSVA
	v_small_performance = getFSVAClassificationFlowPerformanceOnTestSet(tr_set_data,v_small_batches_data,tr.y,v_small_batches.y,known_variables,class1,class2,sigSizes=sigSizes)
	v_large_performance = getFSVAClassificationFlowPerformanceOnTestSet(tr_set_data,v_large_batches_data,tr.y,v_large_batches.y,known_variables,class1,class2,sigSizes=sigSizes)
	return(list(v_small_batches_data=v_small_batches_data,v_small_batches.l=v_small_batches.l,
			v_small_batches.y =v_small_batches.y,v_small_batches.b=v_small_batches.b,
			v_large_batches_data=v_large_batches_data,v_large_batches.l=v_large_batches.l,
			v_large_batches.y =v_large_batches.y,v_large_batches.b=v_large_batches.b,
			v_small_performance=v_small_performance,v_large_performance=v_large_performance,
			known_variables=known_variables))
}

# This function applies fSVA model on a new test set
getFSVAClassificationFlowPerformanceOnTestSet <- function(training,test,y,test_y,known_variables,class1="1",class2="-1",...){
        fsvaobj = getFSVAObject (training,test,known_variables)
        training = fsvaobj$db; otherSamplesData_after_fsva = fsvaobj$new
        return (getClassificationPerformance(training,otherSamplesData_after_fsva,y,test_y,class1,class2,...))
}

# A method that learns different SVM models for different signature sizes
# and returns the performance scores.
getClassificationPerformance <-function(training,otherSamplesData_after_fsva,y,test_y,class1="1",class2="-1", getROC=T,fsmethod='limma',sigSizes = seq(from=10,to=200,10)){
	selection = GeneSelection(t(training),y,method=fsmethod)
	top_list = toplist(selection,k=1000,show=T)
	print (rownames(training)[top_list$index[1:100]])
	te_rocs=c()
	te_accs=c()
	preds = c()
	for (size in sigSizes){
		selectedInds = top_list$index[1:size]
		curr_tr = t(training[selectedInds,])
		curr_te = t(otherSamplesData_after_fsva[selectedInds,])
		model = svm(curr_tr,y, kernel="linear",scale=F,probability=1)
		pred = predict(model,newdata=curr_te,probability=T)
		predClass1Index = which(colnames(attr(pred,"probabilities"))==class1)
		curr_preds = attr(pred,"prob")[,predClass1Index]
		preds = rbind(preds,curr_preds)
		curr_roc=0
		if (getROC){curr_roc = getROCScore(curr_preds,test_y)}
		curr_acc = getAccuracy(curr_preds,test_y,class1,class2)
		te_rocs = c(te_rocs,curr_roc)
		te_accs = c(te_accs,curr_acc)
	}
	rownames(preds) = sigSizes
	results = hash()
	results[["rocs"]]=te_rocs
	results[["accs"]]=te_accs
	results[["size"]]=sigSizes
	results[['preds']] = preds
	return (results)
}
# Auxiliary functions for dealing with ROC curves.
plot_roc<-function(p,y,set_mar=T,...){
	if (set_mar){par(mar=rep(5,4))}
	pred_obj <- ROCR::prediction(p,y)
	perf_roc <- ROCR::performance(pred_obj, measure = "tpr", x.measure = "fpr") 
	ROCR::plot(perf_roc,...)
	roc_score <- ROCR::performance(pred_obj, measure = "auc")
	return(slot(roc_score,"y.values")[[1]])
}

get_auc_pval<-function(p,y){
	x1 = p[y==y[1]]
	x2 = p[y!=y[1]]
	return (wilcox.test(x1,x2)$p.value)
}
##############################################################################################

# Global variables for our IPD-CONTROL analysis
class1 = "IPD";class2 = "CONTROL"
###############################################
# Training + validation vs. the test set
tr_set_data = trv_data_obj[[1]]
tr_set_batch = trv_data_obj[[2]]
tr_set_labels = trv_data_obj[[3]]
test_set_raw = read.table(expression_files[["test_set"]],row.names=1,header=T,check.names=F)
trv_vs_test_analysis = get_analysis_on_a_new_sample_set(tr_set_data,tr_set_batch,tr_set_labels,test_set_raw,cel2labels_file,cel2date_file,runQuantileNorm=T)
# Training vs. validation
tr_set_data2 = tr_data_obj[[1]]
tr_set_batch2 = tr_data_obj[[2]]
tr_set_labels2 = tr_data_obj[[3]]
vali_set_raw = read.table(expression_files[["validation_set"]],row.names=1,header=T,check.names=F)
tr_vs_vali_analysis = get_analysis_on_a_new_sample_set(tr_set_data2,tr_set_batch2,tr_set_labels2,vali_set_raw,cel2labels_file,cel2date_file,runQuantileNorm=T)

te_small_batches_data = trv_vs_test_analysis$v_small_batches_data
te_small_batches.y = trv_vs_test_analysis$v_small_batches.y
te_large_batches_data = trv_vs_test_analysis$v_large_batches_data
te_large_batches.y = trv_vs_test_analysis$v_large_batches.y
sig_size = '100'
te_p_s = c(trv_vs_test_analysis$v_small_performance[['preds']][sig_size,])
te_y_s = c(trv_vs_test_analysis$v_small_batches.y)
te_p_l = c(trv_vs_test_analysis$v_large_performance[['preds']][sig_size,])
te_y_l = c(trv_vs_test_analysis$v_large_batches.y)
te_p = c(te_p_s,te_p_l)
te_y = c(te_y_s,te_y_l)
v_p_s = c(tr_vs_vali_analysis$v_small_performance[['preds']][sig_size,])
v_y_s = c(tr_vs_vali_analysis$v_small_batches.y)
v_p_l = c(tr_vs_vali_analysis$v_large_performance[['preds']][sig_size,])
v_y_l = c(tr_vs_vali_analysis$v_large_batches.y)
v_p = c(v_p_s,v_p_l)
v_y = c(v_y_s,v_y_l)

# Figure 1C
v_roc_score = plot_roc(c(v_p_l,v_p_s),c(v_y_l,v_y_s),col='black',lwd=5,cex.axis=3,cex.lab=2,font.lab=6)
get_auc_pval(c(v_p_l,v_p_s),c(v_y_l,v_y_s))
te_roc_score = plot_roc(te_p,te_y,col='green',lwd=5,cex.axis=3,cex.lab=2,font.lab=6,add=T)
get_auc_pval(te_p,te_y)
abline(0,1,lwd=4,lty=2)
legend(x = 'bottomright',legend = c("Validation set","Test set"),col=c('black','green'),lty=1,lwd=5,cex=2)

# Figure 1D
# Test: new vs old batches
te_old_roc_score = plot_roc(te_p_l,te_y_l,col='red',lwd=5,cex.axis=3,cex.lab=2,font.lab=6)
te_new_roc_score = plot_roc(te_p_s,te_y_s,col='blue',lwd=5,cex.axis=3,cex.lab=2,font.lab=6,add=T)
abline(0,1,lwd=4,lty=2)
legend(x = 'bottomright',legend = c("Test set: old batches","Test set: new batches"),col=c('red','blue'),lty=1,lwd=5,cex=2)

####################################### Analyze samples of other diseases and the selected probe signature ###############################
# Prepare the new data to be tested - the samples should be
# both the excluded IPD-CONTROL test set and the set of samples
# with other NDDs that were not used for training the models.
d1 = trv_data_raw;d2 = test_set_raw
d1_full_labels = getLablesFromFile(d1,cel2labels_file)
# remove IPD and controls from d1
l_to_remove = c("IPD","CONTROL")
to_remove = which(sapply(d1_full_labels,is.element,set=l_to_remove))
d1 = d1[,-to_remove]
new_data_raw = as.matrix(cbind(d1,d2))
new_data_raw = applyQuantileNormalization(new_data_raw)

# Remove low variation probes according to training
probe_list = rownames(trv_data_obj[[1]])
new_data_filter1 = new_data_raw[probe_list,]

# Get the labels of the new samples
new_data_full_labels = getLablesFromFile(new_data_raw,cel2labels_file)

# Do fsva on the traning and apply the model on the test
training = trv_data_obj[[1]]
known_variables = trv_vs_test_analysis$known_variables
test = new_data_filter1
trainMod = model.matrix(~known_variables$disease,data=known_variables)
trainMod0 = model.matrix(~1,data=known_variables)
trainSv = sva(as.matrix(training),trainMod,trainMod0)
fsvaobj = fsva(dbdat=as.matrix(training),mod=trainMod,sv=trainSv,newdat=as.matrix(test),method="exact")
tr_set_data_fsva = fsvaobj$db
new_data_fsva = fsvaobj$new

########################## Feature selection: get the signature #################################
selection = GeneSelection(t(tr_set_data_fsva),trv_data_obj[[3]],method='limma')
top_list = toplist(selection,k=100,show=F)
top_list[,1] = rownames(tr_set_data_fsva)[top_list$index]

# Mapping probes to symbols
probe2entrez = read.table(probe2entrez_EXPANDER_file,row.names=1)
probe2entrez = probe2entrez[top_list[,1],1]
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
entrez2symbol <- unlist(as.list(x[mapped_genes]))
sig_symb = entrez2symbol[as.character(probe2entrez)]
top_list = cbind(top_list[,1],probe2entrez,sig_symb,top_list[,-1])
# Write the results into a file to get Supplementary Table 3
write.table(top_list,file="Supplementary_Table_3.txt",sep="\t",quote=F)

######################### Get the predictions on the test set, separated by classes ##############
class1 = "IPD";size = 100
selectedInds = top_list[,1]
curr_tr = 	t(tr_set_data_fsva[selectedInds,])
curr_te = t(new_data_fsva[selectedInds,])
model = svm(curr_tr,trv_data_obj[[3]], kernel="linear",scale=F,probability=1)
pred = predict(model,newdata=curr_te,probability=T)
predClass1Index = which(colnames(attr(pred,"probabilities"))==class1)
predictions = (attr(pred,"probabilities"))[,predClass1Index]

# L2872 has MSA
new_data_full_labels[which(grepl(names(new_data_full_labels),pattern = "L2872"))] = 'MSA'
new_sample_labels = new_data_full_labels
table(new_sample_labels)
# to make the plot easier to follow we unite some classes here
counts = table(new_sample_labels)
small_classes = names(counts[which(counts<4)])
is_in_small_classes = sapply(new_sample_labels,is.element,set=small_classes)
small_classes_inds = which(is_in_small_classes)
new_sample_labels[small_classes_inds] = "Other NDDs"

classes = unique(new_sample_labels)
# put Control and IPD at the begining
classes = unique(c("CONTROL","IPD",classes))
classes_preds = list()
must_have_classes = c("CONTROL","IPD")
for (cl in classes){
	inds = which (new_sample_labels == as.character(cl))
	if (length(inds)<1 && !is.element(cl,set=must_have_classes)){next}; 
	name = paste(as.character(cl),length(inds),sep=":")
	classes_preds[[name]]=predictions[inds]
}

############### Supplementary Figure 3 ##################
classes_preds[['HD:27']] =  c(classes_preds[["HD_HD_BATCH:8"]],classes_preds[["HD:19"]])
classes_preds[['Atypical PD:19']] = c(classes_preds[["MSA:9"]],classes_preds[["PSP:8"]],classes_preds[["CBD:2"]])
classes_preds[['NDD:48']] = c(classes_preds[['Atypical PD:19']],classes_preds[['HD:27']],classes_preds[["PD_DEMENTIA:2"]])
cols = c("blue","red","white")
boxplot(classes_preds[c("CONTROL:40","IPD:30",'NDD:48')],ylab = "Prob (Sample is IPD)",col=cols,cex.axis = 1,cex.lab=1.2)









