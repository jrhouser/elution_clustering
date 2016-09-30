
import pandas as pd
import numpy as np
import sklearn.mixture as sk
import os
import sys
import pickle
sys.path.append('/media/jrh94/HD1/Documents/protein_complex_maps/protein_complex_maps/')
import complex_comparison as cc
import argparse


def main():


	parser = argparse.ArgumentParser(description="Cluster elution data using Gaussian Mixture Models")
	parser.add_argument("--elution_folder", action="store", dest="folder", required=True, 
		                            help="folder containing raw elution data in tab format with first column containing protein labels")
	
	parser.add_argument("--clustN", action="store", dest="clustN", required=True, 
		                            help="number of clusters (k). need for input to the GMM. for trying multiple cluster numbers seperate by a space (e.g. 10 100 1000)")
	
	parser.add_argument("--standard_filenames", action="store", dest="standard_filenames", required=True, 
		                            help="file names containing clusters for comparison. can handle multiple comparisons. need to be tab seperated (e.g. cluster1.tab cluster2.tab)")


	parser.add_argument("--base_output_filename", action="store", dest="base_output_filename", required=True, 
		                            help="the base output filename for the results. it should not contain an extension (e.g. ~/folder/cluster_results_filename")


	parser.add_argument("--cluster_bayes", action="store_true", dest="cluster_bayes", required=False, default=False, 
		                            help="use this flag if you want to cluster based upon the Baysian GMM method")

	args = parser.parse_args()
	clustN = args.clustN.split()
	#print clustN
	standard_complexes=[]
	
	print args.cluster_bayes

	for standard_filename in args.standard_filenames.split():
		#print standard_filename
		standard_complexes.append(generate_comparison_complexes(standard_filename))

	df_concat = concatenate_elutions(args.folder)

	X=df_concat.fillna(0).values

	metric,cluster_prediction_out = cluster_gmm_core(X,clustN,args.base_output_filename,standard_complexes,args.standard_filenames.split(),df_concat,args.cluster_bayes)

	


def concatenate_elutions(folder):

	files = os.listdir(folder)

	data=[]
	prots=[]
	df_concat=pd.DataFrame()
	for f in files:

	    df=pd.read_table(folder + f ).set_index('Unnamed: 0')

	    df=df.ix[df.sum(axis=1)>0]
	    X=df.fillna(0).values
	    
	    Xnorm=np.max(X,axis=1).reshape((len(X),1))*np.ones((1,len(X[0])))
	    X=np.divide(X,Xnorm)
	    df=pd.DataFrame(X,index=df.index,columns=df.columns)
	    data.append(X)
	    prots.append(df.index)
	    df_concat=pd.concat([df_concat,df],axis=1)
    	
	df_concat=df_concat.fillna(0)
	
	return df_concat


def generate_comparison_complexes(standard_filename):
	

	#kdrew: read gold standard into list
	gstd_file = open(standard_filename,"rb")
	gold_standard_complexes = []
	for line in gstd_file.readlines():
	    gold_standard_complexes.append(line.split())

	return gold_standard_complexes




def format_cluster(y,names):
    
    U=np.unique(y)
    names=np.array(names)
    
    cmplx=[]
    for u in U:
        cmplx.append(names[np.where(y==u)[0]])
        
    return cmplx
        



def cluster_gmm_core(X,clustN,outputfilename,standard_complexes,standard_filenames,df_concat,cluster_bayes):
    cluster_prediction_out=[]
    metric = pd.DataFrame(index=[str(x) for x in clustN])
    

    
    for c in clustN:
     	
	if cluster_bayes:
		model=sk.BayesianGaussianMixture(int(c))
		outputfilename=outputfilename+'_bayesGMM'
	else:
        	model=sk.GaussianMixture(int(c))
        model.fit(np.nan_to_num(X))
        y = model.predict(np.nan_to_num(X))
        #print y
        
        cluster_prediction=format_cluster(y,df_concat.index)
	try:        
		cluster_prediction=cluster_prediction[np.where([len(cp)>1 for cp in cluster_prediction])[0]] # get rid of clusters with only 1 protein as these are the trivial boring clusters.
        except:
		        print 'removing clusters of one failed'
		        
		
	print cluster_prediction
	cluster_prediction_out.append(cluster_prediction)
        pickle.dump(cluster_prediction_out,open(outputfilename+'_clusters.p','wb'))

	for i,gold_standard_complexes in enumerate(standard_complexes):
		cplx_comparison = cc.ComplexComparison(gold_standard_complexes, cluster_prediction) 
		cplx_comparison_normalize = cc.ComplexComparison(gold_standard_complexes, cluster_prediction, normalize_by_combinations=True, pseudocount=0.00001) 
	    	
		compare_type=standard_filenames[i].split('.')[0].split('/')[-1]
		#print cplx_comparison
	    
		metric.loc[str(c)+ ' ' + compare_type,'bic']=model.bic(np.nan_to_num(X))
		metric.loc[str(c) + ' ' + compare_type,'aic '+compare_type]=model.aic(np.nan_to_num(X))
		print c,model.bic(np.nan_to_num(X)),cplx_comparison.acc()
		metric.loc[str(c)+ ' ' + compare_type,'acc'] = cplx_comparison.acc() 
		metric.loc[str(c)+ ' ' + compare_type,'sensitivity'] = cplx_comparison.sensitivity() 
		metric.loc[str(c)+ ' ' + compare_type,'ppv'] = cplx_comparison.ppv() 
		metric.loc[str(c) + ' '+ compare_type,'mmr'] = cplx_comparison.mmr()
		
		try:
		    metric.loc[str(c)+ ' '+compare_type,'precision_recall_product'] = cplx_comparison.precision_recall_product() 
		    ccmm = cplx_comparison.clique_comparison_metric_mean()
		    metric.loc[str(c) + ' '+compare_type,'clique_precision_mean'] = ccmm['precision_mean']
		    metric.loc[str(c) + ' '+compare_type,'clique_recall_mean'] = ccmm['recall_mean']
		
		    ccmm_normalize = cplx_comparison_normalize.clique_comparison_metric_mean()
		    metric.loc[str(c)+' ' +compare_type,'clique_precision_mean_normalize'] = ccmm_normalize['precision_mean']
		    metric.loc[str(c) + ' '+compare_type,'clique_recall_mean_normalize'] = ccmm_normalize['recall_mean']
		    ccmm_normalize_weighted = cplx_comparison_normalize.clique_comparison_metric_mean(weighted=True)
		    metric.loc[str(c)+ ' ' +copmare_type,'clique_precision_mean_normalize_weighted'] = ccmm_normalize_weighted['precision_mean']
		    metric.loc[str(c) + ' ' +compare_type,'clique_recall_mean_normalize_weighted'] = ccmm_normalize_weighted['recall_mean']
		except Exception as e:
		        print e
		        continue
        
       
    
        metric.to_csv(outputfilename+'.csv')
        

    return metric,cluster_prediction_out
    



if __name__ == "__main__": 
	main()



