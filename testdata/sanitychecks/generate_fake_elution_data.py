import pandas as pd
import numpy as np



def main():
	
	#generate 5 clusters of with random elution profiles
 	clusters=[]	
	#for a given cluster generate distributions
	df=pd.DataFrame(columns=range(100))
	q=0
	for i in range(5):
		cluster_n = ''
		mu, sigma = np.random.uniform(20,80), np.random.uniform(0,20) # mean and standard deviation
		for j in range(np.random.randint(5)+1):
			sample = np.random.normal(mu,sigma,500)
			sample = np.clip(sample,0,99)
			sample = np.round(sample)
			unique, counts = np.unique(sample,return_counts=True)
			elution=np.zeros(100)
			for k,u in enumerate(unique):
				elution[u]=counts[k]
			df.ix['p'+str(q)+'_cluster'+str(i)]=elution
			
			cluster_n=cluster_n+'	' + 'p'+str(q)+'_cluster'+str(i)
			q=q+1
		clusters.append(cluster_n)

	pd.DataFrame(clusters).to_csv('./test_clusters.txt',index=False,header=False)
	df.to_csv('./test_elution.tab',sep='\t')
	print df


if __name__ == "__main__":
	main()

