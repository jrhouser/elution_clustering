import pandas as pd
import numpy as np



def main():
	
	#generate 5 clusters of with random elution profiles
	
	#for a given cluster generate distributions
	df=pd.DataFrame(columns=range(100))
	for i in range(5):
	
		mu, sigma = np.random.uniform(0,100), np.random.uniform(0,10) # mean and standard deviation
		for j in range(np.random.randint(5)):
			sample = np.random.normal(mu,sigma,500)
			sample = np.clip(sample,0,99)
			unique, counts = np.unique(sample,return_counts=True)
			elution=np.zeros(100)
			for k,u in enumerate(unique):
				elution[u]=counts[k]
			df.ix['p'+str(j)+'_cluster'+str(i)]=elution


	df.to_csv('./test_elution.tab',sep='\t')
	print df


if __name__ == "__main__":
	main()

