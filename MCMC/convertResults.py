import cPickle as pickle
import csv



if __name__ == '__main__':
	print('Loading pickle.')
	mcmc = pickle.load(open('BEZ252','rb'))
	
	print('Done loading pickle.');
	
	with open('8LV41B.csv', 'wb') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(mcmc.positions)