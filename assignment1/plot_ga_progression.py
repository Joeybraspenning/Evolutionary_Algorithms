import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Latin Modern Roman'

import seaborn as sns

# sns.set()

def load_MC():
	data = []

	#loop over all the csv files
	for i in range(1, 21):
		filename = f'fitness_mc_{i}.csv'
		datastring = ''
		with open(filename, 'r') as f:
			for line in f:
				datastring += line

		data.append(datastring.split(','))

	return data

def load_GA():
	filename = 'Fitness_ga_2.csv'

	data = []
	with open(filename, 'r') as f:
		for line in f:
			data.append(line.split(','))

	return data

def process_data(data):
	data = np.array(data).astype(float)
	#average and standard deviation of this timeseries
	mean_progression = np.mean(data, axis = 0)
	std_progression = np.std(data, axis = 0)

	x_data = np.arange(1, data.shape[1]+1)

	#print statistics
	best_results = np.min(data, axis = 1)
	print(f'Mean: {np.mean(best_results)}')
	print(f'STD: {np.std(best_results)}')

	return mean_progression, std_progression, x_data, data

def plot_MC(mean_progression, std_progression, x_data):
	plt.plot(x_data, mean_progression)
	plt.fill_between(x_data, mean_progression + std_progression, mean_progression - std_progression, alpha = 0.3)

	plt.yscale('log')
	plt.ylim((1e3, 1e4))
	plt.grid(linestyle = '--')

	plt.xlabel('Number of evaluations')
	plt.ylabel('Power loss [kW]')
	plt.title(r'Average progression of Monte Carlo optimizer with 1$\sigma$ error band')

	plt.savefig('MC_progression.png', dpi = 300, bbox_inches = 'tight')
	plt.show()

def plot_GA(mean_progression, std_progression, generations, data):
	#make an array indicating the number of evaluations
	evals = np.linspace(1, 10000, len(generations))

	#plot the best result from the literature
	plt.plot([-10,generations[-1]+10], [869.7271, 869.7271], color = 'red', linewidth = 0.8)

	for i in range(data.shape[0]):
		plt.plot(generations, data[i], linewidth = 0.3, color = 'grey')

	plt.plot(generations, mean_progression)
	plt.fill_between(generations, mean_progression + std_progression, mean_progression - std_progression, alpha = 0.3)

	plt.xlabel('Number of generations')
	plt.ylabel('Power loss [kW]')
	plt.title(r'Average progression of genetic algorithm with 1$\sigma$ error band')

	plt.yscale('log')
	plt.ylim((7e2, 3e3))
	plt.grid(linestyle = '--')

	#add axis
	fig = plt.gcf()
	ax2 = fig.add_axes((0.18,0.1,0.667,0.0))
	ax2.yaxis.set_visible(False)
	ax2.set_xticks(np.linspace(0, 10000, 6))
	ax2.set_xlabel('Number of evaluations')

	fig.subplots_adjust(bottom=0.2)

	plt.savefig('GA_progression_mu_lambda.png', dpi = 300, bbox_inches = 'tight')
	plt.show()

# data = load_MC()
# mean_progression, std_progression, x_data, __ = process_data(data)
# plot_MC(mean_progression, std_progression, x_data)

data = load_GA()
mean_progression, std_progression, generations, data = process_data(data)
plot_GA(mean_progression, std_progression, generations, data)