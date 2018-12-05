import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
# rcParams['font.family'] = 'Latin Modern Roman'
# import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

import seaborn as sns

# sns.set()


def load_single_bbf(bbf, ES_type):
	fopt_history = np.array(pd.read_csv(f'./Data_fopt/{bbf}_fopt_history_{ES_type}.csv'))
	xopts = np.array(pd.read_csv(f'./Data_xopt/{bbf}_xopt_{ES_type}.csv'))

	return fopt_history, xopts

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

	return mean_progression, std_progression, x_data

def makeplot(mean_progression, std_progression, generations, fopt_history, bbf, ES_type, plotloc = './Plots/'):
	#make an array indicating the number of evaluations
	evals = np.linspace(1, 10000, len(generations))

	# gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
	gs = plt.GridSpec(2, 1, height_ratios=[4, 1], hspace = 0)

	ax1 = plt.subplot(gs[0])

	#plot every individual run
	for i in range(fopt_history.shape[0]):
		plt.plot(generations, fopt_history[i], linewidth = 0.3, color = 'grey')

	#plot the average runs
	plt.plot(generations, mean_progression, color = '#2934A3')
	#also plot the standard deviation
	# plt.fill_between(generations, mean_progression + std_progression, mean_progression - std_progression, alpha = 0.3)

	#obtain limits
	xlims = plt.xlim()

	#plot the optimal values per run
	fopt = np.min(fopt_history, axis = 1)
	for i in range(len(fopt)):
		plt.plot([401, 406], [fopt[i], fopt[i]], color = 'red', linewidth = 1)
	#plot also the mean of these optimal values
	mean_fopt = np.mean(fopt)
	plt.plot([410, 415], [mean_fopt, mean_fopt], color = '#2934A3', linewidth = 2)

	plt.ylabel('Function value')
	plt.title(f'Progression of evolutionary algorithm for {bbf.upper()}')
	if int(bbf[3]) < 4:
		plt.yscale('log')
	#fix limits
	plt.xlim(xlims)
	# plt.ylim((7e2, 3e3))
	plt.grid(linestyle = '--')


	#plot the standard deviation in a separate plot
	ax2 = plt.subplot(gs[1], sharex = ax1)
	plt.plot(generations, np.std(fopt_history, axis = 0), color = '#647AC9')

	ylims = plt.ylim()
	maxylim = ylims[1]

	xticksloc = np.linspace(0, maxylim, 4)[:3]
	fnogwat = ticker.ScalarFormatter(useOffset = False, useMathText = True)
	g = lambda x, pos : '${}$'.format(fnogwat._formatSciNotation('%1.1e' % x))
	fmt = ticker.FuncFormatter(g)

	xticks_string = []
	for xtick in xticksloc:
		xticks_string.append('{}'.format(fmt(xtick)))

	plt.yticks(xticksloc, xticks_string)
	# ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0e'))

	plt.xlim(xlims)
	plt.xlabel('Number of generations', labelpad = 1)
	plt.grid(linestyle = '--')

	#add axis
	fig = plt.gcf()
	ax3 = fig.add_axes((0.156,0.1,0.71,0.0))
	ax3.yaxis.set_visible(False)
	ax3.set_xticks(np.linspace(0, 10000, 6))
	ax3.set_xlabel('Number of evaluations', labelpad = 1)

	fig.subplots_adjust(bottom=0.2)

	plt.savefig(f'{plotloc}{bbf}_ES_progression_{ES_type}.png', dpi = 300, bbox_inches = 'tight')
	# plt.show()
	plt.close()


ES_type = '5_25'
for b in range(1,6):
	bbf = f'bbf{b}'
	print('')
	print(bbf)

	fopt_history, xopts = load_single_bbf(bbf, ES_type)
	mean_progression, std_progression, generations = process_data(fopt_history)
	makeplot(mean_progression, std_progression, generations, fopt_history, bbf, ES_type)