import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
# rcParams['font.family'] = 'Latin Modern Roman'
import matplotlib.gridspec as gridspec
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
		plt.plot([generations[-1]+1, generations[-1]+6], [fopt[i], fopt[i]], color = 'red', linewidth = 1)
	#plot also the mean of these optimal values
	mean_fopt = np.mean(fopt)
	plt.plot([generations[-1]+10, generations[-1]+15], [mean_fopt, mean_fopt], color = '#2934A3', linewidth = 2)

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

def make_grid_plots(mean_progression, std_progression, generations, fopt_history, bbf, ES_type, gs, row_main_ylims, row_error_ylims, ES_id):

	bbf_id = int(bbf[3])

	#make an array indicating the number of evaluations
	evals = np.linspace(1, 10000, len(generations))

	##### make the main plot
	ax1 = plt.subplot(gs[0])

	#plot every individual run
	for i in range(fopt_history.shape[0]):
		plt.plot(generations, fopt_history[i], linewidth = 0.3, color = 'grey')

	#plot the average runs
	plt.plot(generations, mean_progression, color = '#2934A3')

	#obtain x limits
	xlims = plt.xlim()

	#plot the optimal values per run
	fopt = np.min(fopt_history, axis = 1)
	for i in range(len(fopt)):
		plt.plot([generations[-1]+1, generations[-1]+6], [fopt[i], fopt[i]], color = 'red', linewidth = 1)
	#plot also the mean of these optimal values
	mean_fopt = np.mean(fopt)
	plt.plot([generations[-1]+10, generations[-1]+15], [mean_fopt, mean_fopt], color = '#2934A3', linewidth = 2)

	#only add y label if we are plotting the first column
	if ES_id == 0:
		plt.ylabel('Function value')
	# plt.title(f'Progression of evolutionary algorithm for {bbf.upper()}')
	if bbf_id < 4:
		plt.yscale('log')
	#fix limits
	plt.xlim(xlims)
	plt.ylim(row_main_ylims[bbf_id-1])

	#disable x ticks but keep grid
	frame1 = plt.gca()
	frame1.axes.xaxis.set_ticklabels([])
	#disable y ticks if needed
	if ES_id > 0:
		frame1.axes.yaxis.set_ticklabels([])
	plt.grid(linestyle = '--')


	##### plot the standard deviation in a separate plot
	ax2 = plt.subplot(gs[1], sharex = ax1)
	plt.plot(generations, np.std(fopt_history, axis = 0), color = '#647AC9')

	#fix the y ticks
	plt.ylim(row_error_ylims[bbf_id-1])
	maxylim = row_error_ylims[bbf_id-1][1]

	yticksloc = np.linspace(0, maxylim, 4)[:3]
	fnogwat = ticker.ScalarFormatter(useOffset = False, useMathText = True)
	g = lambda x, pos : '${}$'.format(fnogwat._formatSciNotation('%1.1e' % x))
	fmt = ticker.FuncFormatter(g)

	yticks_string = []
	for ytick in yticksloc:
		yticks_string.append('{}'.format(fmt(ytick)))
	plt.yticks(yticksloc, yticks_string)

	split_ES_type = ES_type.split('_')

	#disable x ticks but keep grid if needed
	frame1 = plt.gca()
	if bbf_id < 5:
		frame1.axes.xaxis.set_ticklabels([])
	else:
		plt.xticks(np.linspace(0, 10000, 6)/(int(split_ES_type[0]) + int(split_ES_type[1])))
	#disable y ticks if needed
	if ES_id > 0:
		frame1.axes.yaxis.set_ticklabels([])
	plt.grid(linestyle = '--')

	plt.xlim(xlims)
	if bbf_id == 5:
		plt.xlabel('Number of generations', labelpad = 1)

		##### add the bottom extra x axis
		ax3 = ax2.twiny()
		# ax3 = f.add_axes((0.156,0.02,0.71,0.0))
		ax3.yaxis.set_visible(False)
		ax3.xaxis.set_ticks_position('bottom')
		ax3.xaxis.set_label_position('bottom')
		ax3.spines['bottom'].set_position(('outward', 36))
		ax3.set_xticks(np.linspace(0, generations[-1], 6))
		ax3.set_xticklabels(np.linspace(0, 10000, 6))
		ax3.set_xlabel('Number of evaluations', labelpad = 1)

#define the different data sets which we will use for the grid
ES_types = ['3_21', '5_25', '5_35']
# ES_types = ['3_21']
bbfs = []
for b in range(1,6):
	bbfs.append(f'bbf{b}')

#define the gridspec
default_figsize = [6.4, 4.8]
f = plt.figure(figsize = (default_figsize[0]*len(ES_types), default_figsize[1]*len(bbfs)))
gs0 = gridspec.GridSpec(len(bbfs), len(ES_types), figure=f, hspace = 0.05, wspace = 0)

#manually set the y lims for every row
row_main_ylims = [
		[20, 4e5],
		[100, 2e6],
		[30, 3e9],
		[40, 150],
		[67, 69.5]
		]
row_error_ylims = [
		[0, 6e4],
		[0, 2.5e5],
		[0, 3e8],
		[0, 250],
		[0, 0.8]
		]

gs_counter = 0
for bbf in bbfs:
	for ES_id, ES_type in enumerate(ES_types):
		print('')
		print(bbf)

		#make a subgrid for the progression plot and std plot
		gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[4, 1], hspace = 0, subplot_spec=gs0[gs_counter])

		#load and process data
		fopt_history, xopts = load_single_bbf(bbf, ES_type)
		mean_progression, std_progression, generations = process_data(fopt_history)

		make_grid_plots(mean_progression, std_progression, generations, fopt_history, bbf, ES_type, gs00, row_main_ylims, row_error_ylims, ES_id)

		gs_counter += 1

#give the extra x-axes at the bottom more space
f.subplots_adjust(bottom=0.03)

plt.savefig(f'./Plots/ES_progression_comma.png', dpi = 200, bbox_inches = 'tight')
plt.close()