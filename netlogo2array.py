## Analysing Fake Dolphin Data
import sys
from scipy import *
import matplotlib.pyplot as plt

def netlogo2array(filename):
	netlogo = open(filename)
	BIDdata = []
	# Create a Collection of the Paths
	for headingList in netlogo:
		BIDdata.append(headingList)
	for i in range(0,len(BIDdata)):
		BIDdata[i] = BIDdata[i][1:-2]
		BIDdata[i] = BIDdata[i].split(' ')
	# Find the longest path and make them all the same length
	# (copy the final position of all the others)
	shortestPath = min([size(BIDdata[i]) for i in range(0,size(BIDdata))])
	for i in range(0,len(BIDdata)):
		for j in range(shortestPath,len(BIDdata[i])):
			BIDdata[i].pop()
			#(BIDdata[i][-1])

	BIDdata = [map(float,angleList) for angleList in BIDdata]
	BIDdata = array(BIDdata)
	save(filename,BIDdata)

def plotAngles(filename):
	dolphinArray = load(filename)
	xvals = shape(dolphinArray)[1]
	plt.figure()
	[plt.plot(path) for path in dolphinArray]
	plt.axis([0,xvals,0,360])
	plt.xlabel('Time Steps')
	plt.ylabel('Heading (degs)')
	plt.title('Dolphins!')

if __name__ == '__main__':
	filename = str(sys.argv[1])
	netlogo2array(filename)

# #Crazy polar coordinat plot!
# fig = figure()
# ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
# plot(dolphinArray[1][1500:1751]/(2*pi),arange(0,len(dolphinArray[1][1500:1751])))
# show()