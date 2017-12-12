import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,string,numpy

score_file=sys.argv[1]
score_name=sys.argv[2]

scores=numpy.loadtxt(score_file)

plt.hist(scores,bins=100)
plt.xlabel(score_name)
plt.ylabel('P')
plt.title('Histogram of scores')

plt.savefig(score_name+'.jpg',dpi=300)


