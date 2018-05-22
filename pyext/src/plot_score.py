import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,string,numpy

score_file=sys.argv[1]
score_name=sys.argv[2]

scores=numpy.loadtxt(score_file)
print "Mean",score_name,numpy.mean(scores)
print "Standard deviation",score_name,numpy.std(scores)

print "Mean minus std_dev", score_name, numpy.mean(scores) - numpy.std(scores)
print "Mean minus half std_dev", score_name, numpy.mean(scores) - 0.5*numpy.std(scores)

plt.hist(scores,bins=100)
plt.xlabel(score_name)
plt.ylabel('P')
plt.title('Histogram of scores')

plt.savefig(score_name+'.jpg',dpi=300)


