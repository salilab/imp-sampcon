from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,string,pandas

score_file=sys.argv[1]
column=sys.argv[2]

scores=pandas.read_csv(score_file, sep=' ')
show = (len(sys.argv) >= 4 and sys.argv[3] == "True")

if column=="all":
    print("Plotting all columns")
    for k in scores.columns:
        print(k, "|| Mean:", scores[k].mean(), "| SD:", scores[k].std())
        plt.hist(scores[k],bins=100)
        plt.xlabel(k)
        plt.ylabel('P')
        plt.title(k)

        plt.savefig(k+'.jpg',dpi=300)
        plt.clf()
elif column in list(scores.columns.values):

    print(column, "|| Mean:", scores[column].mean(), "| SD:", scores[column].std())
    plt.hist(scores[column], bins=100)
    plt.xlabel(column)
    plt.ylabel('P')
    plt.title(column)
    if show:
        plt.show()
    else:
        plt.savefig(column+".jpg", dpi=300)
else:
    print(column, "is not a valid score parameter. Use one of:")
    print(scores.columns.values)


