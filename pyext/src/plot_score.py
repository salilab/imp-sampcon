from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,string,pandas


def plot_column(scores, column, show):
    print(column, "|| Mean:", scores[column].mean(),
          "| SD:", scores[column].std())
    plt.hist(scores[column], bins=100)
    plt.xlabel(column)
    plt.ylabel('P')
    plt.title(column)
    if show:
        plt.show()
    else:
        plt.savefig(column+".png", dpi=300)
    return plt


score_file=sys.argv[1]
column=sys.argv[2]

scores=pandas.read_csv(score_file, sep=' ')
show = (len(sys.argv) >= 4 and sys.argv[3] == "True")

if column=="all":
    print("Plotting all columns")
    exclude = frozenset(('Model_index', 'Replica_id', 'Frame_id', 'Run_id'))
    for k in scores.columns:
        if k not in exclude:
            plt = plot_column(scores, k, show=False)
            plt.clf()
elif column in scores.columns.values:
    plot_column(scores, column, show=show)
else:
    print(column, "is not a valid score parameter. Use one of:")
    print(scores.columns.values)


