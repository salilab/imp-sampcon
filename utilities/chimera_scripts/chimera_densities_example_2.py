import chimera
from chimera import openModels
from chimera import runCommand
import sys, os
from chimera.tkgui import saveReplyLog

# Names of proteins/domains for which we have created densities
prots = ["A", "B", "C", "D",\
         "E", "F", "G", "H",\
         "I", "J", "K",\
         "L", "M"]

# Volume of each protein/domain
vol = {"A":29058,\
       "B":28001,\
       "C":26871,\
       "D":25900,\
       "E":49517,\
       "F":58818,\
       "G":68819,\
       "H":25964,\
       "I":33266,\
       "J":34865,\
       "K":39706,\
       "L":44009,\
       "M":29922}


stdev = {"A":11,\
         "B":9,\
         "C":9,\
         "D":9,\
         "E":12,\
         "F":12,\
         "G":12,\
         "H":12,\
         "I":12,\
         "J":12,\
         "K":12,\
         "L":12,\
         "M":12}

runCommand('windowsize 1024 1024')
runCommand('set bgcolor white')

num=0
for i in range(len(prots)):
    print i, num
    runCommand('open ' + "./data/" + prots[i]+".situs")
    runCommand('color #646464 #'+str(i))
    if i < 4:
        volSAXS= 2.0 * vol[prots[i]]
    else:
        volSAXS= 2.0 * vol[prots[i]]
    runCommand('volume #'+str(i)+' step 1 transparency 0.50 encloseVolume '+str(volSAXS))
    num += 1

for i in range(len(prots)):
    print i, num
    runCommand('open ' + prots[i]+".mrc")
    runCommand('color #FFCCCC #'+str(num))
    runCommand('volume #'+str(num)+' step 1 transparency 0.50 encloseVolume ' + str(vol[prots[i]]))
    runCommand('move cofr mod #'+str(num))

    runCommand('vop gaussian #' +str(num)+' sd '+str(stdev[prots[i]]))
    num += 1

    runCommand('volume #'+str(num)+' step 1 transparency 0.50 encloseVolume ' + str(vol[prots[i]]))
    runCommand('fitmap #' +str(i)+' #'+str(num)+' metric correlation listFits false search 1000 placement rs')
    runCommand('measure correlation #' +str(i)+' #'+str(num))
    runCommand('measure correlation #' +str(num)+' #'+str(i))

    num += 1
    
saveReplyLog("Results_cross_correlation.txt")
##exit()
