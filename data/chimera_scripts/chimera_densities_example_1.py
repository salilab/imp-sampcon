import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['A', 'B', 'C', 'D',\
         'E', 'F', 'G', 'H']

# Volume of each protein/domain
vol = {'A':86260,\
       'B':103495,\
       'C':187161,\
       'D':135656,\
       'E':113039,\
       'F':147880,\
       'G':121403,\
       'H':127118
       }
# Color of each protein/domain
col =  {'A':"#800080",\
        'B':"#0000FF",\
        'C':"#FF0000",\
        'D':"#6495ED",\
        'E':"#D3D3D3",\
        'F':"#40E0D0",\
        'G':"#FFFF00",\
        'H':"#00FF00"
        }

totvol = 0.0
for p in prots:
    totvol += vol[p]

runCommand('set bgcolor white')
i=0

# EM map
runCommand('open ../Experimental_EM_map.mrc')
runCommand('color #000000 #'+str(i))
runCommand('volume #'+str(i)+' step 1 transparency 0.50 encloseVolume ' + str(totvol))
i+=1

# Read localization densities of all subunits, both samples together
runCommand('open All_subunits.mrc')
runCommand('color #000000 #'+str(i))
runCommand('volume #'+str(i)+' step 1 transparency 0.50 encloseVolume ' + str(totvol))
i+=1    


# Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1')
    runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str( 0.7 - i / 20.0)) 
    i += 1
    
# Read localization density by component for Sample A
for p in prots:
    runCommand('open Sample_1/LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1')
    runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' hide')
    i += 1

# Read localization density by component for Sample B
for num,p in enumerate(prots):
    runCommand('open Sample_2/LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1')
    runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' hide') 
    
    # Calculate cross-correlation between the density maps from Sample A and Sample B
    runCommand('measure correlation #' +str(i)+' #'+str(num))
    runCommand('measure correlation #' +str(num)+' #'+str(i))
    
    
    i += 1

# If the localization densities need to be superposed with the EM map, then get the fitting coordinates to superpose the two.
#runCommand('fitmap #0 #1 metric correlation listFits false search 1000 placement rs')    
#runCommand('open ./all_models.14999/12_3_3210.rmf3')

# Once we have the fitting coordinates, one can use Chimera turn and move commands to superpose them
runCommand('turn -0.74971363,-0.62782964,0.20918754 141.92658853 center 0.00000000,-81.10629115,264.71468669 coord #0 models #0')
runCommand('move -0.74971363,-0.62782964,0.20918754 437.07289348 coord #0 models #0')
runCommand('center')

runCommand('windowsize 1024 1024')
runCommand('turn x 1.0 90')
runCommand('wait 10')
runCommand('move x 64')
runCommand('wait 10')
runCommand('scale 2.5')

saveReplyLog("Results_cross_correlation.txt")

'''
# Run commands to recover same view as first part
runCommand('windowsize 1028 1028')
runCommand('~modeldisp #0')
runCommand('scale 0.99 10')
runCommand('wait 10')
runCommand('turn z 1.0 90')
runCommand('wait 90')
runCommand('turn x -1.0 90')
runCommand('wait 90')
runCommand('clip hither -25.0 30')
runCommand('wait 30')
runCommand('scale 1.02 30')
runCommand('wait 30')
runCommand('select #98-100:.F; ~show sel; ~ribbon sel')

'''
