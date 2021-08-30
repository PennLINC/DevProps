import numpy as np
import glob
import re

# different types of tasks
tasks=['rest','SST','nback','mid']

# parent directory
parentfp='/cbica/projects/abcdfnets/results/wave_output/*/*NLFcmat'

# make a psuedo fc matrix for colnames
pseuFC=np.chararray((17,17),itemsize=11,unicode=True)
# set colnames
for n in range(17):
	for Othern in range(17):
		# +1 because these colnames will be used outside of python
		pseuFC[n,Othern]=str(n+1) + '_' + str(Othern+1)

# get number of subjects from GoodT1s run list
my_file = open("/cbica/projects/abcdfnets/GoodT1s.txt", "r")
content = my_file.read()
content_list = content. split("\n")

# initialize master fc array, 153 fc feats for each task + 1 to cols for subjnames col, +1 to rows for colnames
masterfc=np.empty((len(content_list)+1, (153*4) + 1),dtype=object)
# first column is subj name
masterfc[1:,0]=content_list
# subj name colname
masterfc[0,0]='subjectkey'

# for all tasks
for t in range(len(tasks)):
	print(tasks[t])
	# pattern corresp. to these fcmats
	pattern=parentfp + tasks[t] + '.npy'
	# fc mats corresponding to this task
	fcmatsFns=glob.glob(pattern)
	# initialize output array, +1 for colnames, 153 covers the full fcmat without redundancy
	outArray=np.zeros((len(content_list)+1,153),dtype=object)
	# load in each subject
	for s in range(len(fcmatsFns)):
		thefile=np.load(fcmatsFns[s])
		thefields=re.split('/',fcmatsFns[s])
		thesubj=thefields[6]
		# take upper triangle of their fc matrix
		FCvec=thefile[np.triu_indices(17)]
		# find where this subject occurs in full subj list
		Sindex=content_list.index(thesubj)
		# +1 because first row will be colnames
		# round decimals to reduce filesize
		outArray[Sindex+1,:]=np.round(FCvec,decimals=3)
	# set colnames with same indexing
	ColNames=pseuFC[np.triu_indices(17)]
	# specify task
	ColNames=np.char.add(ColNames,tasks[t])
	outArray[0,:]=ColNames
	# merge into one df, starting past 1st col (subj names)
	masterfc[:,((153*t)+1):((153*(t+1)+1))]=outArray

# remove 0 rows for now (subjs still running)
masterfc = masterfc[~np.all(masterfc[:,1:] == 0, axis=1)]

# save dataframe
saveFn='/cbica/projects/abcdfnets/results/masterfc.csv'
np.savetxt(saveFn,masterfc,delimiter=',',fmt='%s')
