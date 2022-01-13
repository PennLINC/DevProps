import os
import subprocess
import time
import datetime
# grab subjects list
my_file = open("/cbica/projects/pinesParcels/PWs/hcpd_subj_list.txt", "r")
content = my_file.read()
content_list = content. split("\n")
# remove last line (blank)
content_list.pop()
# feed em' in as subjects
subjects = content_list
# while there are more than 0 subjects left to run
while len(subjects)>0:
  # grab qsub info, get number of jobs being run
  qstat = subprocess.check_output(['qstat'],shell=True).decode().split('/bin/python')[0]
  que = len(qstat.split('\n'))-3
  # if we are using less than 6 job slots (one is occupied by this script)
  if que < 6:
    newsub = subjects.pop()
    # test if output file exists
    #OF='/cbica/projects/pinesParcels/results/PWs/Proced/' + str(newsub) + '/' + str(newsub) + '_AngDistMat.mat'
    #if not os.path.exists(OF):
    subprocess.run(["qsub","-l","h_vmem=15G,s_vmem=14G","qsubMatlab_PGG_ANGD.sh",newsub])
