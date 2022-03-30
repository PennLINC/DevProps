import os
import subprocess
import time
import datetime
# grab subjects list
my_file = open("/cbica/projects/pinesParcels/PWs/G600TRs.txt", "r")
content = my_file.read()
content_list = content. split("\n")
# pause to let extract BUTD finish
# time.sleep(20000)
# remove last line (blank)
content_list.pop()
# feed em' in as subjects
subjects = content_list
# while there are more than 0 subjects left to run
while len(subjects)>0:
  # grab qsub info, get number of jobs being run
  qstat = subprocess.check_output(['qstat'],shell=True).decode().split('/bin/python')[0]
  que = len(qstat.split('\n'))-3
  # if we are using less than 7 job slots (one is occupied by this script)
  # ayyyy we done wit 7 jobs atta time we str8 flossin' out here at 35 now 
  if que < 40:
    # see if it is the weekend, 0, 1, 2, 3, and 4 are weekday, 5 and 6 are weekend
    weekno = datetime.datetime.today().weekday()
    # see if it is before 9 or after 5 
    Hour = time.localtime().tm_hour 
    # endgame file to see if subj ran
    # if weekend OR after 6 PM OR before 9 AM
    if weekno > 4 or Hour < 9 or Hour > 17 :
      newsub = subjects.pop()
      # submit job (if conditions are met)
      OpFile='/cbica/projects/pinesParcels/results/PWs/Proced/' + str(newsub) + '/' + str(newsub) + '_CSpunDips4_spoon.csv'
      if not os.path.exists(OpFile):
        print(newsub)
        subprocess.run(["qsub","qsubMatlab.sh",newsub])
      # added this to run 3 subjs (1 slot for this job) during ON hours
    elif que < 40:
      newsub = subjects.pop()
      OpFile='/cbica/projects/pinesParcels/results/PWs/Proced/' + str(newsub) + '/' + str(newsub) + '_CSpunDips4_spoon.csv'
      if not os.path.exists(OpFile):
      # submit job (if conditions are met)
        print(newsub)
        subprocess.run(["qsub","qsubMatlab.sh",newsub])
