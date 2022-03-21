import imageio
import numpy as np    
import sys
subj = sys.argv[1]
# get number of waves
NumW = sys.argv[2]
# filepath for Pwave gifs
wavedir='/cbica/projects/abcdfnets/results/PWplots/' + str(subj) + '/'
# load in a gif object for each wave, using the contentious global method
for x in range(int(NumW)):
	gifFP=wavedir + str(subj) + '_Wave_' + str(x) + '.gif'
	globals()['Gif%s' % x] = imageio.get_reader(gifFP)

# number of frames in current workflow is 25
number_of_frames = 25

#Create writer object
outputfn=wavedir + str(subj) + '_AllWaves.gif'
new_gif = imageio.get_writer(outputfn)
# merge em with numpy
for frame_number in range(number_of_frames):
    # initialize string array of waves
    stringArray=np.array([])
    # repeat globals method to extract the frame from each wave
    for x in range(int(NumW)):
        globals()['img%s' % x]=globals()['Gif%s' % x].get_next_data()
        stringArray=np.append(stringArray,['img' + str(x) + ', '])
    # remove , from last entry
    lastEnt=stringArray[x]
    stringArray[x]=lastEnt[:-2]
    # behold, the ugliest line of code in the whole project
    new_image = np.hstack(((eval(eval(str(stringArray))[0]))))
    new_gif.append_data(new_image)

# one more global variable loop
for x in range(int(NumW)):
    globals()['Gif%s' % x].close()

new_gif.close()
