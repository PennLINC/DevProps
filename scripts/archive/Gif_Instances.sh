subj=$1
waveInst=$2
count=0
for f in `ls -v /cbica/projects/abcdfnets/results/PWplots/${subj}/${subj}_PWinst_${waveInst}_frame_*.png | head -n +25`
do
printf -v counts "%02d" $count
ln -s $f /cbica/projects/abcdfnets/results/PWplots/${subj}/frame_${counts}.png
count=`expr $count + 1`
done

ffmpeg \
  -framerate 5 \
  -i  /cbica/projects/abcdfnets/results/PWplots/${subj}/frame_%02d.png \
  -vf scale=600:-1 \
  /cbica/projects/abcdfnets/results/PWplots/${subj}/${subj}_Wave_${waveInst}.gif \
;

rm /cbica/projects/abcdfnets/results/PWplots/${subj}/frame_*.png -f
