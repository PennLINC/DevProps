% read in mean directions difference from PGG per face: circ mean over time then over subjs to yield these
% from scripts "MeanAngle_EachFace_xTrsxSubjs_L" and _R
left=readtable('~/data/MeanMeanPGGDif_L.csv','ReadVariableNames',false);
right=readtable('~/data/MeanMeanPGGDif_R.csv','ReadVariableNames',false);
% translate and reflect angles to yield vectors of aligned vs orthogonal (0-90 degrees)
writetable(table((abs(abs(str2double(table2array(left(2,:))))-pi)*-1)+pi),'~/data/MeanMeanReflect90L.csv')
writetable(table((abs(abs(str2double(table2array(right(2,:))))-pi)*-1)+pi),'~/data/MeanMeanReflect90R.csv')
% let me explain: let's first strip out the matlab nonsense (data type conversions)

% (abs(abs(left(2,:))-pi)*-1)+pi

% next, let's consider the starting values. They are in the range of -2pi to 2pi. It's mean angle - pgg for each face.
% our goal here is to consider their alignment versus orthogonality.
% Perfect alignment occurs at -2pi, 0, and 2pi.
% "Perfect" orthogonality occurs at -pi and pi.
% so we have to "translate" 4 segments of this range to the same space. -2pi to -pi = 0 to -pi = 0 to pi = 2pi to pi 
% note in each defined segment, alignment is the left value, and orthogonality is the right value

% now lets work from the insider of the ()'s out:
% left(2,:) just selects the angular values. heed not the 2.
% abs(left(2,:)) translates -pi:0 to 0:pi, maintaining 0 = aligned, pi = orthogonal.
% it also translates -2pi:-pi into pi:2pi, maintaining orthogonal = low and aligned = high.
% so half our segments are halfway translated: they share appropraite vales with one OTHER segment.
% but 0:pi and pi:2pi are not equivalent. Higher values are more orthogonal for 0:pi, but more aligned for pi:2pi.
% So instead of 4 sets of angular descriptors we have 2. All in the range of 0 to 2pi, rather than -2pi to 2pi.

% the next ()'s after abs(left(2,:)) are abs(   abs(left(2,:))   -pi).
% This is a little less intuitive: we are explicitly subtracting pi so that we may zero the data at the TRUE VALUE of pi
% this is because we need to reflect half of the angular distribution across pi, so that all ascending numbers mean greater orthogonality. 
% So to reflect, that is, to mirror the distribution, the distribution is centered/zeroed on pi (rather than 0) and abs()'ed again. 

% now numbers that were from 0:pi are mirrored: whereas > # meant > orthogonality, it now means > alignment.
% AND it holds true for numbers that were pi to 2pi. > # now = means > alignment across the board
% however, we lost the true center of the distribution by playing with centering/zero.

% the last step corrects this. (    abs(abs(left(2,:))-pi)    *-1)+pi. By Reflecting ALL numbers (*-1), we maintain their now-shared
% interpretation of ascent/descent. So no matter where we are on the number line, numerical ascent and descent have the same interpretation across.
% We can make this uniform interpretation with one more step. +pi takes 0, which really means "pi away from 0" prior to our translation-for-reflection, and makes it pi. Because the highest values were set to -pi through the *-1 reflection, these now become 0.

% This final makeup is most desireable, because the x-axis now means "radians away from the PGG, irrespective of whether we are talking about hierarchical ascent or descent

% thank you for coming to my ted talk 
