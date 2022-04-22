% read in mean directions difference from PGG per face: circ mean over time then over subjs to yield these
% from scripts "MeanAngle_EachFace_xTrsxSubjs_L" and _R
left720=readtable('~/data/MeanMeanPGGDif_L.csv','ReadVariableNames',false);
right720=readtable('~/data/MeanMeanPGGDif_R.csv','ReadVariableNames',false);

% our goal here is to consider their alignment versus orthogonality.
% Perfect alignment occurs at -2pi, -pi, 0, pi, and 2pi.
% "Perfect" orthogonality occurs at -3/2 pi, -1/2 pi, 1/2 pi and 3/2 pi.
% so we have to "translate" 8 segments of this range to the same space. -2pi to -3/2 pi = 0 to 1/2 pi, etc.
% the final "language" desired here is "low values = aligned with PGG, high values = orthogonal to it"

% translating -2pi:0 to 0:2pi is the easiest. Lets start there.
left360=abs(str2double(table2array(left(2,:))));
right360=abs(str2double(table2array(right(2,:))));
% great: now there are only 4 kinds of segments to translate rather than 8. What was -2pi:-3/2pi is now equiv to 3/2pi:2pi, etc.
% our range is also halved: now 0 to 2pi.

% this gets a little bit trickier now, because we need to center the distribution on 0 a few times in order to reflect it with abs()
% reflecting allows us to convert segments of the circular distribution into equivalence, but requires correction for these transformations

% so the first time we prep to reflect, we will subtract pi from every value so as to center the distribution on pi for reflecting
leftAbs_minpi=left360-pi;
rightAbs_minpi=right360-pi;

% then we reflect it with abs()
leftAbs_minpi_abs=abs(leftAbs_minpi);
rightAbs_minpi_abs=abs(rightAbs_minpi);

% we've successfully collapsed half of our distribution again, but now our distribution is backwards: lowest values are high
leftAbs_minpi_abs_neg=leftAbs_minpi_abs*-1;
rightAbs_minpi_abs_neg=rightAbs_minpi_abs*-1;

% now our distribution is collapsed appropriately and in the right order, but on the wrong side of the tracks. restore the initial pi we subtracted.
left180=leftAbs_minpi_abs_neg+pi;
right180=rightAbs_minpi_abs_neg+pi;

% tight. Now we only have two segments left to translate, and our data is in the range 0 to pi. 0:1/2pi is an aligned-orthogonal ascent, but 1/2pi:pi is an orthogonal-aligned ascent. One more reflection needed, across 1/2pi rather than pi this time.

% bring to reflection point
left180_minhpi=left180-(.5*pi);
right180_minhpi=right180-(.5*pi);

% reflect
left_minhpi_abs=abs(left180_minhpi);
right_minhpi_abs=abs(right180_minhpi);

% de-backwards the distribution
left_minhpi_abs_neg=left_minhpi_abs*-1;
right_minhpi_abs_neg=right_minhpi_abs*-1;

% translate it back to the origin
left90=left_minhpi_abs_neg+(.5*pi);
right90=right_minhpi_abs_neg+(.5*pi);

%the values now mean "radians away from the PGG", irrespective of whether we are talking about hierarchical ascent or descent

% write tables
writetable(table(left90),'~/data/left90.csv')
writetable(table(right90),'~/data/right90.csv')

% print face values

% translate and reflect angles to yield vectors of aligned vs orthogonal (0-90 degrees)
%writetable(table((abs(abs(str2double(table2array(left(2,:))))-pi)*-1)+pi),'~/data/MeanMeanReflect90L.csv')
%writetable(table((abs(abs(str2double(table2array(right(2,:))))-pi)*-1)+pi),'~/data/MeanMeanReflect90R.csv')
% let me explain: let's first strip out the matlab nonsense (data type conversions)
%writetable(table((abs(abs(str2double(table2array(left(2,:))))-pi)*-1)+pi),'~/data/MeanMeanReflect90L.csv')
%writetable(table(((((abs(abs(abs(str2double(table2array(right(2,:))))-pi)*-1)+pi)-(.5*pi)))*-1)+(.5*pi)),'~/data/MeanMeanReflect90R.csv')

%(& ($ abs (@   (abs(abs(left(2,:))-pi)*-1)+pi - (.5*pi) @) $) * - 1 &) + (.5*pi)


% now lets work from the insider of the ()'s out:
% left(2,:) just selects the angular values. heed not the 2.
% abs(left(2,:)) translates -pi:0 to 0:pi, maintaining 0 = aligned, pi = orthogonal.
% it also translates -2pi:-pi into pi:2pi, maintaining orthogonal = low and aligned = high.
% so half our segments are halfway translated: they share appropraite vales with one OTHER segment.
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
