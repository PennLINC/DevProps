function [U, S, V] = plotcsvd(vf, nmodeplot, realTime, useComplexSVD, vectorScale,plot_name,timevo_name)
% PLOTCSVD performs a singular vector decomposition of a vector field.
%   Plots most dominant spatial modes and their time courses, and the
%   trajectory of the top 3 most dominant modes.
%
% INPUTS:
% VF, a complex velocity field., Must be a 3D or 4D vector containing
%   [x,y,time,trial].
% NMODESPLOT, the number of spatial modes to plot (default 5).
% REALTIME, vector containing the real time in seconds of all time points
%   in VF. If REALTIME is empty or zero, only spatial modes will be 
%   plotted (default 1:NT, i.e. sampling frequency 1).
% USECOMPLEXSVD, flag to decompose vectors as complex values instead of
%   separating x- and y-coordinates (default false).
% VECTORSCALE, number to increase or decrease vector lengths in plots
%   (default 1, corresponding to MATLAB default size).
%
% OUTPUTS:
% U, 
%
% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au
% adjusted by apines on Oct 22nd to save rather than display plots
% Use cat(3,vf{1},...vf{n}) to aggregate vector field across contin. segments prior to passing to this script
%% Set default values
if ~exist('nmodeplot', 'var')
    nmodeplot = 5;
end
if ~exist('realTime', 'var')
    realTime = 1:size(vf,3);
end
if ~exist('useComplexSVD', 'var')
    useComplexSVD = false;
end
if ~exist('vectorScale', 'var')
    vectorScale = 1;
end

% Check whether to plot time courses or not
if isempty(realTime) || (~isvector(realTime) && realTime == 0)
    onlySpatial = true;
else
    onlySpatial = false;
end

[nr, nc, nt, ~] = size(vf);


%% Perform SVD with complex velocity field
svdmat = reshape(vf, nr*nc, [])';
    
if ~useComplexSVD
    % Separate x and y components into separate variables for SVD
    svdmat = cat(2, real(svdmat), imag(svdmat));
end

[U,S,V] = svd(svdmat, 0);

if ~useComplexSVD
    V = V(1:nr*nc,:) - 1i*V(nr*nc+(1:nr*nc),:);
end

% Calculate total energy then crop to specified number of modes
eigVals = diag(S);
totVar = sum(eigVals.^2);
prctVar = 100 * eigVals.^2 / totVar;
U = U(:,1:nmodeplot);
S = S(1:nmodeplot, 1:nmodeplot);
V = V(:,1:nmodeplot);

% Make time modes as positive as possible
isNegative = mean(real(U), 1) < 0;
U(:, isNegative) = -U(:, isNegative);
V(:, isNegative) = -V(:, isNegative);

% Find average trajectory across all trials
reUav = zeros(nt, nmodeplot);
imUav = reUav;
absUav = imUav;
for it = 1:nt
    Uav = U(it + (0:nt:size(U,1)-nt), :);
    reUav(it, :) = mean(real(Uav),1);
    imUav(it, :) = mean(imag(Uav),1);
    absUav(it, :) = mean(abs(Uav),1);
end
if useComplexSVD
    UavLims = [min([reUav(:); imUav(:); absUav(:)]), ...
        max([reUav(:); imUav(:); absUav(:)])];
else
    UavLims = [min(reUav(:)), max(reUav(:))];
end

%% Plot spatial modes containing most energy
for imode=1:nmodeplot
    % Plot spatial mode
    if onlySpatial
        subplot(2, ceil(nmodeplot/2), imode)
    else
        subplot(2,nmodeplot,imode)
    end
    
    thisMode = reshape(V(:,imode), nr, nc);
    % convert to saveout rather than display
    quivPlot=quiver(real(thisMode), imag(thisMode), vectorScale);
    %quivPlot.Visible='off';
    set(gca,'YDir','reverse', 'XTick', [], 'YTick', []);
    axis(0.5+[0 nc 0 nr])
    title(sprintf('Mode %i, Var = %0.1f%%', imode, prctVar(imode)))
    %title(sprintf('Var = %0.1f%%, %0.2g%%', prctVar(imode), ...
    %    sum(prctVar(1:imode))))
    
    % Plot time course of spatial mode
    if ~onlySpatial
        subplot(2,nmodeplot,nmodeplot+imode)
        if useComplexSVD
            plot(realTime, reUav(:,imode), realTime, imUav(:,imode), ...
                realTime, absUav(:,imode))
        else
            timePlot=plot(realTime, reUav(:,imode));
            %timePlot.Visible='off'
        end
        ylim(UavLims)
        xlim([realTime(1), realTime(end)])
        hold on
        line(xlim, [0 0], 'Color', 'k')
        hold off
        xlabel('Time')
        ylabel('Component score')
    end
end
suptitle('Top SVD modes')
saveas(timePlot,plot_name)

%% OPTIONAL: Also plot trajectory of top 3 modes
a=figure
a=plot3(smooth(reUav(:,1)), smooth(reUav(:,2)), smooth(reUav(:,3)))
grid on
xlabel('Mode 1')
ylabel('Mode 2')
zlabel('Mode 3')
saveas(a,timevo_name)
