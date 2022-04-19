function out_arrowhandles =...
    quiver3addarrowheads(in_quivhandle,in_arrowheadlength,in_arrowtipangle)

% https://www.mathworks.com/matlabcentral/answers/354324-how-to-make-the-quiver-arrow-head-size-fixed

% Adds arrow heads with constant size to quiver3 plot.
% Arrow heads have the same inclination relative to the z plane as the vectors.
%
% Input arguments:
% in_quivhandle: Handle of quiver plot to be appended
% in_arrowheadlength: Desired arrow head length in vector length units
% in_arrowtipangle: Desired arrow head tip angle in degrees (°)
%
% Output arguments:
% out_arrowhandles: Handles to arrow head lines
X = reshape(in_quivhandle.XData,1,[]);
Y = reshape(in_quivhandle.YData,1,[]);
Z = reshape(in_quivhandle.ZData,1,[]);
U = reshape(in_quivhandle.UData,1,[]);
V = reshape(in_quivhandle.VData,1,[]);
W = reshape(in_quivhandle.WData,1,[]);
aux_Xend = X + U;
aux_Yend = Y + V;
aux_Zend = Z + W;
aux_orthvectors = cross([U;V;W],[U;V;W+1]);
aux_orthvectors = aux_orthvectors ./ vecnorm(aux_orthvectors);
if ~any(aux_orthvectors,1)
    aux_orthvectors(:,~any(aux_orthvectors,1)) = [1;0;0];
end
aux_arrowtips1 = in_arrowheadlength * (-[U;V;W] ./ vecnorm([U;V;W]) -...
    tand(in_arrowtipangle)*aux_orthvectors);
aux_arrowtips2 = aux_arrowtips1 +...
    2*in_arrowheadlength*tand(in_arrowtipangle)*aux_orthvectors;
aux_arrowhandle1 = quiver3(in_quivhandle.Parent,aux_Xend,aux_Yend,aux_Zend,...
    aux_arrowtips1(1,:),aux_arrowtips1(2,:),aux_arrowtips1(3,:),...
    'LineWidth',in_quivhandle.LineWidth,'Color',in_quivhandle.Color,...
    'ShowArrowHead','off','AutoScale','off');
aux_arrowhandle2 = quiver3(in_quivhandle.Parent,aux_Xend,aux_Yend,aux_Zend,...
    aux_arrowtips2(1,:),aux_arrowtips2(2,:),aux_arrowtips2(3,:),...
    'LineWidth',in_quivhandle.LineWidth,'Color',in_quivhandle.Color,...
    'ShowArrowHead','off','AutoScale','off');
out_arrowhandles = [aux_arrowhandle1 aux_arrowhandle2];
end
