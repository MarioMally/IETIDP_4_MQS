function mynrbplot (nurbs, subd, color)
% 
% NRBPLOT: Plot a NURBS curve or surface, or the boundary of a NURBS volume.
% 
% Calling Sequence:
% 
%   nrbplot (nrb, npnts)
%   nrbplot (nrb, npnts, p, v)
% 
% INPUT:
% 
%   nrb		: NURBS curve, surface or volume, see nrbmak.
% 
%   npnts	: Number of evaluation points, for a surface or volume, a row 
%       vector with the number of points along each direction.
% 
%   [p,v]       : property/value options
%
%               Valid property/value pairs include:
%
%               Property        Value/{Default}
%               -----------------------------------
%               light           {off} | on
%               colormap        {'copper'}
%
% Example:
%
%   Plot the test surface with 20 points along the U direction
%   and 30 along the V direction
%
%   nrbplot(nrbtestsrf, [20 30])
%
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
%    Copyright (C) 2012 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% convert the number of subdivisions in number of points
subd = subd+1;

% plot the curve or surface
if (iscell (nurbs.knots))
 if (size (nurbs.knots,2) == 2) % plot a NURBS surface
  knt = nurbs.knots;
  order = nurbs.order;
  p = nrbeval (nurbs, {linspace(knt{1}(order(1)),knt{1}(end-order(1)+1),subd(1)) ...
                       linspace(knt{2}(order(2)),knt{2}(end-order(2)+1),subd(2))});
  surf(squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)),'FaceColor',color,'LineStyle','-','FaceAlpha',1);
 elseif (size (nurbs.knots,2) == 3) % plot the boundaries of a NURBS volume
  bnd = nrbextract (nurbs);
  hold_flag = ishold;
  mynrbplot (bnd(1), subd(2:3), color);
  hold on
  mynrbplot (bnd(2), subd(2:3), color);
  mynrbplot (bnd(3), subd([1 3]), color);
  mynrbplot (bnd(4), subd([1 3]), color);
  mynrbplot (bnd(5), subd(1:2), color);
  mynrbplot (bnd(6), subd(1:2), color);
  
  if (~hold_flag)
    hold off
  end
 
 else
  error ('mynrbplot: some argument is not correct')
 end
else
  % plot a NURBS curve
  order = nurbs.order;
  p = nrbeval (nurbs, linspace (nurbs.knots(order), nurbs.knots(end-order+1), subd));

  if (any (nurbs.coefs(3,:)))
    % 3D curve
    plot3 (p(1,:), p(2,:), p(3,:)); 
    grid on;
  else
    % 2D curve
    plot (p(1,:), p(2,:));
  end
end
axis equal;

end
