function plot_prim(col, fig_num, origin, inert)
%plot_prim(col, fig_num)
%Plots primary body in CR3BP
%Inputs:
%   col     =   (string)    surface color (e.g 'b', 'r') or
%               (1x3)       array describing color (e.g [1,0,0] for blue)
%   fig_num =   (scalar)    number of figure where plot is desired

global PRIM mu RUNIT

if nargin < 4;                  inert = 0;                   end
if nargin < 3;                  origin = "BARY";             end
if nargin < 2;                  fig_num = get(gcf,'Number'); end
if fig_num == 0;                fig_num = get(gcf,'Number'); end 
if nargin < 1 || isempty(col);  col = [0.8,0.5,0.3];         end
if isfield(PRIM,'color')
    col = PRIM.color;
end

if upper(origin) == "BARY";     offset = [0;0;0];       r = mu;
elseif upper(origin) == "PRIM"  offset = [mu;0;0];      r = 0;
elseif upper(origin) == "SEC"   offset = [-1+mu;0;0];   r = 1;
end

plot_sphere(PRIM.radius/RUNIT, [-mu;0;0] + offset, fig_num, col);

if inert;   plot_circle(r, col, ':', 2, 'off'); end

grid on
axis equal
end

%Changelog
%Date               Programmer              Action
%August 8, 2019     Jared T. Blanchard      Code written