function plot_prims(col1, col2, fig_num, origin, inert)
%plot_prims(col, fig_num)
%Plots primary and secondary body in CR3BP
%Inputs:
%   col1    =   (string)    surface color of primary (e.g 'b', 'r') or
%               (1x3)       array describing color (e.g [1,0,0] for blue)
%   col2    =   (string)    surface color of secondary (e.g 'b', 'r') or
%               (1x3)       array describing color (e.g [1,0,0] for blue)
%   fig_num =   (scalar)    number of figure where plot is desired

if nargin < 5;      inert = 0;                      end
if nargin < 4;      origin = "BARY";                end
if nargin < 3;      fig_num = get(gcf,'Number');    end
if fig_num == 0;    fig_num = get(gcf,'Number');    end 
if nargin < 2;      col2 = [];                      end
if nargin < 1;      col1 = [];                      end

ax = axis;
hold on
plot_prim(col1, fig_num, origin, inert);
plot_sec(col2, fig_num, origin, inert);
axis([min(-1.2,ax(1)),max(1.2,ax(2)),min(-1.2,ax(3)),max(1.2,ax(4))]);
axis equal %maybe don't need this
grid on
camlight
end

%Changelog
%Date               Programmer              Action
%08/08/2019         Jared T. Blanchard      file created
