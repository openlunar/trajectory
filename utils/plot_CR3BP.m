function plot_CR3BP(col1, col2, col_l, idx, fig_num)
%plot_CR3BP(col1, col2, fig_num)
%Plots primary and secondary body in CR3BP
%Inputs:
%   col1    =   (string)    surface color of primary (e.g 'b', 'r') or
%               (1x3)       array describing color (e.g [1,0,0] for blue)
%   col2    =   (string)    surface color of secondary (e.g 'b', 'r') or
%               (1x3)       array describing color (e.g [1,0,0] for blue)
%   col_l   =   (string)    marker color of lagrange points (e.g 'b') or
%               (1x3)       array describing color (e.g [1,0,0] for blue)
%   idx     =   (1xN)       index of Lagrange points to plot (e.g. [1,2,3]
%   fig_num =   (scalar)    number of figure where plot is desired

if nargin < 5; fig_num = get(gcf,'Number'); end
if nargin < 4; idx   = 1:5;                 end
if nargin < 3; col_l = 'k';                 end
if nargin < 2; col2 = [];        end
if nargin < 1; col1 = [];        end


hold on

lpts = CR3BPLpts();

plot_prims(col1, col2, fig_num);
plot_lpts(lpts);
view(2)
camlight
grid on
axis equal
end

%Changelog
%Date               Programmer              Action
%August 8, 2019     Jared T. Blanchard      Code written