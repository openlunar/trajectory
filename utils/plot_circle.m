function plot_circle(r, col, style, width, handlevis, numpts)
%  plot_circle(r, col, style, width, handlevis, numpts)
if nargin < 6;  numpts = 100;       end
if nargin < 5;  handlevis = 'off';  end
if nargin < 4;  width = 1;          end
if nargin < 3;  style = '-';        end
if nargin < 2;  col = 'k';       end
if nargin < 1;  r = 1;              end

theta = linspace(0,2*pi, numpts); 
plot(r*cos(theta),r*sin(theta),'Color',col,'LineStyle',style,'LineWidth',width,'HandleVisibility',handlevis);
end

%Changelog
%Date           Programmer              Action
%08/23/2019     Jared T. Blanchard      File created