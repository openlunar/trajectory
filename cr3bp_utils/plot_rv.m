function h = plot_rv(rv,col,proj,linewidth)
% h=plot_rv(rv,col,proj,linewidth)
% Plots a point or trajectory
% Inputs:
%     rv = [6xN] state vectors [rv1, rv2] where rv = [rx;ry;rz;vx;vy;vz]
%     col = (char) color of line ['k']
%     proj = (scalar) projection type [4]
%     linewidth = (scalar) width of plotted line [0.5]
%       1  xy
%       2  xz
%       3  yz
%       4  xyz
%       5  all

if nargin < 4;  linewidth=0.5;  end
if nargin < 3;  proj=4;         end
if nargin < 2;  col='k';        end

[n,m] = size(rv);
if n ~= 6
    if m == 6
        rv = rv';
        m = n;
    else
        error("rv must be 6xN")
    end
end

if proj <= 1
   plot(rv(1,:),rv(2,:),col,'LineWidth', linewidth); hold on
   xlabel('X [NON]');ylabel('Y [NON]')
elseif proj==2
   plot(rv(1,:),rv(3,:),col,'LineWidth', linewidth); hold on
   xlabel('X [NON]');ylabel('Z [NON]')
elseif proj==3
   plot(rv(2,:),rv(3,:),col,'LineWidth', linewidth); hold on
   xlabel('Y [NON]');ylabel('Z [NON]')
elseif proj==4
   plot3(rv(1,:),rv(2,:),rv(3,:),col,'LineWidth', linewidth); hold on
   xlabel('X [NON]');ylabel('Y [NON]');zlabel('Z [NON]')
   view(2)
elseif proj==5
   subplot(2,2,1), plot(rv(1,:),rv(2,:),col,'LineWidth', linewidth);grid;axis('equal');
   xlabel('X [NON]');ylabel('Y [NON]');
   hold on;
   subplot(2,2,2), plot(rv(1,:),rv(3,:),col,'LineWidth', linewidth);grid;axis('equal');
   xlabel('X [NON]');ylabel('Z [NON]');
   hold on;
   subplot(2,2,3), plot(rv(2,:),rv(3,:),col,'LineWidth', linewidth);grid;axis('equal');
   xlabel('Y [NON]');ylabel('Z [NON]');
   hold on;
   subplot(2,2,4), plot3(rv(1,:),rv(2,:),rv(3,:),col,'LineWidth', linewidth);
   xlabel('X [NON]');ylabel('Y [NON]');zlabel('Z [NON]');axis('equal');
   hold on;
end

grid on
axis equal

if nargout > 0; h = gcf;    end

%Changelog
%Date           Programmer              Action
%08/14/2019     Jared T. Blanchard      file created from plt.m
%08/20/2019     Jared T. Blanchard      added 'LineWidth', linewidth capability