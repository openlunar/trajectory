function plot_lpts(lpts, idx, col)

if nargin < 3; col = "k";
if nargin < 2; idx = 1:5; end

[n, m] = size(lpts);
if n ~= 3 && m == 3 
    lpts = lpts';
end

for i = idx
    plot3(lpts(1,i), lpts(2,i), lpts(3,i), col+"x", "HandleVisibility","off")
end
end