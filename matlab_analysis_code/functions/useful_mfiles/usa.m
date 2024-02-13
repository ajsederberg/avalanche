function cmap = usa(N)
% American flag colors colormap

if nargin == 0
    N = 100;
end

cmap = zeros(2*N,3);
cmap(:,1) = [linspace(0,1-1/N,N) linspace(1,0.8784,N)];
cmap(:,2) = [linspace(0.0032,1-1/N,N) linspace(1,0.0863,N)];
cmap(:,3) = [linspace(0.6471,1-1/N,N) linspace(1,0.1686,N)];