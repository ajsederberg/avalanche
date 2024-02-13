function h_im = imagescWNans(varargin)
% sets NaN pixels to white by adding white to the colormap given    
    h_im = imagesc(varargin{:});
    %%
    cmap = colormap(gca);
    new_cmap = [1 1 1; cmap];
    colormap(gca, new_cmap)
end