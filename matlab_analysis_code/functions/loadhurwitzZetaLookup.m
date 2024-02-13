function [glu_hurwitzZeta, glu_dhurwitzZeta] = loadhurwitzZetaLookup()

x_fn = load('avalanches_matlab_code/analysis_results/hurwitzLookup_20220118', 'spl_hurwitzZeta', ...
    'spl_dhurwitzZeta');
%%
% global glu_hurwitzZeta
glu_hurwitzZeta = @(a, x) x_fn.spl_hurwitzZeta(a, x);
% global glu_dhurwitzZeta
glu_dhurwitzZeta = @(a, x) x_fn.spl_dhurwitzZeta(a, x);
