function [I_si_h, I_ind_i, I_multi] = computeMultiSpinsInfo(eta_vals, eps_vals, j_vals)
% computes the information that the set of spins {couplings j_vals} has about the
% latent field
% warning: number of terms computed goes at 2^length(j_vals)

if length(j_vals) > 15
    I_si_h = nan;
    I_ind_i = nan;
    I_multi = nan;
    disp('will not attempt for length(j_vals) > 15')
    return
end
I_ind_i = zeros(length(eps_vals), length(eta_vals), length(j_vals));
for ii = 1:length(j_vals)
    I_ind_i(:, :, ii) = computeSingleSpinInfo(eta_vals, eps_vals, j_vals(ii));
end
I_ind = sum(I_ind_i, 3);
% compute multi information among spins with coupling values j_vals

% this gives the rate function for a spin with coupling J_i as a function
% of field h
r_h = @(h, eta_h, eps_h, J_i) (1 + exp(-eta_h*h*J_i + eps_h)).^(-1);

I_multi = 0*I_ind;

%% construct binary state vecs
num_states = 2^length(j_vals);
s_i_include = zeros(num_states, length(j_vals));
for ii = 1:length(j_vals)
    s_i_include(:, ii) = mod(floor(((1:num_states)-1)/2^(ii-1)), 2);
end
%%
for ii = 1:length(eps_vals)
    for jj = 1:length(eta_vals)

        eta0 = eta_vals(jj);
        eps0 = eps_vals(ii);

        % compute <s_i> for each j_i
        r_bar_i = 0*j_vals;
        for kk = 1:length(j_vals)
            r_bar_i(kk) = integral(@(x) r_h(x, eta0, eps0, j_vals(kk)).*exp(-x.^2/2)/sqrt(2*pi), -5, 5);
        end

        % compute joint distribution : integral (over h) of p(h) p(s1 | h) p(s2 | h)
        P_si_joint = zeros(num_states, 1);
        % 0,0
        P_ind_s_i = ones(num_states, 1);
        for i_state = 1:num_states
            %% compose the function P(s_i = 0,1 | h) for all 0,1 combinations
            r_fun = @(h) 1;
            for kk = 1:length(j_vals)
                % construct joint probabiltiy term
                new_r_fun = @(h) (1-s_i_include(i_state, kk)) + ...
                    (2*s_i_include(i_state, kk) - 1)*r_h(h, eta0, eps0, j_vals(kk));
                r_fun = @(h) new_r_fun(h).*r_fun(h);

                % calculate independent probability term
                P_ind_s_i(i_state) = P_ind_s_i(i_state)*((1-s_i_include(i_state, kk)) + ...
                    (2*s_i_include(i_state, kk) - 1)*r_bar_i(kk));
            end
%%
            %% evaluate integral over h
            P_si_joint(i_state) = integral(@(x) normpdf(x).*r_fun(x), -5, 5);

        end
        
        %
        % compute information 
%         P_ind_s1s2 = [(1 - r_bar1)*(1 - r_bar2) (1-r_bar1)*(r_bar2); ...
%             (r_bar1)*(1-r_bar2) (r_bar1)*(r_bar2)];
       
        % 
%         r_bar1_onespin(ii, jj) = r_bar1;
%         r_bar2_onespin(ii, jj) = r_bar2;
        
%         info_r = integral(@(x) normpdf(x).*info_fun(x), -5, 5);
        I_multi(ii, jj) = sum(P_si_joint.*log(eps + P_si_joint./P_ind_s_i));
    end
end

I_si_h = I_ind - I_multi;




end