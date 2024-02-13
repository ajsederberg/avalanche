function y = generateOUProcess(tau, dt, T)
% y = generateOUProcess(tau, dt, T) simulates a mean-zero OU process with
% variance 1 and autocorrelation timescale tau. Simulation resolution is
% dt.

N = ceil(T/dt);
y = zeros(N, 1);

for i_t = 2:N
    y(i_t) = (1 - dt/tau)*y(i_t-1) + sqrt(2*dt/tau)*randn;
end

end