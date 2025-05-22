% MAIN script for discretization scheme implementation--Lifted model
% Markovian approximation for rough Heston model
clear; clc;
rng(1);

%% Parameters
%  Rough Heston model:  (8.1.3)-(8.1.4) in (Jaber, 2018)
H = 0.2;
theta = 0.05;
rho = -0.1;
lambda = 0.1;
nu = 0.3;
r = 0; d = 0;
V0 = 0.05;
eps = 0.0001;

%  Lifted Heston model: (8.2.1)-(8.2.3) in (Jaber, 2018)
nf = 100; % number of factors (50 would yield reasonable results for 1 day to maturity)
alpha = H + 0.5;
rn = 1 + 10 * nf^(-0.9);
cs = ones(1,nf) .* ((rn^(1-alpha) - 1)*rn^((alpha-1)*(1+nf/2))) ./ (gamma(alpha) * gamma(2-alpha));
xs = ones(1,nf) .* ((1-alpha)/(2-alpha)*(rn^(2-alpha)-1)/(rn^(1-alpha)-1));
for i = 1:nf
    cs(i) = cs(i) * rn^((1-alpha)*i);
    xs(i) = xs(i) * rn^(i-1-nf/2);
end

% time grid definition
dt_maturity = 4/12;  % gap between maturities in IV surface
num_maturity = 1;    % number of maturities for IV surface
num_moneyness = 41;
step_maturity = 200;  % steps in-fill (100 for reasonable results for 1/10 day to maturity)
num_steps_mc = num_maturity * step_maturity; % number of steps in total
num_path_mc = 3e5;    % Monte Carlo paths
dt = dt_maturity / step_maturity; 
J = 2; % hybrid scheme #(Gaussian noise)

%% Kernel coefficients and random matrices
Sigma = zeros(J+2,J+2);
j_vec = 2:J+1; k_vec = 2:J+1;
Sigma(1,3:J+2) = ((j_vec-1).^(H+1/2) - (j_vec-2).^(H+1/2)) / (H+1/2) * dt^(H+1/2); % dt term absorbed in path generation in reference code
F1mat = hypergeom([-H+1/2,1],H+3/2,[1,j_vec]'./[1,k_vec].*([1,j_vec]'<[1,k_vec]));
Sigma(3:J+2,3:J+2) = ((j_vec'<k_vec).*((j_vec-1)'.^(H+1/2).* (k_vec-1).^(H-1/2) .* F1mat(j_vec-1,k_vec-1) ... 
    -(j_vec-2)'.^(H+1/2).* ([1,1:J-1]).^(H-1/2) .* F1mat([1,1:J-1],[1,1:J-1]))/(H+1/2))* dt^(2*H);
Sigma = Sigma + Sigma';
Sigma(logical(eye(J+2))) = [dt, dt, ((j_vec-1).^(2*H) - (j_vec-2).^(2*H)) / (2*H) * dt^(2*H)]';

% generate noise matrix (W_i, W_i', W_{i,1},...,W_{i,J}), i = 0,1,...,nT,
% n = 1/dt, T = num_maturity * dt_maturity
noise_mat_3d_mc = normrnd(0,1,[num_path_mc,num_steps_mc,J+2]);
L = chol(Sigma)'; % lower triangular matrix; **numerical problem here if H=1/2**
% L = cholcov(Sigma);
% disp(L);

for j = J+2:-1:1
    noise_mat_3d_mc(:,:,j) = L(j,j) * noise_mat_3d_mc(:,:,j);
    for i = 1:j-1
        noise_mat_3d_mc(:,:,j) = noise_mat_3d_mc(:,:,j) + L(j,i) * noise_mat_3d_mc(:,:,i);
    end
end

%% Implicit-explicit scheme lifted Heston (Markov Approximation)
u_current = zeros(1,nf);
num_steps_mc = num_maturity * step_maturity;
forward_curve = lambda * theta * sum(cs./xs .* (1-exp(-(1:num_steps_mc)'*dt*xs)),2)'; % func g_0(t)
tic;
[iv_mat_1, tau_vec_1, moneyness_mat_1] = iv_simu_markov_approx(V0, num_maturity, dt_maturity,...
    num_moneyness, num_path_mc, step_maturity, u_current, forward_curve, cs, xs, ...,
    theta, rho, lambda, nu, H, noise_mat_3d_mc, r, eps);
t1 = toc;

%% Hybrid scheme rough Heston
forward_curve_hyb = lambda * theta * ((1:num_steps_mc)*dt).^(H+1/2) / (H+1/2) / gamma(H+1/2);
tic;
[iv_mat_2, tau_vec_2, moneyness_mat_2] = iv_simu_hybrid_scheme(V0, num_maturity, dt_maturity,...
    num_moneyness, num_path_mc, step_maturity, forward_curve_hyb,...
    theta, rho, lambda, nu, H, noise_mat_3d_mc, J, r, eps);
t2=toc;

%% Inverse transform (benchmark)
addpath('Fourier_Transform_Inversion');
n = 1000; N = 20000;
left = 1e-12; right = 200;
iv_mat_3 = zeros(size(iv_mat_2));
tic;
log_strikes = linspace(-2*sqrt(V0*dt_maturity),2*sqrt(V0*dt_maturity),num_moneyness);
%iv_mat_3(1,:) = IVCurveArray(log_strikes, H, rho, nu, lambda, theta, V0, dt_maturity, r, d, n, left, right, N);
for k = 1:num_maturity
    log_strikes = linspace(-2*sqrt(V0*dt_maturity*k),2*sqrt(V0*dt_maturity*k),num_moneyness);
    iv_mat_3(k,:) = IVCurveArray(log_strikes, H, rho, nu, lambda, theta, V0, dt_maturity, r, d, n, left, right, N);
    %disp(k);
end
t3=toc;

%% Acceleration method
u_current = zeros(1,nf);
num_steps_mc = num_maturity * step_maturity;
% forward_curve = lambda * theta * sum(cs./xs .* (1-exp(-(1:num_steps_mc)'*dt*xs)),2)'; % func g_0(t)
forward_curve_hyb = lambda * theta * ((1:num_steps_mc)*dt).^(H+1/2) / (H+1/2) / gamma(H+1/2);
tic;
[iv_mat_4, tau_vec_4, moneyness_mat_4] = iv_simu_pc(V0, num_maturity, dt_maturity,...
    num_moneyness, num_path_mc, step_maturity, u_current, forward_curve_hyb, cs, xs, ...,
    theta, rho, lambda, nu, H, noise_mat_3d_mc, J, r, eps);
t4=toc;

%% Comparison
figure;
% indices = [1,2,3,4];
indices = 1:2;
for i = 1:1
    % subplot(1, 4, i);  % Creates a 2x2 grid and selects the i-th subplot.
    hold on;
    plot(log(moneyness_mat_1(indices(i),:)), iv_mat_1(indices(i),:), "DisplayName", "Markovian Approximation Method");
    plot(log(moneyness_mat_2(indices(i),:)), iv_mat_2(indices(i),:), "DisplayName", "Hybrid Scheme");
    plot(log(moneyness_mat_2(indices(i),:)), iv_mat_3(indices(i),:), "DisplayName", "Fourier Transform Inversion");
    %plot(log(moneyness_mat_4(indices(i),:)), iv_mat_4(indices(i),:), "DisplayName", "Acceleration");
    hold off;
    
    legend('Location','best');
    title(sprintf('Implied Vol Curve for t = 4/12'));
end

%% Report
err_exact_1 = sqrt(mean((iv_mat_1 - iv_mat_3).^2,'all')); err_exact_2 = sqrt(mean((iv_mat_2 - iv_mat_3).^2,'all')); err_exact_3 = 0; err_exact_4 = sqrt(mean((iv_mat_4 - iv_mat_3).^2,'all'));
err_hyb_1 = sqrt(mean((iv_mat_1 - iv_mat_2).^2,'all')); err_hyb_2 = 0; err_hyb_3 = sqrt(mean((iv_mat_3 - iv_mat_2).^2,'all')); err_hyb_4 = sqrt(mean((iv_mat_4 - iv_mat_2).^2,'all'));
fprintf("Parameters: nf = %d, steps_per_maturity = %d.\n", nf, step_maturity);
fprintf("Markov Approximation: cpu time: %5.2f, err-to-exact: %3.2E, err-to-hybrid: %3.2E.\n", t1, err_exact_1, err_hyb_1);
fprintf("       Hybrid Scheme: cpu time: %5.2f, err-to-exact: %3.2E, err-to-hybrid: %3.2E.\n", t2, err_exact_2, err_hyb_2);
fprintf("             Inverse Transform: cpu time: %5.2f, err-to-exact: %3.2E, err-to-hybrid: %3.2E.\n", t3, err_exact_3, err_hyb_3);
fprintf("        Acceleration: cpu time: %5.2f, err-to-exact: %3.2E, err-to-hybrid: %3.2E.\n", t4, err_exact_4, err_hyb_4);
