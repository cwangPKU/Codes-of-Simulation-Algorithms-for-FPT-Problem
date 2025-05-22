% MAIN script for discretization scheme implementation--Lifted model
% Markovian approximation for rough Heston model
clear; clc;
rng(1);

%% Parameters
%  Rough Heston model:  (8.1.3)-(8.1.4) in (Jaber, 2018)
H = 0.1;
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
t = 1/12;
num_steps_mc = 200; % number of steps in total
num_path_mc = 3e5;    % Monte Carlo paths
dt = t / num_path_mc; 
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
forward_curve = lambda * theta * sum(cs./xs .* (1-exp(-(1:num_steps_mc)'*dt*xs)),2)'; % func g_0(t)
tic;
logS_markov_approx = logS_simu_markov_approx(V0, num_steps_mc, t,...
    num_path_mc, u_current, forward_curve, cs, xs, ...,
    theta, rho, lambda, nu, H, noise_mat_3d_mc, r, eps);
t1 = toc;

%% Hybrid scheme rough Heston
forward_curve_hyb = lambda * theta * ((1:num_steps_mc)*dt).^(H+1/2) / (H+1/2) / gamma(H+1/2);
tic;
logS_hybrid_scheme = logS_simu_hybrid_scheme(V0, num_steps_mc, t,...
    num_path_mc, forward_curve_hyb, theta, rho, lambda, nu, H, noise_mat_3d_mc, J, r, eps);
t2=toc;
%%
a = -5e-3;
b = 1e-3;
[prob_a_markov_approx, prob_b_markov_approx, prob_in_markov_approx, logS_final_markov_approx] = path_prob(a, b, logS_markov_approx);
[prob_a_hybrid_scheme, prob_b_hybrid_scheme, prob_in_hybrid_scheme, logS_final_hybrid_scheme] = path_prob(a, b, logS_hybrid_scheme);

disp([prob_a_markov_approx, prob_b_markov_approx, prob_in_markov_approx; prob_a_hybrid_scheme, prob_b_hybrid_scheme, prob_in_hybrid_scheme])

%%
W_T = logS_final_markov_approx;
prob_in = prob_in_markov_approx;
[counts, edges] = histcounts(W_T);
binWidth = edges(2) - edges(1); % 分箱宽度（假设等宽）
originalArea = sum(counts) * binWidth;
scalingFactor = prob_in / originalArea;
scaledCounts = counts * scalingFactor;

% 绘制直方图
figure;
h = bar(edges(1:end-1), scaledCounts, 1, 'EdgeColor', 'none', 'DisplayName', 'Histogram');

% 获取直方图数据
x = h.XData;                   % 柱子中心坐标 (1×N)
y = h.YData;                   % 柱子高度 (1×N)
w = diff(edges)/2;             % 柱子半宽 (1×N)

% 生成顶部边界路径（仅含柱子的左右上顶点）
boundary_x = [x - w; x + w];   % 左/右边缘 x 坐标 (2×N)
boundary_y = [y; y];           % 对应高度 (2×N)

% 将路径展开为连续点序列
boundary_x = boundary_x(:);    % 列向量: [左1,右1,左2,右2,...]
boundary_y = boundary_y(:);    % 列向量: [y1,y1,y2,y2,...]

% 绘制红色虚线（仅顶部）
hold on;
plot(boundary_x, boundary_y, '--r',...
    'LineWidth', 1.8,...
    'DisplayName', 'Conditional Probability Density');

% 图形装饰
xlabel('Position (x)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Conditional Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
title('Conditional Probability Density for Markov Approximation', 'FontSize', 14);
grid on;
box on;
legend('show');

%%
W_T = logS_final_hybrid_scheme;
prob_in = prob_in_hybrid_scheme;
[counts, edges] = histcounts(W_T);
binWidth = edges(2) - edges(1); % 分箱宽度（假设等宽）
originalArea = sum(counts) * binWidth;
scalingFactor = prob_in / originalArea;
scaledCounts = counts * scalingFactor;

figure;
h = bar(edges(1:end-1), scaledCounts, 1, 'EdgeColor', 'none', 'DisplayName', 'Histogram');

% 获取直方图数据
x = h.XData;                   % 柱子中心坐标 (1×N)
y = h.YData;                   % 柱子高度 (1×N)
w = diff(edges)/2;             % 柱子半宽 (1×N)

% 生成顶部边界路径（仅含柱子的左右上顶点）
boundary_x = [x - w; x + w];   % 左/右边缘 x 坐标 (2×N)
boundary_y = [y; y];           % 对应高度 (2×N)

% 将路径展开为连续点序列
boundary_x = boundary_x(:);    % 列向量: [左1,右1,左2,右2,...]
boundary_y = boundary_y(:);    % 列向量: [y1,y1,y2,y2,...]

% 绘制红色虚线（仅顶部）
hold on;
plot(boundary_x, boundary_y, '--r',...
    'LineWidth', 1.8,...
    'DisplayName', 'Conditional Probability Density');

% 图形装饰
xlabel('Position (x)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Conditional Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
title('Conditional Probability Density for Hybrid Scheme', 'FontSize', 14);
grid on;
box on;
legend('show');