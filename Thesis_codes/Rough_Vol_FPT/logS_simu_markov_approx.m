function logS = logS_simu_markov_approx(spotV, num_steps_mc, t,...
    num_path_mc, u_current, forward_curve, cs, xs, ...,
    theta, rho, lambda, nu, H, noise_mat_3d_mc, r, eps)
% IV_SIMU_MARKOV_APPROX Implicit-explicit scheme for lifted Heston model IV
% simulation under risk-neutral measure
% NOTE: U_t(i) are Markov processes, previous information is included in u_current
    %num_steps_mc = num_maturity * steps_per_maturity;
    dt = t / num_steps_mc;
    %num_factor = size(cs, 2);
    u_current_mc = ones(num_path_mc,1) * u_current;
    V_mat_mc = spotV * ones(num_path_mc, num_steps_mc+1) + [0, forward_curve]; % initialize variance paths
    logS = zeros(num_path_mc,num_steps_mc+1);
    %logS_mc_current = zeros(num_path_mc,1);
    %option_price_mat = zeros(num_maturity, num_moneyness);
    %moneyness_mat = zeros(num_maturity, num_moneyness);
    %tau_vec = zeros(num_maturity,1);

    % disp("Path simulation started>>>");
    % tic;

    for i = 1:num_steps_mc
        f_mc_current = max(V_mat_mc(:,i), eps);
        delta_u = -lambda*V_mat_mc(:,i)*dt + nu*sqrt(f_mc_current).*noise_mat_3d_mc(:,i,1);
        % u_current_mc = (u_current_mc + delta_u) ./ (1 + xs*dt);
        u_current_mc = (u_current_mc + delta_u) .* exp(-xs*dt);
        V_mat_mc(:,i+1) = V_mat_mc(:,i+1) + sum(cs .* u_current_mc,2); % DOUBLE CHECKED here
        logS(:, i+1) = logS(:, i) + (r - f_mc_current/2) * dt + ...
            sqrt(f_mc_current) .* (rho*noise_mat_3d_mc(:,i,1) + sqrt(1-rho^2)*noise_mat_3d_mc(:,i,2));
        % Option price calculation
        %{
        if mod(i, steps_per_maturity) == 0
            index_maturity = i / steps_per_maturity;
            tau_vec(index_maturity) = dt_maturity * index_maturity;
            moneyness_mat(index_maturity, :) = exp(linspace(-2*sqrt(spotV)*sqrt(tau_vec(index_maturity)),...
                2*sqrt(spotV)*sqrt(tau_vec(index_maturity)), num_moneyness));
            option_price_mat(index_maturity, :) = mean( (exp(logS_mc_current) - moneyness_mat(index_maturity, :))...
                .* (exp(logS_mc_current) - moneyness_mat(index_maturity, :) >0), 1 );
        end
        %}
    end

    % elapsed_time = toc;
    % disp(elapsed_time);
    % figure;
    % hold on;
    % for k = 1:10000
    %     plot(V_mat_mc(k,:));
    % end
    % hold off;
    
    % calculate iv surface
    %{
    iv_mat = zeros(size(option_price_mat));
    for i = 1:size(option_price_mat,1)
        for j = 1:size(option_price_mat,2)
            iv_mat(i, j) = blsimpv(1, moneyness_mat(i, j), r, tau_vec(i), option_price_mat(i, j),...
                'Class',{'call'});
        end
    end
    %}
end

