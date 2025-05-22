function logS = logS_simu_hybrid_scheme(spotV, num_steps_mc, t,...
    num_path_mc, forward_curve, theta, rho, lambda, nu, H, noise_mat_3d_mc, J, r, eps)
%IV_SIMU_HYBRID_SCHEME Hybrid scheme iv simulation under risk-neutral
%measure (benchmark)
    dt = t / num_steps_mc;
    V_mat_mc = spotV * ones(num_path_mc, num_steps_mc+1) + [0, forward_curve]; % initialize variance paths
    logS = zeros(num_path_mc, num_steps_mc+1);
    %logS_mc_current = zeros(num_path_mc,1); % current stock price = 1
    %option_price_mat = zeros(num_maturity, num_moneyness);
    %moneyness_mat = zeros(num_maturity, num_moneyness);
    %tau_vec = zeros(num_maturity,1);
    
    kernel_length = num_steps_mc*2;
    kernel_dt = (linspace(1,kernel_length,kernel_length).^(H+1/2) - ...
        linspace(0,kernel_length-1,kernel_length).^(H+1/2)) / (H+1/2);
    gamma_val = gamma(H+1/2);

    for i = 1:num_steps_mc
        f_mc_current = sqrt(max(V_mat_mc(:,i), eps));
        logS(:, i+1) = logS(:, i) + (r - f_mc_current.^2/2) * dt...
            + f_mc_current .* (rho * noise_mat_3d_mc(:,i,1) + sqrt(1-rho^2) * noise_mat_3d_mc(:,i,2));
        %logS(:, i+1) = logS_mc_current;
        % update volatility path
        j_vec = i+1 : num_steps_mc+1;
        adj_J = min(J, num_steps_mc+1-i);
        tmpd = nu / gamma_val * f_mc_current;
        tmpc = -dt^(H+1/2) * lambda / gamma_val * V_mat_mc(:,i); 
        V_mat_mc(:,j_vec) = V_mat_mc(:,j_vec) + tmpc * kernel_dt(j_vec-i)...
            + tmpd * dt^(H-1/2) .* [zeros(1,adj_J), kernel_dt(J+1:num_steps_mc+1-i)] .* noise_mat_3d_mc(:,i,1); 
        for j = 1:adj_J
            V_mat_mc(:,i+j) = V_mat_mc(:,i+j) + tmpd .* noise_mat_3d_mc(:,i,j+2); % dt^(H-1/2) not applicable here
        end
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
        disp(i);
    end
    %{
    % calculate iv surface
    iv_mat = zeros(size(option_price_mat));
    for i = 1:size(option_price_mat,1)
        for j = 1:size(option_price_mat,2)
            iv_mat(i, j) = blsimpv(1, moneyness_mat(i, j), r, tau_vec(i), option_price_mat(i, j),...
                'Class',{'call'});
        end
    end
    %}
end

