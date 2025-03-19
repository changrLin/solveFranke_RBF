% 生成Chebyshev点
N_values = [5, 10, 20];
for k = 1:length(N_values)
    N = N_values(k);
    x_cheb = cos((2*(1:N) - 1)*pi/(2*N));
    
    % 计算函数值作为观测数据
    y_obs = exp(2) + sin(2*pi*x_cheb);
    % 生成一组预测点
    if N == 5 
        x_pred = [-0.8,-0.2,0.1,0.5,0.7];
    else
        x_pred = linspace(0, 1, N);
    end
    % 进行Kriging插值
    sigma2 = 1;
    l = 0.5;
    K = pdist2(x_pred',x_pred');
    K_pred_obs = sigma2 * exp(-pdist2(x_pred',x_cheb').^2/(2*l^2));
    K_pred_pred = sigma2 * exp(-pdist2(x_pred',x_pred').^2/(2*l^2));
    K_obs_obs = sigma2 * exp(-pdist2(x_cheb',x_cheb').^2/(2*l^2));
    pinv_K_obs_obs = pinv(K_obs_obs);
    % 计算预测值和预测方差
    y_pred = K_pred_obs * pinv_K_obs_obs * y_obs';
    sigma2_pred = diag(K_pred_pred - K_pred_obs * pinv_K_obs_obs * K_pred_obs');
    % 计算置信区间
    z_alpha_2 = 1.96;
    lower = y_pred - z_alpha_2 * sqrt(sigma2_pred);
    upper = y_pred + z_alpha_2 * sqrt(sigma2_pred);
    % 绘制真实函数和置信区间
    figure
    hold on;
    fill([x_pred, fliplr(x_pred)], [lower', fliplr(upper')], [0.4, 0.6, 0.9]);
    plot(x_pred, exp(2) + sin(2*pi*x_pred), 'r-*', 'LineWidth', 0.5);
    legend('置信区间','真实函数值');
    title(['N = ', num2str(N)]);
    hold off;
end