% 定义参数范围 (改为5个点)
sigma_min = 0.1;
sigma_max = 2;
delta_min = -2;
delta_max = -0.1;
n_points = 5;

% 创建sigma和delta的网格
sigma = linspace(sigma_min, sigma_max, n_points);
delta = linspace(delta_min, delta_max, n_points);
[sigma_grid, delta_grid] = meshgrid(sigma, delta);

% 定义税率t的值
t_values = [0.10, 0.5, 1];

% 预分配存储结构
results = struct();

for i = 1:length(t_values)
    t = t_values(i);
    
    % 避免除以零和负底数
    epsilon = 1e-6;
    denominator = sigma_grid - delta_grid;
    denominator(denominator == 0) = epsilon;
    
    % 计算指数
    exponent_A = (sigma_grid .* (1 + delta_grid)) ./ denominator;
    exponent_B = (delta_grid .* (1 + sigma_grid)) ./ denominator;
    exponent_T = ((1+sigma_grid) .*delta_grid)./ denominator;
    
    % 计算各部分结果
    term_CS = (1 + t).^(exponent_A ) - 1;
    term_CS = term_CS ./ (delta_grid + 1 + epsilon);
    
    term_PS = 1 - ( (1 + t).^(exponent_B) );
    term_PS = term_PS ./ (sigma_grid + 1 + epsilon);
    
    tax_revenue = t .* ( (1 + t).^exponent_T );
    DWL = term_CS + term_PS - tax_revenue;
    
    % 存储结果
    results(i).t = t;
    results(i).sigma = sigma;
    results(i).delta = delta;
    results(i).term_CS = term_CS;
    results(i).term_PS = term_PS;
    results(i).tax_revenue = tax_revenue;
    results(i).DWL = DWL;
end

% 格式化输出结果
for i = 1:length(results)
    fprintf('==================== t = %.2f ====================\n', results(i).t);
    
    % 输出sigma值
    fprintf('sigma = [');
    fprintf('%.4f ', results(i).sigma);
    fprintf(']\n\n');
    
    % 输出delta值
    fprintf('delta = [');
    fprintf('%.4f ', results(i).delta);
    fprintf(']\n\n');
    
    % 输出消费者剩余变化
    disp('Δ Consumer Surplus:');
    disp(array2table(results(i).term_CS,...
        'RowNames', arrayfun(@(x)sprintf('δ=%.3f',x),results(i).delta,'UniformOutput',false),...
        'VariableNames', arrayfun(@(x)sprintf('σ=%.3f',x),results(i).sigma,'UniformOutput',false)));
    
    % 输出生产者剩余变化
    disp('Δ Producer Surplus:');
    disp(array2table(results(i).term_PS,...
        'RowNames', arrayfun(@(x)sprintf('δ=%.3f',x),results(i).delta,'UniformOutput',false),...
        'VariableNames', arrayfun(@(x)sprintf('σ=%.3f',x),results(i).sigma,'UniformOutput',false)));
    
    % 输出税收收入
    disp('Tax Revenue:');
    disp(array2table(results(i).tax_revenue,...
        'RowNames', arrayfun(@(x)sprintf('δ=%.3f',x),results(i).delta,'UniformOutput',false),...
        'VariableNames', arrayfun(@(x)sprintf('σ=%.3f',x),results(i).sigma,'UniformOutput',false)));
    
    % 输出无谓损失
    disp('Deadweight Loss:');
    disp(array2table(results(i).DWL,...
        'RowNames', arrayfun(@(x)sprintf('δ=%.3f',x),results(i).delta,'UniformOutput',false),...
        'VariableNames', arrayfun(@(x)sprintf('σ=%.3f',x),results(i).sigma,'UniformOutput',false)));
    
    fprintf('\n\n');
end