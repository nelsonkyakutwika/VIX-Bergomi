

tic;
% Call the function that computes prices via quantization 
quantization_prices = quantization_calls();
total_time_quantization = toc;

tic;
% Call the function that computes prices by exact quadrature 
exact_prices = exact_calls();
total_time_exact = toc;

% Percentage relative absolute error for the puts and calls 
percentage_error = abs(exact_prices - quantization_prices)./exact_prices * 100;

% VIX call strikes, 18 strikes
strikes = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];



figure; 
plot(strikes, exact_prices, 'o', 'MarkerFaceColor', 'blue', 'MarkerSize', 4,  LineStyle='-');
hold on 
plot(strikes, quantization_prices, 'diamond', 'MarkerFaceColor', 'red', 'MarkerSize', 4,  LineStyle='--');
hold off 
legend('Exact quadrature', 'Quantization', 'Location', 'northeast', 'FontSize', 13)
legend('boxoff')
title('VIX call prices in the one-factor model, VIX future = 15.29', 'FontSize', 12);
xlabel('Strike', 'FontSize', 15)
ylabel('VIX call price', 'FontSize', 15)

% Ensure 'images' folder exists
images_folder = fullfile(pwd, 'images');
if ~exist(images_folder, 'dir')
    mkdir(images_folder);
end

% Export figure
full_path = fullfile(images_folder, 'call-prices-one-factor.pdf');
exportgraphics(gcf, full_path);

figure; 
plot(strikes, percentage_error, 'o', 'MarkerFaceColor', 'blue', 'MarkerSize', 4, LineStyle='--');
title('Relative errors in the one-factor model', 'FontSize', 13, ...
    'FontWeight', 'bold', 'Color', 'black');
xlabel('Strike', 'FontSize', 15)
ylabel('Relative error in %', 'FontSize', 15)

folder_path = 'images';  
file_name = 'RE-one-factor-calls.pdf'; 
full_path = fullfile(folder_path, file_name);
exportgraphics(gcf, full_path);


function [prices_quantization] = quantization_calls()

% Model parameters
gamma = 0.61; omega1 = 5.53; omega2 = 0.69; xi0 = 0.03; k = 1;

% VIX maturity

t = 3/12;
T1 = t;
T2 = T1 + years(days(30));

% Quantizer and probabilities
X_quantizer = load("qpoints_1000");
X_quantizer(end, :) = [];

% Extract probabilities and quantization points
Probs = X_quantizer(:, 1);
quant_points = X_quantizer(:, 2);

% Standard deviation of X
SD_X = sqrt((1 - exp(-2 * k * t)) / (2 * k));

% h(t, T)
h_tT = @(T) exp(-2 * k * (T - t)) * SD_X^2;

% x(t, T), scalar-valued, with quantization points looped inside
x_tT = @(T, X) exp(-k * (T - t)) * SD_X * X;

% Compute f(t, X)
f_tT = @(T, X) (1 - gamma) * exp(omega1 * x_tT(T,X) - 0.5 * omega1^2 * h_tT(T)) + ...
    gamma * exp(omega2 * x_tT(T,X) - 0.5 * omega2^2 * h_tT(T));

% Compute V(t, T1, T2) with integral over each quant point
V_t_T1_T2 = 10000*xi0 * arrayfun(@(X) ...
    integral(@(T) f_tT(T, X), T1, T2, 'ArrayValued', true, 'AbsTol', 1e-10, 'RelTol', 1e-6), ...
    quant_points') / (T2 - T1);

% Model future price
quant_future_price = sum(Probs' .* sqrt(V_t_T1_T2)); 
fprintf('Quantization future price is: %.4f\n', quant_future_price);

% Model call prices
% VIX call strikes 
strikes_calls = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];

call_prices_paths = sqrt(V_t_T1_T2) - strikes_calls';
call_prices_paths(call_prices_paths < 0) = 0;       % Set negatives to zero

% Model call prices
prices_quantization = sum(call_prices_paths .* repmat(Probs', size(call_prices_paths, 1), 1), 2); 

end


function [exact_call_prices] = exact_calls()

% Model parameters
gamma = 0.61; omega1 = 5.53; omega2 = 0.69; xi0 = 0.03; k = 1;

% VIX maturity
t = 3/12;
T1 = t;
T2 = T1 + years(days(30));

% Standard deviation of X
SD_X = sqrt((1 - exp(-2 * k * t)) / (2 * k));

% Standard normal distribution density
density = @(z_x) (1 / sqrt(2 * pi)) * exp(-0.5 * z_x.^2);

% h(t, T)
h_tT = @(T) exp(-2 * k * (T - t)) .* SD_X.^2;

% x(t, T)
x_tT = @(T, z_x) exp(-k * (T - t)) * SD_X * z_x;

% f(t, x)
f_t_ZX = @(T, z_x) (10000*xi0/(T2-T1)) .* ((1 - gamma) * ...
    exp(omega1 * x_tT(T, z_x) - 0.5 * omega1^2 * h_tT(T)) + ...
    gamma * exp(omega2 * x_tT(T, z_x) - 0.5 * omega2^2 * h_tT(T)));

%Integrate f(T, z_x) with respect to time over [T1, T2]
f_integrated_t = @(z_x) arrayfun(@(x) integral(@(T) f_t_ZX(T, x), T1, T2, ...
    'AbsTol', 1e-10, 'RelTol', 1e-6), z_x);

% f_integrated_t multiplied by the density
integrand_future = @(z_x) f_integrated_t(z_x).^0.5 .* density(z_x);

% Compute the single integral over z_x 
future_price = integral(integrand_future, -10, 10, 'AbsTol', 1e-10, 'RelTol', 1e-6);
fprintf('Exact future price is: %.4f\n', future_price);

% VIX call strikes
strikes_calls = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];

% Initialize call prices array
exact_call_prices = zeros(length(strikes_calls), 1);

%Loop over call strikes
for k = 1:length(strikes_calls)
    strike = strikes_calls(k);

    %Define the integrand for the current strike
    integrand_call = @(z_x) max(f_integrated_t(z_x).^0.5 - strike, 0) .* density(z_x);

    %Compute the call price for the current strike using a single integral
    exact_call_prices(k) = integral(integrand_call, -10, 10, 'AbsTol', 1e-10, 'RelTol', 1e-6);
end

end
