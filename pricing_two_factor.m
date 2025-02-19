
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
plot(strikes, exact_prices, 'o', 'MarkerFaceColor', 'blue', 'MarkerSize', 4, LineStyle='-');
hold on 
plot(strikes, quantization_prices, 'diamond', 'MarkerFaceColor', 'red', 'MarkerSize', 4,  LineStyle='--');
hold off 
legend('Exact quadrature', 'Quantization', 'Location', 'northeast', 'FontSize', 13)
legend('boxoff')
title('VIX call prices in the two-factor model, VIX future = 15.40', 'FontSize', 12, ...
    'FontWeight', 'bold', 'Color', 'black');
xlabel('Strike', 'FontSize', 15)
ylabel('VIX call price', 'FontSize', 15)

% Ensure 'images' folder exists
images_folder = fullfile(pwd, 'images');
if ~exist(images_folder, 'dir')
    mkdir(images_folder);
end

% Export figure
full_path = fullfile(images_folder, 'call-prices-two-factor.pdf');
exportgraphics(gcf, full_path);


figure; 
plot(strikes, percentage_error, 'o', 'MarkerFaceColor', 'blue', 'MarkerSize', 4, LineStyle='--');
title('Relative errors in the two-factor model', 'FontSize', 13, ...
    'FontWeight', 'bold', 'Color', 'black');
xlabel('Strike', 'FontSize', 15)
ylabel('Relative error in %', 'FontSize', 15)

folder_path = 'images';  
file_name = 'RE-two-factor-calls.pdf'; 
full_path = fullfile(folder_path, file_name);
exportgraphics(gcf, full_path);


function [exact_call_prices] = exact_calls()

% Vol-of-vol parameters
k1 = 7.54; k2 = 0.24; rho = 0.7; theta = 0.23; % Set III in Bergomi 2016

gamma = 0.60; omega1 = 9.12; omega2 = 1.10; xi0 = 0.03;

% Calculation of alpha_theta^2
alphaThetaSquared = 1 / ((1-theta)^2 + theta^2 + 2*rho*theta*(1-theta));

% VIX maturity
t = 3/12;
T1 = t;
T2 = T1 + years(days(30));

% Correlation between X1 and X2
nume = rho * ((1 - exp( -(k1 + k2) * t)) / (k1 + k2));
SD_X1 = sqrt((1 - exp(-2 * k1 * t)) / (2 * k1));
SD_X2 = sqrt((1 - exp(-2 * k2 * t)) / (2 * k2));
rho12 = nume / (SD_X1 * SD_X2);


% Joint PDF for the standard bivariate normal distribution (zero mean, std = 1)
density = @(z_x1, z_x2) (1 / (2 * pi * sqrt(1 - rho12^2))) * exp(-1 / (2 * (1 - rho12^2)) * ...
    (z_x1.^2 + z_x2.^2 - 2 * rho12 * (z_x1 .* z_x2)));

% The function h(t,T)
h_tT = @(T) alphaThetaSquared * ((1-theta)^2 .* exp(-2 * k1 * (T - t)) .* SD_X1.^2 + ...
    theta^2 .* exp(-2 * k2 * (T - t)) .* SD_X2.^2 + ...
    2 * theta * (1-theta) .* exp(-(k1 + k2) * (T - t)) .* nume);

% The term x(t, T)
x_tT = @(T, z_x1, z_x2)sqrt(alphaThetaSquared) * ((1-theta) .* exp(-k1 * (T - t)) .* SD_X1 .* z_x1 + ...
    theta .* exp(-k2 * (T - t)) .* SD_X2 .* z_x2);


% The function f(t, X1, X2) multiplied by xi0/(T2-T1)
f_t_ZX1_ZX2 = @(T, z_x1, z_x2) (10000*xi0/(T2-T1)).*( (1-gamma) .* ...
    exp(omega1 .* x_tT(T, z_x1, z_x2) - 0.5 * omega1^2 .* h_tT(T)) + ...
    gamma .* exp(omega2 .* x_tT(T, z_x1, z_x2) - 0.5 * omega2^2 .* h_tT(T)) ) ;


% Integrate f(T, z_x1, z_x2) with respect to time over [T1, T2]
% Since T is scalar, we integrate for each pair of z_x1 and z_x2
% With MATLAB's default Tols
f_integrated_t = @(z_x1, z_x2) arrayfun(@(x1, x2) integral(@(T) f_t_ZX1_ZX2(T, x1, x2), T1, T2, ...
    'AbsTol', 1e-10, 'RelTol', 1e-6), z_x1, z_x2);

% f_integrated_t multiplied by the joint density
integrand_future = @(z_x1, z_x2) f_integrated_t(z_x1, z_x2).^0.5 .* density(z_x1, z_x2);

% Use integral2 to compute the double integral over z_x1 and z_x2
% To approximate infinity, we use a large range, say [-10, 10]
model_future_price = integral2(integrand_future, -10, 10, -10, 10, 'AbsTol', 1e-10, 'RelTol', 1e-6);
fprintf('Exact future price is: %.4f\n', model_future_price );

% VIX call strikes
strikes_calls = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];

% Initialize call prices array
exact_call_prices = zeros(length(strikes_calls),1);

% Loop over call strikes
for k = 1:length(strikes_calls)
    strike = strikes_calls(k);

    % Define the integrand for the current strike
    integrand_call = @(z_x1, z_x2) max(f_integrated_t(z_x1, z_x2).^0.5 - strike, 0) ...
        .* density(z_x1, z_x2);

    % Compute the call price for the current strike using integral2
    exact_call_prices(k) = integral2(integrand_call, -10, 10, -10, 10, 'AbsTol', ...
        1e-10, 'RelTol', 1e-6);
end

end


function [prices_quantization] = quantization_calls()

% Vol-of-vol parameters
k1 = 7.54; k2 = 0.24; rho = 0.7; theta = 0.23; % Set III in Bergomi 2016

gamma = 0.60; omega1 = 9.12; omega2 = 1.10; xi0 = 0.03;

% Calculation of alpha_theta^2
alphaThetaSquared = 1 / ((1-theta)^2 + theta^2 + 2*rho*theta*(1-theta));

% VIX maturity
t = 3/12;
T1 = t;
T2 = T1 + years(days(30));

% Correlation between X1 and X2
nume = rho * ((1 - exp( -(k1 + k2) * t)) / (k1 + k2));
SD_X1 = sqrt((1 - exp(-2 * k1 * t)) / (2 * k1));
SD_X2 = sqrt((1 - exp(-2 * k2 * t)) / (2 * k2));
rho12 = nume / (SD_X1 * SD_X2);

% Probabilities and quantizers
quantizer = load('qpoints_1450');
quantizer(end, :) = [];
Probs =  quantizer(:, 1);
X1_quantizer = transpose(quantizer(:, 2));
X2_quantizer_temp =  transpose(quantizer(:, 3));

%Adjust X2 quantizer
X2_quantizer = rho12.* X1_quantizer + sqrt(1 - rho12^2) .* X2_quantizer_temp;

% h(t, T)
h_tT = @(T) alphaThetaSquared * ((1-theta)^2 .* exp(-2 * k1 * (T - t)) .* SD_X1.^2 + ...
    theta^2 .* exp(-2 * k2 * (T - t)) .* SD_X2.^2 + ...
    2 * theta * (1-theta) .* exp(-(k1 + k2) * (T - t)) .* nume);

% x(t, T), scalar-valued, with quantization points looped inside
x_tT = @(T, X1, X2) sqrt(alphaThetaSquared) * ...
    ((1-theta) .* exp(-k1 * (T - t)) .* SD_X1 .* X1 + ...
    theta .* exp(-k2 * (T - t)) .* SD_X2 .* X2);

% f(t, X1, X2), scalar-valued
f_tT = @(T, X1, X2) (1 - gamma) * exp(omega1 * x_tT(T, X1, X2) - 0.5 * omega1^2 * h_tT(T)) + ...
    gamma * exp(omega2 * x_tT(T, X1, X2) - 0.5 * omega2^2 * h_tT(T));

% Compute V(t, T1, T2) with integral over each quantizer
V_t_T1_T2 = 10000*xi0 *arrayfun(@(X1, X2) ...
    integral(@(T) f_tT(T, X1, X2), T1, T2, 'ArrayValued', true, 'AbsTol', 1e-10, 'RelTol', 1e-6), ...
    X1_quantizer, X2_quantizer) / (T2 - T1);

% Model future price
model_future_price = sum(Probs' .* sqrt(V_t_T1_T2));
fprintf('Quantization future price is: %.4f\n', model_future_price );

% VIX call strikes
strikes_calls = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];

call_prices_paths = sqrt(V_t_T1_T2) - strikes_calls';
call_prices_paths(call_prices_paths < 0) = 0; % Set negatives to zero

% Model call prices
prices_quantization = sum(call_prices_paths .* repmat(Probs', size(call_prices_paths, 1), 1), 2); 

end















