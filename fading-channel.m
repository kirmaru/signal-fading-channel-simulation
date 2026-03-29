clc;
clear all;
close all;
q = 16;
T = 1e-6;
f0 = 40e6;
fs = 60e6;
t = 0:1/fs:T-1/fs;
N = length(t);
E = 1;
eps = 1.0;
delta_f = 1/T;
f_i = f0 + (0:q-1) * delta_f;
phi_cos = sqrt(2/T) * cos(2*pi*f_i' * t);
phi_sin = sqrt(2/T) * sin(2*pi*f_i' * t);
s_cos = sqrt(2*E/T) * cos(2*pi*f_i' * t);
s_sin = sqrt(2*E/T) * sin(2*pi*f_i' * t);
gamma_dB = -10:2:11;
gamma = 10.^(gamma_dB / 10);
Nerr_max = 50;
Pe_exp = zeros(size(gamma_dB));
for k = 1:length(gamma_dB)
    disp = E / (8 * gamma(k) * (1/fs));
    sigma = sqrt(disp);
    Nerr = 0;
    Ntest = 0;
    while Nerr < Nerr_max
        i = randi([1, q]);
        theta = 2 * pi * rand();
        x = sqrt(eps/2) + sqrt((1 - eps)/2) * randn();
        y = sqrt(eps/2) + sqrt((1 - eps)/2) * randn();
        nu = sqrt(x^2 + y^2);
        r = nu * (cos(theta) * s_cos(i, :) + sin(theta) * s_sin(i, :)) + 2 *
        sigma * randn(1, N);
        max_val = -inf;
        i_hat = 0;
        for j = 1:q
            rci = sum(r .* phi_cos(j, :)) * (1/fs);
            rsi = sum(r .* phi_sin(j, :)) * (1/fs);
            energy = rci^2 + rsi^2;
            if energy > max_val
                max_val = energy;
                i_hat = j;
            end
        end
        if i_hat ~= i
            Nerr = Nerr + 1;
        end
        Ntest = Ntest + 1;
    end
    Pe_exp(k) = Nerr / Ntest;
end
Pe_theory = zeros(size(gamma));
for l = 1:q-1
    C = nchoosek(q-1, l);
    coeff = (-1)^(l + 1);
    denom = 1 + l + l * (1 - eps) * gamma;
    exponent = - (l * eps * gamma) ./ denom;
    term = coeff * (1 ./ denom) .* exp(exponent);
    Pe_theory = Pe_theory + C * term;
end
figure;
semilogy(gamma_dB, Pe_exp, 'bo-', 'LineWidth', 2);
hold on;
semilogy(gamma_dB, Pe_theory, 'go-', 'LineWidth', 2);
grid on;
xlabel('Отношение SNR (дБ)');
ylabel('Вероятность ошибки');
legend('Экспериментальная вероятность', 'Теоретическая вероятность');
title(['Модуляция с ', num2str(q), ' ортогональными сигналами, \epsilon = ', num2str(eps)]);
