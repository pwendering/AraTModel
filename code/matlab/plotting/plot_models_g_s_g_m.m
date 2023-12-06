% Compare and plot different temperature models for mesophyll conductance
% (gm) and stomatal conductance (gs)
close all

%% Data from von Caemmerer and Evans (2015), Plant, Cell and Environment
T = [15 20 25 30 35 40];
g_m = [0.13 0.175 0.223 0.225 0.195 0.129];
g_s = [0.22 0.25 0.275 0.28 0.274 0.302];

%% Mesophyll conductance
figure
subplot(1,2,1)
plot(T,g_m, '-vk')
hold on
ylim([0 0.4])
ylabel('Mesophyll conductance [mol bar^{-1} m^{-2} s^{-1}]')
xlabel('Leaf temperature [°C]')

g_m_25 = g_m(T==25);
f_buckley = @(d,T) g_m_25*exp(-d*(log(T/27.5)).^2);
[buckley_fit,~,~,~,MSE_buckley] = nlinfit(T,g_m,f_buckley, 1);
line(T,f_buckley(buckley_fit,T), 'Color', [0 .8 0], 'LineStyle', '--',...
    'LineWidth', 1.5)

f_pol_2 = @(F,T) F(1)*T.^2 + F(2)*T + F(3);
[pol_2_fit,~,~,~,MSE_pol_2] = nlinfit(T,g_m,f_pol_2, [1 1 1]);
line(T,f_pol_2(pol_2_fit,T), 'Color', [.8 0 .3], 'LineStyle', '--',...
    'LineWidth', 1.5)

f_pol_3 = @(F,T) F(1)*T.^3 + F(2)*T.^2 + F(3)*T + F(4);
[pol_3_fit,~,~,~,MSE_pol_3] = nlinfit(T,g_m,f_pol_3, [1 1 1 1]);
line(T,f_pol_3(pol_3_fit,T), 'Color', [.3 0 .8], 'LineStyle', '--',...
    'LineWidth', 1.5)

f_gauss = @(F,T) F(1)*exp(-((T-F(2))/F(3)).^2);
[gauss_fit,gof_gauss] = fit(T',g_m', 'gauss1');
MSE_gauss = gof_gauss.rmse^2;
Y_pred = f_gauss([gauss_fit.a1 gauss_fit.b1, gauss_fit.c1],T);
line(T,Y_pred, 'Color', [1 .6 .6], 'LineStyle', '--',...
    'LineWidth', 1.5)

legend({...
    'Original data (von Caemmerer & Evans, 2015)';...
    ['g_{m25} \cdot e^{-d \cdot ln(T/T_{opt})^2} (MSE = ' num2str(MSE_buckley,'%.4g') ')'];...
    ['aT^{2} + bT + c (MSE = ' num2str(MSE_pol_2,'%.4g') ')'];...
    ['aT^{3} + bT^{2} + cT + d (MSE = ' num2str(MSE_pol_3,'%.4g') ')'];...
    ['Gaussian (MSE = ' num2str(MSE_gauss,'%.4g') ')']},...
    'LineWidth', .8)
text(16,.37,'A','FontSize', 16, 'FontWeight', 'bold')

set(gca, 'LineWidth', 1.5, 'FontSize', 12)


%% Stomatal conductance
lw = 1.5;
t_range = [0 40];
% subplot(1,2,2)
figure
hold on
ylim([0 0.4])
ylabel('Stomatal conductance (mol m^{-2} s^{-1})')
xlabel('Leaf temperature (°C)')

lm = fitlm(T,g_s);
f = @(x) lm.Coefficients.Estimate(2)*x+lm.Coefficients.Estimate(1);
fplot(f, t_range, 'LineWidth', lw)

f_pol_2 = @(F,T) F(1)*T.^2 + F(2)*T + F(3);
[pol_2_fit,~,~,~,MSE_pol_2] = nlinfit(T,g_s,f_pol_2, [1 1 1]);
f = @(x)f_pol_2(pol_2_fit, x);
fplot(f, t_range, 'LineWidth', lw)

f_pol_3 = @(F,T) F(1)*T.^3 + F(2)*T.^2 + F(3)*T + F(4);
[pol_3_fit,~,~,~,MSE_pol_3] = nlinfit(T,g_s,f_pol_3, [1 1 1 1]);
f = @(x)f_pol_3(pol_3_fit, x);
fplot(f, t_range, 'LineWidth', lw)

f_gauss = @(F,T) F(1)*exp(-((T-F(2))/F(3)).^2);
[gauss_fit,gof_gauss] = fit(T',g_s', 'gauss1');
MSE_gauss = gof_gauss.rmse^2;
f = @(x)f_gauss([gauss_fit.a1 gauss_fit.b1, gauss_fit.c1], x);
fplot(f, t_range, 'LineWidth', lw)

plot(T, g_s, 'o', 'MarkerSize', 8,...
    'MarkerFaceColor', [.7 .7 .7],...
    'MarkerEdgeColor', [.4 .4 .4])

legend({...
    ['aT + b (MSE = ' num2str(lm.MSE, '%.2e') ')'];...
    ['aT^{2} + bT + c (MSE = ' num2str(MSE_pol_2,'%.2e') ')'];...
    ['aT^{3} + bT^{2} + cT + d (MSE = ' num2str(MSE_pol_3,'%.2e') ')'];...
    ['Gaussian (MSE = ' num2str(MSE_gauss,'%.2e') ')'];...
    'von Caemmerer & Evans (2015)'},...
    'Location', 'southeast', 'LineWidth', .8,...
    'box', 'off')

set(gca, 'LineWidth', 1.3, 'FontSize', 12, 'box', 'off')

set(gcf, 'OuterPosition', 1000*[-1.1977    0.1230    0.5187    0.4740])

saveas(gcf, 'ath-conductances-model-fits-Caemmerer-Evans-2015.png')