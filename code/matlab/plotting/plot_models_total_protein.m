% compare different models to describe the temperature dependence of the
% total protein content on temperature

%% T-dependent model for protein content
clear;clc;close all
protContentTable = readtable(config('proteinContentFile'),...
    'ReadVariableNames',false,'FileType','text');
X = table2array(protContentTable(1,:));
Y_raw = table2array(protContentTable(2:end,:));
Y_raw = Y_raw/config('DW2FW')/1000;
Y = mean(Y_raw, 'omitnan');

% temperature range to plot functions
tempRange = linspace(0,40,100);

% initialize RMSE array and function definition array
rmse = nan(7,1);
functionDefinitions = {
    '$f_1(T) = aT + b$';
    '$f_2(T) = aT^2 + bT + c$';
    '$f_3(T) = aT^3 + bT^2 + cT + d$';
    '$f_4(T) = aT^4 + bT^3 + cT^2 + dT + e$';
    '$f_5(T) = \frac{a}{1 + exp(1+\frac{b}{T})} + d$';
    '$f_6(T) = a \cdot (1 + \frac{b + T}{b-c}) \cdot -\frac{T}{b}^{\frac{b}{b-c}}$';
    '$f_7(T) = \frac{1}{a^b} \cdot T^{b-c} \cdot e^{-\frac{T}{a}}$'
    '$\textrm{experimental data}$'
};

fig = figure;
hold on
lw = 1.5;

% color blind friendly palette
colors = [[51 34 136]; [17 119 51]; [136 204 238]; [221 204 119]; ...
    [204 102 119]; [170 68 153]; [136 34 85]]/255;

p = zeros(8,1);

% simple linear model
lm = fitlm(X,Y);
p(1) = plot(tempRange,lm.predict(tempRange'), 'LineWidth', lw, ...
    'Color', colors(1, :));
rmse(1) = sqrt(lm.MSE);

% polynomial model, coefficient 2
modelFunction = @(F,T) F(1).*(T.^2) + F(2).*T + F(3);
[nlm,~,~,~,MSE] =  nlinfit(X,Y,modelFunction,[1 1 1]);
p(2) = plot(tempRange,modelFunction(nlm, tempRange), 'LineWidth', lw, ...
    'Color', colors(2, :));
rmse(2) = sqrt(MSE);

% polynomial model, coefficient 3
modelFunction = @(F,T) F(1).*(T.^3) + F(2).*(T.^2) + F(3).*T + F(4);
[nlm,~,~,~,MSE] = nlinfit(X,Y,modelFunction,[1 1 1 1]);
p(3) = plot(tempRange,modelFunction(nlm,tempRange), 'LineWidth', lw, ...
    'Color', colors(3, :));
rmse(3) = sqrt(MSE);

% polynomial model, coefficient 4
modelFunction = @(F,T) F(1).*(T.^4) + F(2).*(T.^3) + F(3).*(T.^2) + F(4).*T + F(5);
[nlm,~,~,~,MSE] = nlinfit(X,Y,modelFunction,[1 1 1 1 1]);
p(4) = plot(tempRange,modelFunction(nlm,tempRange), 'LineWidth', lw, ...
    'Color', colors(4, :));
rmse(4) = sqrt(MSE);

% sigmoid function
modelFunction = @(F,T) F(1) ./ (1+exp(1+F(2)./T)) + F(3);
[nlm,~,~,~,MSE] = nlinfit(X,Y,modelFunction,[max(max(Y_raw)) 1 1 min(min(Y_raw))]);
p(5) = plot(tempRange,modelFunction(nlm,tempRange), 'LineWidth', lw, ...
    'Color', colors(5, :));
rmse(5) = sqrt(MSE);

% beta growth function
modelFunction = @(F,T) F(1)*(1+(F(2)+T)./(F(2)-F(3))).*(-T./F(2)).^(F(2)/(F(2)-F(3)));
[nlm,~,~,~,MSE] = nlinfit(X,Y,modelFunction,[max(max(Y_raw)) -10 1]);
p(6) = plot(tempRange,modelFunction(nlm,tempRange), 'LineWidth', lw, ...
    'Color', colors(6, :));
rmse(6) = sqrt(MSE);

% similar to gamma PDF
modelFunction = @(F,T) (1/F(1).^F(2)) .* T.^(F(2)-F(3)) .* exp(-T/F(1));
[nlm,~,~,~,MSE] = nlinfit(X,Y,modelFunction,[1 1 1]);
p(7) = plot(tempRange,modelFunction(nlm,tempRange), 'LineWidth', lw, ...
    'Color', colors(7, :));
rmse(7) = sqrt(MSE);

% experimental data
plot(X, Y_raw', 'ko', 'MarkerSize',8,...
    'MarkerFaceColor', [.7 .7 .7],...
    'MarkerEdgeColor', [.4 .4 .4]);
p(8) = plot(NaN, NaN, 'ko', 'MarkerSize',8,...
    'MarkerFaceColor', [.7 .7 .7],...
    'MarkerEdgeColor', [.4 .4 .4]);
legend_text = functionDefinitions;
legend(p, legend_text,...
    'LineWidth',1,...
    'Location', 'northwest',...
    'Interpreter','latex',...
    'FontSize', 12,...
    'Box', 'off',...
    'NumColumns', 2)
xlabel('Temperature (°C)')
ylabel('Total protein content (g gDW^{-1})')
ylim([0 .5])
set(gca, 'FontSize', 12, 'FontName', 'Arial')
set(fig, 'OuterPosition', [353.6667   93.6667  742.6667  605.3333])
saveas(fig,'ath-total-protein-models.png')