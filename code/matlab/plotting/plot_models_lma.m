% compare the quality of different curve fits to model the temperature
% dependence of the leaf mass per area (LMA)

close all
tempRange = celsius2kelvin(linspace(10,40,100))';
R = config('R');
lmaTable = readtable(config('lmaFile'), 'ReadVariableNames',false,'FileType','text');
X = celsius2kelvin(repmat(table2array(lmaTable(1,:)),1,size(lmaTable,1)-1));
Y = reshape(table2array(lmaTable(2:end,:))',...
    1,numel(table2array(lmaTable(2:end,:))));

figure
mainAxes = axes;
lw = 1.5;
xlabel('Temperature (°C)')
ylabel('Leaf dry mass per area (g m^{-2})')

colors = lines(4);

% Linear model
lm = fitlm(X,Y);
line(kelvin2celsius(tempRange), lm.predict(tempRange),...
    'LineWidth', lw, 'Color', colors(1, :))

% Exponential model
f_exp = @(F,T) F(1)*exp(-F(2)*(298.15-T)./(R*T*298.15));
[f_fit_exp,~,~,~,MSE_exp] = nlinfit(X,Y,f_exp,[lm.predict(298.15),-40000]);
hold on
line(kelvin2celsius(tempRange),f_exp(f_fit_exp,tempRange),...
    'LineWidth', lw, 'Color', colors(2, :))

% sigmoid function
f_sig = @(F,T) F(1) ./ (1+exp(F(2)-F(3)./T)) + F(4);
[f_fit_sig,~,~,~,MSE_sig] = nlinfit(X,Y,f_sig,[max(Y) 1 1 min(Y)]);

hold on
line(kelvin2celsius(tempRange),f_sig(f_fit_sig,tempRange),...
    'LineWidth', lw, 'Color', colors(3, :))

% gaussian function
f_gauss = @(F,T) F(1)*exp(-((T-F(2))/F(3)).^2);
X = X(~isnan(Y));
Y = Y(~isnan(Y));
[gauss_fit,gof_gauss] = fit(X',Y', 'gauss1');
line(kelvin2celsius(tempRange),...
    f_gauss([gauss_fit.a1 gauss_fit.b1, gauss_fit.c1],tempRange),...
    'LineWidth', lw, 'Color', colors(4, :));

scatter(mainAxes,kelvin2celsius(X),Y,50,'o','filled',...
    'MarkerFaceColor', [.7 .7 .7],...
    'MarkerEdgeColor', [.4 .4 .4])

legend({...
    '$f_1(T) = a \cdot T + b$',...
    '$f_2(T) = a \cdot e^{-b \cdot \frac{298.15-T}{RT \cdot 298.15}}$',...
    '$f_3(T) = \frac{a}{1+e^{b-\frac{c}{T}}} + d$',...
    '$f_4(T) = a \cdot e^{{-(\frac{T-b}{c})}^2}$',...
    '$\textrm{experimental data}$'},...
    'FontSize', 12,...
    'FontName', 'Arial',...
    'Interpreter', 'latex',...
    'Box', 'off')
ylim([0 max(get(gca,'YLim'))])
set(gca,...
    'LineWidth',1.5,...
    'FontSize', 12,...
    'FontName', 'Arial')

% print RMSE values
fprintf('linear function:\tRMSE = %.3f\n', lm.RMSE)
fprintf('exponential function:\tRMSE = %.3f\n', sqrt(MSE_exp))
fprintf('sigmoid function:\tRMSE = %.3f\n', sqrt(MSE_sig))
fprintf('Gaussian function:\tRMSE = %.3f\n', gof_gauss.rmse)

print('ath-lma-models.png', '-painters', '-dpng')