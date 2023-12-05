clear;  close all;  clc;

X = readtable('data_US.csv');
x = table2array(X(209:887,'A_H1N1_pdm09'))'; % 4-6, 10-12
Date_Time = linspace(datetime(2009,1,1), datetime(2021,12,29), 679);

% X = readtable('data_UK.csv');
% x = table2array(X(731:1409,"A_H1N1_pdm09"))'; % 4-6, 10-12
% Date_Time = linspace(datetime(2009,1,1), datetime(2021,12,29), 679);

% X = readtable('data_France.csv');
% x = table2array(X(731:1409,"A_H1N1_pdm09"))'; % 4-6, 10-12
% Date_Time = linspace(datetime(2009,1,1), datetime(2021,12,29), 679);

% X = readtable('data_China.csv');
% x = table2array(X(731:1409,"A_H1N1_pdm09"))'; % 4-6, 10-12
% Date_Time = linspace(datetime(2009,1,1), datetime(2021,12,29), 679);


% yy = table2array(X(223:887,'Year'));
% ww = table2array(X(223:887, 'Week'));
% ind = not (yy > 2010 & ( ww < 41 & ww > 20));


% x = table2array(X(562:887,"B_Yamagata_lineage"))'; % 10, 14
% Date_Time = linspace(datetime(2015,10,1), datetime(2021,12,29), 326);

% ind = x > 10;
% x = x(ind);
% Date_Time = Date_Time(ind);

% x = rmmissing(x);
x(isnan(x)) = 0;
% sum(isnan(x))
n = length(x);
% plot(x);
size_delta = 5;

[sAA, BB, regret] = online_learning(x, size_delta, 0.5);
%  ====================================
min(find(BB==1))
max(find(BB==1))

figure;  hold on;
plot(Date_Time, x); 
plot(Date_Time, AA, 'b.', 'LineWidth', 2);  
plot(Date_Time, sAA, 'rx', 'LineWidth', 1);
% xlabel("n");
% legend('Data', '\delta');
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'AH1N1jump-US.pdf') %Save figure
% saveas(gcf, 'AH1N1jump-UK.pdf') %Save figure
% saveas(gcf, 'AH1N1jump-France.pdf') %Save figure
% saveas(gcf, 'AH1N1jump-China.pdf') %Save figure

% saveas(gcf, 'Byamagatajump.pdf') %Save figure

figure;  
plot(Date_Time,BB, 'bx', 'LineWidth', 2);  %title('B');
ylim([0 1.5]);
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'AH1N1B-US.pdf') %Save figure
% saveas(gcf, 'AH1N1B-UK.pdf') %Save figure
% saveas(gcf, 'AH1N1B-France.pdf') %Save figure
% saveas(gcf, 'AH1N1B-Chian.pdf') %Save figure

% saveas(gcf, 'ByamagataB.pdf') %Save figure

% prop_BB = zeros(1,354);
% for jj = 1:354 
%     nB  = ceil(jj*0.3);
%     t_B = BB(nB:jj);
%     prop_BB(jj) = mean(t_B);
% end 
% 
% sum(BB(156:end))


% figure; hold on;
% ref = 1./sqrt(1:n);
% plot(Date_Time, ref,'r-', 'LineWidth', 2); 
% plot(Date_Time, regret, 'b-', 'LineWidth', 2);
% % title('Regret ');
% set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'AH1N1regret.pdf') %Save figure
% saveas(gcf, 'Byamagataregret.pdf') %Save figure

