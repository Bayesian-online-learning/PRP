 clear; close all; 
 clc;

 rng(123); 
% Generate Data
%this program is used to generate the raw sequence for the simulation in
%section 3.3 (manuscript)

%% Models
Ex = 1;
if (Ex == 1)
    % Example 3.1
    n = 1000; % total length of the sequence
    Data = zeros(n, 1);
    PP = sort(randsample(20:950, 9)); %sort(randi([20 950],1,9));
    disp(PP)
    Data(1:(PP(1)-1)) = 0.5*rand(PP(1)-1,1); %nomal, stable ;
    Data(PP(1):(PP(2)-1))=2.5+0.5*rand(PP(2)-PP(1),1); % first jump;
    Data(PP(2):(PP(3)-1)) = 1 + 0.5*rand(PP(3)-PP(2),1);
    Data(PP(3):(PP(4)-1)) = -0.5 + 0.5*rand(PP(4)-PP(3),1);
    Data(PP(4):(PP(5)-1)) = 3.5 + 0.5*rand(PP(5)-PP(4),1);
    Data(PP(5):(PP(6)-1)) = 2 + 0.5*rand(PP(6)-PP(5),1);
    Data(PP(6):(PP(7)-1)) = 0.5 + 0.5*rand(PP(7) - PP(6),1);
    Data(PP(7):(PP(8)-1)) = 2.5 + 0.5*rand(PP(8)-PP(7),1);
    Data(PP(8):(PP(9)-1)) = 0.5*rand(PP(9) - PP(8),1);
    Data(PP(9):n) = 1.5 + 0.5 *rand(n - PP(9)+1,1);
    JJ = [2.5,-1.5,-1.5,4,-1.5,-1.5,2,-2,1.5];
    x = Data';
    size_delta = 1;
elseif (Ex==2)
    % Example 3.2
    n = 1000; % total length of the sequence
    change_cutoff=10;
    Data = zeros(n, 1);
    PP = sort(randsample(50:600, 4)); %sort(randi([20 950],1,9));
    disp(PP)
    % plans: jump1 at 21 (up for 30), then down at 50,
    % jump2 at 100 (up to 20), then down at 180.
    Data(1:(PP(1)-1)) = 0.5*rand(PP(1)-1,1); %nomal, stable ;
    Data(PP(1):(PP(2)-1))=2.5+0.5*rand(PP(2)-PP(1),1); % first jump;
    Data(PP(2):(PP(3)-1)) = 0.5*rand(PP(3)-PP(2),1);
    Data(PP(3):(PP(4)-1)) = 1.25+0.5*rand(PP(4)-PP(3),1);
    Data(PP(4):n) = 0.5*rand(n-PP(4)+1,1);
    JJ = [2.5, -2.5,1.25,-1.25];
    x = Data';
    size_delta = 1;
end
%  ====================================

[sAA, BB, regret] = online_learning(x, size_delta, 0.5);

% Component A sequence
figure;  hold on;
plot(x, '*'); 
% plot(AA, 'r.', 'LineWidth', 2);  
plot(sAA, 'rx', 'LineWidth', 2);
xlabel("Stage N");
ylabel("Y");
legend('Data', '\delta');

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'Ex1.pdf') %Save figure

% Component B sequence
figure;  
plot(BB);  
ylim([0 1.5]);
xlabel("Stage N");
title('B');

% Regret sequence
figure; hold on;
ref = 2./sqrt(1:n);
plot(ref,'r-', 'LineWidth', 2); 
plot(regret, 'b-', 'LineWidth', 2);
xlabel("Stage N");
title('2/sqrt(n) and regret ');

%out = reshape(regret, 1, n);
%name = "sim1_n" + n + ".txt";
%dlmwrite(name, out, '-append');