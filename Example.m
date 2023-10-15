 clear; close all; 
 clc;

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

iL = 1;   %  L-end-point of interval 
iR = 2;   %  R-end-point of interval

ef = 0;   %  counter of jump found
eft_win = zeros(n, 5);
rac_C   = 0.25;  % ratio for window size

AA = zeros(n, 1);
BB = AA;  Loss = AA;  CL = AA;  Rn = AA;  regret = AA;  DE = AA; 

deltaold = 0;
tic;
while iR < (n+1)
      if  (iR-iL==1)  
          difLR = x(iR) - x(iL);
          if (abs(difLR) < size_delta)
              AA(iR) = AA(iL); %% JW change
              nB  = ceil(iR*rac_C);
              t_B = AA(nB:iR); 
              dB  = bootstrp(1, @mean, t_B);
              if (dB == 0.5)
                  BB(iR) = binornd(1, 0.5);
              else          
                  BB(iR) = (dB > 0.5);
              end
              DE(iR) = DE(iL);
              iR = iR+1;  % Not Jump              
          else 
              disp('Big Jump ')
              ef = ef + 1; 
              iJ = iL;
              eft_win(ef, :) = [ef, iL, iJ, iR, difLR];
              DE(iJ:iR) = deltaold + difLR;
              deltaold = deltaold + difLR; 
              AA(iR) = ( difLR > 0 );  % 1 or 0
              nB  = ceil(iR*rac_C);
              t_B = AA(nB:iR); 
              dB  = bootstrp(1, @mean, t_B);
              if (dB == 0.5)
                  BB(iR) = binornd(1, 0.5);
              else          
                  BB(iR) = (dB > 0.5);
              end
              Loss(iL:iR) = abs( AA(iL:iR)-BB(iL:iR) ); %%% JW: iL: IR
              for ii = iL:iR 
                  CL(ii) = mean(Loss(1:ii));
                  Rn(ii) = min( mean(AA(1:ii)), 1-mean(AA(1:ii)) ); 
              end
              regret(iL:iR) = abs( CL(iL:iR)-Rn(iL:iR) ); 
              iL = iR;
              iR = iR + 1;
          end
      else 
          x_itv = x(iL : iR);  
          [est_tao, est_del] = Est_Tao_Delta(x_itv, iR-iL+1); 
          if (abs(est_del) < size_delta )
              AA(iR) = AA(iL);
              nB  = ceil(iR*rac_C);
              DE(iR) = DE(iL);
              t_B = AA(nB:iR);
              dB  = bootstrp(1, @mean, t_B); 
              if (dB == 0.5)
                  BB(iR) = binornd(1, 0.5);
              else          
                  BB(iR) = (dB > 0.5);
              end
              iR = iR + 1;   % Not-Jump 
          else 
              ef =  ef + 1;  % Jump 
              iJ = iL + est_tao;
              eft_win(ef, :) = [ef, iL, iJ, iR, est_del];
              DE(iJ:iR) = deltaold + est_del;
              deltaold = deltaold + est_del; 
              AA(iJ:iR) = (est_del > 0);
              for jj = iJ:iR 
                  nB  = ceil(jj*rac_C);
                  t_B = AA(nB:jj);
                  dB  = bootstrp(1, @mean, t_B); 
                  if (dB == 0.5)
                     BB(jj) = binornd(1, 0.5);
                  else          
                     BB(jj) = (dB > 0.5);
                  end
              end
              Loss(iL:iR) = abs(AA(iL:iR)-BB(iL:iR)); %% JW change
              for ii = iL: iR
                  CL(ii) = mean(Loss(1:ii));
                  Rn(ii) = min( mean(AA(1:ii)), 1-mean(AA(1:ii)) ); 
                  regret(ii) = abs( CL(ii)-Rn(ii) );   
              end
              % Move Forward
              iL = iJ; 
              iR = iR+1;
          end
      end
end

eft_win = eft_win(1:ef, :); 
sAA = AA;
disp(' ');
sAA(eft_win(1, 2):(eft_win(1, 4)-1)) = mean(x(eft_win(1, 2):(eft_win(1, 4)-1)));
for ii = 1:ef
    fprintf('%4d   %5d  %5d  %5d  %6.2f\n', eft_win(ii, 1:4), eft_win(ii, 5));
    if (ii < ef)
        seq = eft_win(ii+1, 2):(eft_win(ii+1, 4)-1);
    else 
        seq = eft_win(ii, 4):n;
    end 
    sAA(seq) = mean(x(eft_win(ii, 2):(eft_win(ii, 4)-1))) + eft_win(ii,5);
end
del = eft_win(:, 5);
toc;

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