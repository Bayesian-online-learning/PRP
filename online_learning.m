function [AA, sAA, BB, regret] = online_learning(x, size_delta, rac_C)
    iL = 1;   %  L-end-point of interval 
    iR = 2;   %  R-end-point of interval
    n = length(x);
    ef = 0;   %  counter of jump found
    eft_win = zeros(n, 5);
    % rac_C   = 0.25;  % ratio for window size

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