function [est_tao, est_del] = Est_Tao_Delta_2(x, N)

% N < 250, otherwise  overflow

% (2.1.9)
cumx  = cumsum(x);
barx  = cumx./(1:N);
rx    = x(N:-1:1);
barx0 = cumsum(rx)./(1:N);
bbarx = barx0(N:-1:1);
deltatao = bbarx(2:N) - barx(1:N-1); 
% Dtao2    = deltatao.^2;

% (2.1.11)
% HXT   = sum( (x - barx(N)).^2 );   % disp(HXT);
Ntao  = ( (1:N-1).*(N-1:-1:1) );
% Htao  = HXT - (Ntao).*Dtao2/N;     disp([max(Htao), min(Htao)]);

Htao2 = zeros(1, N-1);
for tt = 1: N-1
    Htao2(tt) = sum( (x(1:tt)-barx(tt)).^2) + sum( (x(tt+1:N) - bbarx(tt)).^2);
end

% disp([max(Htao2), min(Htao2)]);

% (2.2.1)
LNtao  = log(N)-log(Ntao);
LHtao  = log(Htao2)*(N-2);
LPtao  = (LNtao - LHtao)/2;

[~, est_tao] = max(LPtao);
 
% ================ 

nd = 20; % # of delta 
esd = deltatao(est_tao);

% disp('delta at est-tao'); disp(esd);

%delta_set  = linspace(esd-100, esd+100, nd);
delta_set  = linspace(esd*0.9, esd*1.1, nd);


% nd = length(delta_set);
Qdelta(nd) = 0;
for ii = 1:nd
    dt  = delta_set(ii);
    Qdelta_00 =   Htao2 + Ntao.*(dt-deltatao).^2/N;
    Qdelta(ii) = sum( Qdelta_00.^( 1-N/2 ) );
end
[~, id_delta] = max(Qdelta);
est_del = delta_set(id_delta);  
