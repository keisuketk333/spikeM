% $Author: Yasuko Matsubara 
% $Date: 2014-04-30

%--- input ------------------------------%
% T:   duration of the sequence
% params: SpikeM parameters 
%   params(1): N
%   params(2): beta*N
%   params(3): slope
%   params(4): nb  (starting time of breaking news)
%   params(5): Sb  (strength of external shock)
%   params(6): bgn (background noise)
%   params(7): Pp  (period)
%   params(8): Pa  (strength of periodicity
%   params(9): Ps  (phase Ps of periodicity
%   params(10):B0  (default, B0=0;)
%----------------------------------------%

%--- output ------------------------------%
% idx:	index (size: Tx1)
% dB:    count of informed bloggers (size: Tx1)
% U:     count of un-informed bloggers (size: Tx1)
%----------------------------------------%

%----------------------------------%
%% give the number of awake people at time-tick n
%----------------------------------%
function [idx, dB, U] = M_spikeM(T, params)
  
    params=abs(params);
    %
    N       =  params(1);
    betaN   =  params(2);
    beta0 =  betaN/N;
    slope   = -abs(params(3));
    %
    nc      =  round(params(4));
    Sc      =  params(5);
    bgn     =  params(6);
    %
    Pp   =  params(7); 
    Pa   =  params(8); 
    Ps  =   params(9); 
    
    %--- for multi-wdsize (please ignore)
    wd=1; %7; 
    T=T*wd; nc=nc*wd;
    bgn=bgn/wd; Sc=Sc/wd;
    %--- for multi-wdsize (please ignore)

	U  =(N)*ones(1, T); 
	dB =(0)*ones(1, T); 
	B  =(0)*ones(1, T); 
   
    %% if init value (B0) is not zero
    if(length(params)==10)
        B0=params(10);
    else
        B0=0; % default setting
    end    
    
    %% init
    U(1)  = N-B0;
    dB(1) = B0;
    B(1)  = B0;

    %--------------------------------------%
    for n=0 : T-2
        dsum  = 0; % ~ number of sneezes
        for i=nc : n
            bi=dB(i+1);
            Si=exo(i, nc, Sc);
            dsum = dsum + (bi+Si) * decay_pl(n+1-i, beta0, slope);
        end

        P_n1 = period(n+1, Pp, Ps, Pa);
        dB(n+2) = P_n1 * ( U(n+1)*dsum + bgn );

        % upper-bound b(n+1)
        if(dB(n+2) > U(n+1)); dB(n+2) = U(n+1); end
  
        U(n+2) = U(n+1) - dB(n+2); 
        B(n+2) = B(n+1) + dB(n+2); 

        % avoid over fitting
        if(abs(B(n+2) + U(n+2) - N) > 0.001 || U(n+2)==0)
            dB=NaN*zeros(1, T); B=NaN*zeros(1, T); U=NaN*zeros(1, T); 
            break;
        end

    end
    %--------------------------------------%

    %--- for multi-wdsize (please ignore)
    b_tmp=zeros(1,T);
    for j=1:floor(T/wd)
        for jj=1: wd
            b_tmp(j) = b_tmp(j) + dB((j-1)*wd+jj);
        end
    end
    dB=b_tmp; T=floor(T/wd); %nc=floor(nc/wd);
    %--- for multi-wdsize (please ignore)

    % fixed duration
    idx=1:T; dB=dB(1:T); B=B(1:T); U=U(1:T);
    idx=idx'; dB=dB'; B=B'; U=U';


end
%----------------------------------%
%% power law decay, with slope (e.g., -1.5)
%----------------------------------%
function [decay] = decay_pl(n, beta, p)
    if(n<1)
       decay = 0;
    else
       decay = beta*(n^p);
    end
end

%----------------------------------%
%% exogenous shock
%----------------------------------%
function [value] = exo(t, nc, Sc)
    if(t==nc)
      value=  Sc;
    else
      value = 0;
    end
end
%----------------------------------%
%% sinusoidal periodicity
%----------------------------------%
function [value] = period(t, Pp, Ps, Pa)
    if(Pa==0 || Pp <= 0)
        value=1;
    else
        %value = 1 + Pa * cos((2*pi/Pp)*(t-Ps));
        % previous version (sin)
        value = sin((2*pi/Pp)*(t+Ps));
        value = (value + 1) * 0.5;
        value = (1 - (value*Pa));       
    end
end

