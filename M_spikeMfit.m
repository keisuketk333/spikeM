
% $Author: Yasuko Matsubara 
% $Date: 2014-04-30

%--- input ------------------------------%
% dat:   input sequence 
% pfreq: period of the cycle (e.g., pfreq=24 hours)
%        (if you don''t need freq, set pfreq=-1;)
% outfn: output file name
% iter:  # of max iteration
% wantPlot:  (=1)if you need GUI plot, (=0) else
%----------------------------------------%


function [RSE, params] = M_spikeMfit( dat, pfreq, outfn, iter, wantPlot )
    %% fitting SpikeM
    [RSE, params] = LMFit(dat, pfreq, iter, wantPlot);
    %% plot fitting results 
    plotsRNF(dat, params, outfn);
    %% save spikeM parameters
    save([outfn,'.param'], 'params', '-ascii');
end


function [RSE, params] = LMFit( dat, pfreq, iter, wantPlot)

    %% parameter settings
    % duration of sequence
    T=length(dat);
    % init parameters
    params=init_params(dat,pfreq);
    [LBs, UBs]=init_const(dat, pfreq, T);
    % params settings
    options=setOpt(); 
    % default
    order=[];
    order=[order, 4 5]; % nc Sc
    order=[order, 2 1]; % BetaN, N_max   
    if(pfreq~=-1)
        order=[order, 9 8]; % pshift, prate
    end
    order=[order, 6]; % background noise 
    %%

    %% start fitting
    RSEP=inf;
    for i=1: iter
        %type='lin';
        %% if you want "tail-part-sensitive" fitting, please try below
        %if(i<iter/2 || i==iter); type='log'; else type='lin'; end
        if(i<iter/2); type='log'; else type='lin'; end
        params0=params;
        disp(['iter=', num2str(i), '/', num2str(iter)]);  
        %% if you want to plot
        if(wantPlot); plotsRNF(dat, params, []); end
        %% for each param
        for loc=1: length(order)
            lo=order(loc);
            if(lo==6&&strcmp(type,'lin')==1); type='log'; end
            % nb
            if(lo==4)
                [params(lo)] = FD_search(dat, params, lo, type);
            else
                try
                    [params(lo)] = lsqnonlin(...
                        @(x)F_RNF(dat, params, lo, x, type), ... 
                        params(lo), [], [], options); 
                catch
                    params=params0;
                end
                params = const(params, LBs,UBs);
            end
        end
        % compute RSError
        RSE = printRNF(dat, params,1);
        TH=0.0001;
        if( (abs(RSEP-RSE)<TH)  || (RSE < 0.8*10^-3) )
            break;
        end
        RSEP = RSE;
    if(isnan(RSEP)); params=zeros(1,9); end
    end  
    
end

function [params] = init_params(dat,Pp)
    %% init params
    % RNF-base
    N_max = sum(dat); 
    betaN = 1.0;
    slope = -1.5;
    % RNF-X
    nc=0;
    Sc=0.1;
    bgnoise=0.01; 
    % RNF-P
    Pa=0.1;
    Ps=0;
    %
    params=zeros(1,9);
    params(1) = N_max;
    params(2) = betaN;
    params(3) = slope;
    params(4) = nc;
    params(5) = Sc;
    params(6) = bgnoise;
    params(7) = Pp;
    params(8) = Pa;
    params(9) = Ps;
end

%%
function [LBs, UBs]=init_const(dat, pfreq, T)
    LB_base   = [sum(dat)       0.01         -1.5];  
    UB_base   = [inf            2.0          -1.5];  
    LB_X      = [0              0.0         0];
    UB_X      = [T/2            inf         inf]; %0.1];
    LB_P      = [pfreq          0.05        0];
    UB_P      = [pfreq          1           pfreq];
    LBs       = [LB_base,       LB_X,       LB_P];
    UBs       = [UB_base,       UB_X,       UB_P];
end

%%
function [params] = const(params, LB, UB)
    params = abs(params);
    %% pshift
    params(9) = mod(params(9), params(7));
    %% L & U bounding
    params(params < LB) = LB(params < LB);
    params(params > UB) = UB(params > UB);
end

function [options] = setOpt()
    options=optimset('MaxIter', 100, ...
    'Algorithm','levenberg-marquardt',...
    'display','off');
end

function [X] = removeSparse(X, wd, th)
    wd=ceil(wd/2);
    n=length(X);
    for t=1:n
        st=t-wd;ed=t+wd;
        if(st<1); st=1;end;
        if(ed>n); ed=n;end;
        counts=sum((X(st:ed)<th));
        len=ed-st;
        if(counts>len/2); X(t)=0; end
    end
end
% discrete fitting
function [estimate] = FD_search(dat, params, loc, scale)
    %%
    th=1.0; 
    % if starting point is too sparce, then, ignore the point
    spWD=4; dat=removeSparse(dat, spWD, th);
    loclist=find(dat>th);
    if(isempty(loclist)); st=1; 
    else st=loclist(1)-1; end
    if(st<0); st=0; end
    th=max(dat);
    loclist=find(dat==th);
    ed=loclist(1);
    %%
    idxlist=st:ed;
    sselist=zeros(1,length(idxlist));
    for i=1: length(idxlist);
        params(loc) = idxlist(i);
        sselist(i)=F_RNF(dat, params, -1, -1, scale);
    end
    minlst=find(sselist==min(sselist));
    if(isempty(minlst)); 
        estimate=1; 
    else
        estimate=idxlist(minlst(1));
    end
end
%----------------------------------%
% Rise and Fall fitting
%----------------------------------%
function sse=F_RNF(dat, params, loc, x, scale)
    if(loc~=-1)
        params(loc)=x;
    end
    T=length(dat);
    [idx,b,u]=M_spikeM(T, params);
    if(strcmp(scale,'lin')==1)
        % do nothing
    elseif(strcmp(scale,'log')==1)
        b=log(b+1);
        dat=log(dat+1);
    elseif(strcmp(scale,'R5')==1)
        b=b.^(1/5);
        dat=dat.^(1/5);   
    end
    sse = sqrt(mean((b - dat).^2));
end


%% for visualization 
function [RSE_LIN] = printRNF(dat, params, fid)
    RSE_LIN=F_RNF(dat, params, -1,-1, 'lin');
    %% output parameters
    fprintf(fid, '===================================\n');
    fprintf(fid, ['N       = ',  num2str(params(1), '%.0f'), '\n']);
    fprintf(fid, ['beta*N  = ',  num2str(params(2), '%.2f'), '\n']);
    fprintf(fid, ['slope   = ',  num2str(params(3), '%.1f'), '\n']);
    %
    fprintf(fid, ['nc      = ',  num2str(params(4), '%.0f'), '\n']);
    fprintf(fid, ['Sc      = ',  num2str(params(5), '%.2f'), '\n']);
    fprintf(fid, ['bgnoise = ',      num2str(params(6), '%.2f'), '\n']);
    %
    fprintf(fid, ['pcycle (Pp, Pa, Ps) = ', ...
                                num2str(params(7)), ' ', ...
                                num2str(params(8)), ' ', ...
                                num2str(params(9)), '\n']);
    fprintf(fid, '-----------------------------------\n');
    fprintf(fid, ['error (LIN)  = ',   num2str(RSE_LIN, '%.2e'), '\n']);        
    fprintf(fid, '===================================\n');
end

%% for visualization (LOG & LIN scale)
function [RSE] = plotsRNF(dat, params, outfn)

    T=length(dat);
    [n,b,u]=M_spikeM(T, params);
    RMSE=F_RNF(dat, params, -1,-1, 'lin');

    %% --- linear plot --- %%
    if(~isempty(outfn)); 
        figure(1); 
    else
        figure(10)
        set(gcf, 'Position', [0 0 700 250])
        subplot(121)
    end
    plot(n,dat, 'o', 'color',[.6, .6, .6]);
    hold on
    plot(n, b, 'r-');
    legend('Original','\DeltaB(n)');
    xlim([min(n),max(n)])
    xlabel('Time (n)');
    ylabel('Value (lin-lin)');
    hold off
    %
	title([ 'N =' , num2str(params(1), '%.0f'), ...
        ', \beta*N=',  num2str(params(2), '%.2f')]);    

    %% --- log plot --- %%
    if(~isempty(outfn)); 
        figure(2); 
    else
        subplot(122)
    end
    b=(b+1);
    u=(u+1);
    dat=(dat+1);
    loglog(n,dat, 'o', 'color',[.6, .6, .6])
    hold on
    loglog(n, b, 'r-'); 
    loglog(n, u, '--', 'color',[.4, .9, .4]); 
    legend('Original','\DeltaB(n)','U(n)');
    xlim([params(4), length(b)]);
    xlabel('Time (n)');
    ylabel('Value (log-log)');
    hold off
    %
    %% save figure 
    if(~isempty(outfn))
          disp(['save as: ', outfn])
          figure(1)
          saveas(gcf, [outfn,'LIN'], 'fig');
          saveas(gcf, [outfn,'LIN'], 'png');
          figure(2)
          saveas(gcf, [outfn,'LOG'], 'fig');
          saveas(gcf, [outfn,'LOG'], 'png');
    end
end
