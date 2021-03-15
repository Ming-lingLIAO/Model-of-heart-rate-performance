%%%% To Do List: (1)Th = Topt + dT
%%%%             (2)Catch the chol exceptions and by-pass such cases



! copy bi_functn.m bi_functn.m

Files = dir('./*.csv');
NumOfFiles = length(Files);
R1 = [0.25,0.5,0.75,0.95];
% R1 = 1;
% R2 = R1;
R = length(R1);
nopt = 11;	%%% NUMBER OF PARAMETERS TO BE OPTIMISED!
MaxMargin = 30;
MindT = 1;

Transform.Aopt1 = 0;   Transform.Aopt2 = 0;
Transform.aopt1 = 0;   Transform.aopt2 = 0;
Transform.Ah1 = 0;     Transfomr.Ah2 = 0;
Transform.ah1 = 0;     Transform.ah2 = 0;
CASE = [0];   %%% 0: uniform; 1: log-transformed; 2: exp-transformed
CnF = [0.5];
for CN = 1:length(CnF)
    cnF = CnF(CN);
    for CC = 1:length(CASE)   %%% disregard the last case, which we consider as not of great importance
    Case = CASE(CC);

    for Nf = 1:NumOfFiles   %%%%
        DATA = readtable(Files(Nf).name);
        NumOfSamples = width(DATA)/2;
%		WinL = size(DATA,1)*0.1;
%		WinL = size(DATA(:,2*j-1:2*j),1)*0.1;
%       WinL = 10;
        Title = Files(Nf).name(1:end-4);
        IDX = find(Title=='_');
        IDX = IDX(1);
        Acronym = [Title(1),upper(Title(IDX+1:IDX+2))];

%       FdlName = ['./Result_WinL=length0.2_TMV',Acronym,'_Case',num2str(Case),'_cn',num2str(cnF*10)];
		FdlName = ['Bi_Result_', Title, '_WinL=length0.1'];
%		FdlName = ['Bi_Result_', Title, '_WinL=10','_cn', num2str(cnF*10)];

        if ~exist(FdlName)
            mkdir(FdlName);
        end
        
        ResultFileName = [Title,'_result','.csv'];
        StatFileName = [Title,'_stat','.csv'];
        X = ones(NumOfSamples,nopt+2)*nan;
        RMS = ones(NumOfSamples,1);
        ACCmean = zeros(NumOfSamples,1);
        ACCmaxi = ACCmean;
        ACCmini = ACCmaxi;
        meanPAR = zeros(NumOfSamples,nopt);
        stdPAR = meanPAR;
        skewPAR = meanPAR;
        kurtPAR = meanPAR;
          for j = 1:NumOfSamples
              LDATA = size(DATA(:,2*j-1:2*j));
			  WinL = size(DATA(:,2*j-1:2*j),1)*0.1;
              if LDATA(1)>=3*WinL
                    x = ones(1,nopt)*nan;
                    close all
                    MAX = -1e+012;
                    data = table2array(DATA(:,2*j-1:2*j));
                    I = find(~isnan(data(:,1)));
                    data = data(I,:);
                    T = data(:,1); 
                    obs = smoothdata(data(:,2),'movmean',WinL);
                    obs_raw = data(:,2);
        %            obs = data(:,2);
                    T = T+273.15;
                    [~,ia,~] = unique(T); T = T(ia); obs = obs(ia);

        %			plot(T,obs)
                    %%%% copy and rename your objective function as functn.m

                    xdata = T; ydata = obs;
                    Range = max(xdata) - min(xdata);      %%%% Assume: the sample dies at the very end of the experiment so the feasible temperature range has been THOROUGHLY EXPLORED

                 for k1 = 1:R
                     r = R1(k1);
    %                 for k2 = k1+1:R-1
    %                     r2l = R1(k2);
    %                     r2u = R1(k2+1);
                        disp(['Now computing ', Acronym, ' sample ' ,num2str(j),' CnFac = ',num2str(cnF),' Case ',num2str(Case) ' r = ', num2str(r)]);
                        % Define the lower/upper bounds of parameters of the objective function
                        %%% UNQUOTE TO GET 12 OR 14 PARAMETERS MODEL
                    if Case == 0   %%% baseline of x7: min(xdata)
                        bl= [log(0.1), (r-0.25)*Range,           log(10),      log(10),      0.01,...
                            log(0.1), r*Range, log(10),      log(10),      0.01, 1e-6];
                             %% Topt2 = Topt1+dT
                        bu= [log(1000),    r*Range,      log(5000000), log(5000000), 30,...
                            log(1000),    Range,         log(5000000), log(5000000), 30, min(ydata)-1e-6];
                    elseif Case ==1   %%% base of x7: min(xdata)
                        bl= [log(0.1), log((r-0.25)*Range+0.0001),           log(10),      log(10),      log(0.01),...
                            log(0.1), log(r*Range), log(10),      log(10),      log(0.01), 1e-6];
                             %% Topt2 = Topt1+dT
                        bu= [log(1000), log(r*Range),      log(5000000), log(5000000), log(30),...
                            log(1000),  log(Range),         log(5000000), log(5000000), log(30), min(ydata)-1e-6];
                    elseif Case == 2   %%% base of x7: x2
                         bl= [log(0.1), log((r-0.25)*Range+0.0001),           log(10),      log(10),      log(0.01),...
                            log(0.1), log(0.01), log(10),      log(10),      0.01, 1e-6];
                             %% Topt2 = Topt1+dT
                         bu= [log(1000),    log(r*Range),      log(5000000), log(5000000), log(30),...
                            log(1000),    log(MaxMargin),         log(5000000), log(5000000), log(30), min(ydata)-1e-6];
                    end
                        s = 2000;
                        q = 20;
                        L = (s/q)/10;
                        T = 1e+06;        
                        cn = 2.4/sqrt(nopt)*cnF;
                        M = 10000;
                        %cn = 0.4;

                        [bestf,bestx,NIter,D,Seq,acc,flag,diffpct,diffabs,GR,orig,new,SeqExc] = SCEM(nopt,s,q,bl,bu,xdata,ydata,Case,Transform,L,T,cn,M);
                        %%%% If SEM is terminated due to chol failure, then
                        %%%% bestf = -1e+012;
                        disp(['GR = ',num2str(GR)]);
                        MAXDIFF = max(abs(diffpct(:)));     
                        MAXDIFF_abs = max(abs(diffabs(:)));
                        disp(['max(DIFF) = ',num2str(MAXDIFF),'%'])
                        disp(['max(DIFFabs) = ',num2str(MAXDIFF_abs)])
                        if bestf > MAX
                            MAX = bestf;
                            x = bestx;
                            SeqN = Seq;

                            NITER = NIter;
                            ACC = acc;

                            ACCmean(j) = mean(acc/NITER);   %%% mean/max/mini across different sequences
                            disp(['mean acc rate:',num2str(ACCmean(j))]);
                            ACCmaxi(j) = max(acc/NITER);
                            disp(['max acc rate:',num2str(ACCmaxi(j))]);
                            ACCmini(j) = min(acc/NITER);
                            disp(['min acc rate:',num2str(ACCmini(j))]);


                            DiffPct = diffpct;
                            DiffAbs0 = diffabs;
                            Padding = zeros(1,nopt);
                            Padding(1) = MAXDIFF_abs;
                            DiffAbs = [Padding;DiffAbs0];
                        end
    %                 end
                  end
                    PARAMS=[];
                    offset = 200;
                    for i = 1:q
                        PARAMS = [PARAMS;SeqN(10000*(q-1)+offset+1:10000*(q-1)+NITER,:)];
                    end

                    meanPAR(j,:) = nanmean(PARAMS(:,1:end-1));
                    stdPAR(j,:) = nanstd(PARAMS(:,1:end-1));
                    skewPAR(j,:) = skewness(PARAMS(:,1:end-1));
                    kurtPAR(j,:) = kurtosis(PARAMS(:,1:end-1));
                    if Case == 0
                       OrigPar2 = x(2);

                       OrigPar5 = x(5);

                       OrigPar7 = x(7);
                       X7_Baseline = min(xdata); 

                       OrigPar10 = x(10);
                       
                       OrigPar11 = x(11);
                    else                 
                       OrigPar2 = exp(x(2));
                       PARAMS(:,2) = exp(PARAMS(:,2));

                       OrigPar5 = exp(x(5));
                       PARAMS(:,5) = exp(PARAMS(:,5));

                       OrigPar7 = exp(x(7));
                       PARAMS(:,7) = exp(PARAMS(:,7));
                       if Case == 1
                           X7_Baseline = min(xdata);
                       else
                           X7_Baseline = min(xdata)+OrigPar2;
                       end

                       OrigPar10 = exp(x(10));
                       PARAMS(:,10) = exp(PARAMS(:,10));
                       
                       OrigPar11 = x(11);
                       PARAMS(:,11) = PARAMS(:,11);
                    end
                    x1 = exp(x(1));
                    x2 = min(xdata)+OrigPar2;
                    x3 = exp(x(3));
                    x4 = exp(x(4));
                    x5 = x2+OrigPar5;   %%% exp(x(5)): delta_T

                    x6 = exp(x(6));
                    x7 = X7_Baseline+OrigPar7;
                    x8 = exp(x(8));
                    x9 = exp(x(9));
                    x10 = x7+OrigPar10;
                    
                    x11 = OrigPar11;

                    F1 = x1 * xdata/x2 .* exp(x3*(1/x2-1./xdata))./(1+exp(x4*(1/x5-1./xdata)));
                    F2 = x6 * xdata/x7 .* exp(x8*(1/x7-1./xdata))./(1+exp(x9*(1/x10-1./xdata)));
                    MDL = F1+F2+x11;         



                    if x5 < x10
                       res = [x1,x6,x2-273,x7-273,x3,x8,x4,x9,x5-273,x10-273,x11,GR,MAXDIFF_abs];
                    else
                       res = [x6,x1,x7-273,x2-273,x8,x3,x9,x4,x10-273,x5-273,x11,GR,MAXDIFF_abs];
                    end

                    X(j,:) = res;
    %                 X_raw(j,:) = [x1,x6,x2-273,x7-273,x3,x8,x4,x9,x5-273,x10-273,flag];
                    RMS(j) = rms(MDL-obs);
                    close all;

                    figure(1);
                    histogram(min(xdata)+PARAMS(:,2)-273,20,'Normalization','probability'); hold on; 
                    xline(x2-273,'--b','linewidth',1.4); hold on;
                    histogram(X7_Baseline+PARAMS(:,7)-273,20,'Normalization','probability'); hold on;
                    xline(x7-273,'--r','linewidth',1.4);
                    xlabel('Topt');
                    title([Title,'Topt Case',num2str(Case),',cn',num2str(10*cnF)]);
                    FileName = [FdlName,'/',Title,'_',num2str(j),'Topt'];
                    saveas(gcf, FileName, 'pdf');

                    figure(2);
                    histogram(min(xdata)+PARAMS(:,2)+PARAMS(:,5)-273,20,'Normalization','probability'); hold on;
                    xline(x5-273,'--b','linewidth',1.4);
                    histogram(X7_Baseline+PARAMS(:,7)+PARAMS(:,10)-273,20,'Normalization','probability'); hold on;
                    xline(x10-273,'--r','linewidth',1.4);
                    xlabel('T_H');
                    title([Title,'TH Case',num2str(Case),',cn',num2str(10*cnF)]);
                    FileName = [FdlName,'/',Title,'_',num2str(j),'T_H'];
                    saveas(gcf, FileName, 'pdf');

                    figure(3);
                    h1=plot(xdata-273,MDL-x11); hold on; h2=plot(xdata-273,obs-x11);hold on;
                    plot(xdata-273,F1,'k--'); hold on; plot(xdata-273,F2,'m--'); hold on;
                    plot(data(:,1),obs_raw-x11,'r--','linewidth',1.4);
                    xlabel('Temp (Celsius)');
                    ylabel('HR (bp/min.)');
                    title([Title,'Rslt Case',num2str(Case),',cn',num2str(10*cnF)]);
                    legend([h1,h2],{'Model','Obs.'},'location','best');           
                    FileName = [FdlName,'/',Title,'_',num2str(j),'_MDL'];
                    saveas(gcf, FileName, 'png');
                    if max(abs(DiffPct(:)))>0
                        DIFFpctName = [Title,'_',num2str(j),'_DiffPct','.csv'];
                        writematrix(DiffPct,[FdlName,'/',DIFFpctName]);
                    end
                    if max(abs(DiffAbs(:)))>0
                        DIFFabsName = [Title,'_',num2str(j),'_DiffAbs','.csv'];
                        writematrix(DiffAbs,[FdlName,'/',DIFFabsName]);
                    end
              end
          end

%          RESULTS = array2table([X,RMS,ACCmean,ACCmaxi,ACCmini],'VariableNames',...
          RESULTS = array2table([X(:,1:10),RMS,X(:,11:13),ACCmean,ACCmaxi,ACCmini],'VariableNames',...
          {'Rho_ref_1', 'Rho_ref_2',...
                  'T_ref_1', 'T_ref_2',...
                  'dHa_R_1', 'dHa_R_2',...
                  'dH_high_R_1', 'dH_high_R_2',...
                  'T_high_1', 'T_high_2',...
				  'rms',...
				  'Bp0','GR','MaxAbsErr',...
                  'ACCmean','ACCmaxi','ACCmini'});
          writetable(RESULTS,[FdlName,'/',ResultFileName]);

          STATS = array2table([meanPAR(:,2),stdPAR(:,2),skewPAR(:,2),kurtPAR(:,2),...
                   meanPAR(:,5),stdPAR(:,5),skewPAR(:,5),kurtPAR(:,5),...
                   meanPAR(:,7),stdPAR(:,7),skewPAR(:,7),kurtPAR(:,7),...
                   meanPAR(:,10),stdPAR(:,10),skewPAR(:,10),kurtPAR(:,10)],'VariableNames',...
                   {'meanP2','stdP2','skewP2','kurtP2',...
                   'meanP5','stdP5','skewP5','kurtP5',...
                   'meanP7','stdP7','skewP7','kurtP7',...
                   'meanP10','stdP10','skewP10','kurtP10'});
          writetable(STATS,[FdlName,'/',StatFileName]);

    % 	  RESULTS_raw = array2table([X_raw,RMS],'VariableNames',...
    % 	  {'Rho_ref_1', 'Rho_ref_2',...
    % 	          'T_ref_1', 'T_ref_2',...
    %               'dHa_R_1', 'dHa_R_2',...
    %               'dH_high_R_1', 'dH_high_R_2',...
    %               'T_high_1', 'T_high_2',...
    % 			  'flag', 'rms'});
    %       writetable(RESULTS_raw,['./Results2/',ResultFileName_raw]);
    % 	  
    end
    end
end
