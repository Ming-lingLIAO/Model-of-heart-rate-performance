%%%% To Do List: (1)Th = Topt + dT
%%%%             (2)Catch the chol exceptions and by-pass such cases



! copy uni_functn_V2.m functn.m

Files = dir('./*.csv');
NumOfFiles = length(Files);
nopt = 6;	%%% NUMBER OF PARAMETERS TO BE OPTIMISED!

for Nf = 1:NumOfFiles
    DATA = readtable(Files(Nf).name);
    NumOfSamples = width(DATA)/2;
%	WinL = size(DATA,1)*0.1;
%	WinL = 10;
    Title = Files(Nf).name(1:end-4);
	
	FdlName = ['Uni_Result_', Title, '_WinL=length0.1'];
%	FdlName = ['Uni_Result_', Title, '_WinL=10'];
	if ~exist(FdlName)
	   mkdir(FdlName);
    end
	
    ResultFileName = [Title,'_result','.csv'];
	ResultFileName_raw = [Title,'_result','_raw','.csv'];
%    X = ones(NumOfSamples,nopt)*nan;
    X = ones(NumOfSamples,nopt+1)*nan;   % 2019-10-25 LGY
    RMS = ones(NumOfSamples,1);
      for j = 1:NumOfSamples
            x = ones(1,nopt)*nan;
            close all
            MAX = -1e+012;
            data = table2array(DATA(:,2*j-1:2*j));
			WinL = size(DATA(:,2*j-1:2*j),1)*0.1;
            I = find(~isnan(data(:,1)));
            data = data(I,:);
            T = data(:,1); 
            obs = smoothdata(data(:,2),'movmean',WinL);
			obs_raw = data(:,2);	
            T = T+273.15;
            [~,ia,~] = unique(T); T = T(ia); obs = obs(ia);
            
            xdata = T; ydata = obs;
            Range = max(xdata) - min(xdata);      %%% Assume: the sample dies at the very end of the experiment so the feasible temperature range has been THOROUGHLY EXPLORED
		
			% Define the lower/upper bounds of parameters of the objective function
%             bl= [log(0.1), log(0.1), log(10), log(10), log(0.1)];	% Tref = RANGEmin+dT
%             bu= [log(1000), log(Range), log(5000000), log(5000000), log(30)];

            bl= [log(0.1), 0.01, log(10), log(10), log(0.1),0];	% Tref = RANGEmin+dT
            bu= [log(1000),Range, log(5000000), log(5000000), log(30),min(ydata)];

            s = 2000;
            q = 20;
            L = (s/q)/10;
            T = 1e+06;
            cn = 2.4/sqrt(nopt);   %%% CnF?
            M = 10000;

            [bestf,bestx,NIter,D,Seq,acc,flag] = SCEM(nopt,s,q,bl,bu,xdata,ydata,L,T,cn,M);
            %%% If SEM is terminated due to chol failure, then
            %%% bestf = -1e+012;
                    
            if bestf > MAX
               MAX = bestf;
               x = bestx;
               SeqN = Seq;
               NITER = NIter;
            end

            PARAMS=[];
            offset = 200;
            for i = 1:q
                PARAMS = [PARAMS;SeqN(10000*(q-1)+offset+1:10000*(q-1)+NITER,:)];
            end
            
            x1 = exp(x(1));      
            x2 = min(xdata)+x(2); 
            x3 = exp(x(3));
            x4 = exp(x(4));
            x5 = x2+x(5);
            x6 = x(6);

            F1 = x1 * xdata/x2 .* exp(x3*(1/x2-1./xdata))./(1+exp(x4*(1/x5-1./xdata))) + x6;
            MDL = F1;         

			res = [x1,x2-273,x3,x4,x5-273,x6,flag];   % 2019-10-25 LGY
			X(j,:) = res;
			RMS(j) = rms(MDL-obs);
			
            figure(1);
            histogram(min(xdata)+exp(PARAMS(:,2))-273,20,'Normalization','probability'); hold on; 
            xlabel('Topt');
            title([Title,'Topt']);
%            FileName = ['./Results2/',Title,'_',num2str(j),'_Topt'];
            FileName = [FdlName,'/', Title,'_',num2str(j),'_Topt'];
            saveas(gcf, FileName, 'pdf');

            figure(2);
            histogram(min(xdata)+exp(PARAMS(:,2))+exp(PARAMS(:,5))-273,20,'Normalization','probability'); hold on; 
            xlabel('T_H');
            title([Title,'T_H']);
            FileName = [FdlName,'/',Title,'_',num2str(j),'_T_H'];
            saveas(gcf, FileName, 'pdf');

            figure(3);
            h1=plot(xdata-273,MDL); hold on; h2=plot(xdata-273,obs); hold on;
			plot(data(:,1),obs_raw,'r--','linewidth',1.4);
%            plot(xdata-273,F1,'k--'); hold on; plot(xdata-273,F2,'m--');
            xlabel('Temp (Celsius)');
            ylabel('HR (bp/min.)');
            legend([h1,h2],{'Model','Obs.'},'Location','southwest');           
            FileName = [FdlName,'/',Title,'_',num2str(j),'_MDL'];
%            saveas(gcf, FileName, 'pdf');
            saveas(gcf, FileName, 'png');

%			FileName = [FdlName,'/',Title,'_',num2str(j),'.mat'];
%			save (FileName)
      end


 
      RESULTS = array2table([X(:,1:5),RMS,X(:,6:7)],'VariableNames',...
	  {'Rho_ref_2',...
	          'T_ref_2', ...
              'dHa_R_2',...
              'dH_high_R_2',...
              'T_high_2',...
			  'rms',...
			  'interc','flag'});
      writetable(RESULTS,[FdlName,'/',ResultFileName]);
end


