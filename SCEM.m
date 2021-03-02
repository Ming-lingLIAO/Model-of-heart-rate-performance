function [bestf,bestx,NIter,D,Seq,acc,flag] = SCEM(nopt,s,q,bl,bu,xdata,ydata,L,T,cn,M)
%%%%%%% Gaoyang Li, 2nd, Jul., 2019

%%%% nopt: dimension of the problem
%%%% s: TOTAL population size
%%%% q: number of complexes
%%%% bl: lower bound of parameters; bu: upper bound
%%%% L: Number of iters in SEM;
%%%% T: Threshold ratio for updating
%%%% cn: jump-rate

flag = 0;
m = s/q;   %%% population of each complex
NIter = 0;
% M = 10000;   %%%% Max. possible number of iterations; can be a user-specified value! (CAUTION: may induce memory overflow!)
fvalue = zeros(1,M);
GR_crit = 1.2;
GR = 10;
offset = 200;
acc = zeros(1,q);
%%%%%% STEP1: LATIN SUPERCUBE SAMPLING FROM A NON-INFORMATIVE PRIOR DISTRIBUTION
%%%%%% STEP2: EVALUATE FUNCTIONS AND SORT
%%%% D, Seq, C: the last column stores the corresponding posterior density
%%%% value
Rnd = zeros(nopt,s);
for i = 1:nopt
    rng(i+1230);
    Rnd(i,:) = randperm(s);
end
Rnd = Rnd-1;
D = zeros(s,nopt+1);
rng(1230);
Rnd2 = rand(s,nopt);
for j = 1:s
    for i = 1:nopt
%         if i ~= 5 && i ~= 10
%        if i~= 7
         D(j,i) = bl(i) + Rnd(i,j)/s*(bu(i)-bl(i)) + Rnd2(j,i)*(bu(i)-bl(i))/s;
%        else
%           rr = min(max(xdata) - (min(xdata)+exp(D(j,2))),bu(i));
%           D(j,i) = rand*rr;
%        end
    end
    %%%%%%% EVALUATE FUNCTIONS AND SORT THE POINTS IN D  !! REMEMBER TO
    %%%%%%% DEFINE functn.m!!!
    D(j,nopt+1) = functn(nopt,D(j,1:nopt),xdata,ydata);    
end
[~,I] = sort(D(:,end),'descend');
D = D(I,:);

%%%%%%% STEP3: INITIALISE q SEQUENCES
Seq = zeros(q*M,nopt+1);
Seq(1:M:(q-1)*M+1,:) = D(1:q,:);
%%%%%
C = zeros(s,nopt+1);
delta = 100; cc=1;
while  (GR > GR_crit || NIter < M/3) && NIter < M-L
    %%%%%% STEP4: PARTITION D INTO q COMPLEXES
    for i = 1:q
        IND_start = m*(i-1);
        for j = 1:m
           C(IND_start+j,:) = D(q*(j-1)+i,:);
        end
    end
    %%%%%% STEP5: EVOLVE EACH SEQUENCE AND COMPLEX
    [NIter,C,Seq,acc,flag] = SEM(NIter,q,nopt,C,Seq,bl,bu,xdata,ydata,m,M,L,T,cn,acc);  %%% NIter is the "g" in the Vrugt03 paper
    %%%%%%
    %%%%%% STEP6: UPDATE D AND SORT
    D = C;
    [~,I] = sort(D(:,end),'descend');
    D = D(I,:);
    fvalue(cc) = D(1,end); 
    if cc > offset
%         fvalue(cc-offset+1:cc)
        delta = abs((max(fvalue(cc-offset+1:cc)) - min(fvalue(cc-offset+1:cc)))/min(fvalue(cc-offset+1:cc)));
    end
    cc=cc+1;
%     D
    %%%%%
    %%%%% STEP7: COMPUTE GR AND CHECK CONVERGENCE...
    if NIter > offset+1
       GR = CheckConverge(NIter,nopt,offset,q,M,Seq);
       display(['GR = ', num2str(GR)]);
       display(['NIter = ', num2str(NIter)]);
    end
    %%%%%
    if flag == 1
       bestx = D(1,1:nopt);
       bestf = D(1,end);
       return
    end
end
   bestx = D(1,1:nopt);
   bestf = D(1,end);
end