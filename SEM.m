function [NIter,C,Seq,acc,flag] = SEM(NIter,q,nopt,C,Seq,bl,bu,xdata,ydata,m,M,L,T,cn,acc)
    beta = 0;
    flag = 0;
%     eps = 1e-12;
    while beta <= L-1
       beta = beta+1;
       NIter = NIter+1;
    %%%%% STEP I: SORT Ck, COMPUTE COVARIANCE, RETRIEVE THETA_T AND OTHER PARAMETERS;
       for k = 1:q
           Ck = C((k-1)*m+1:k*m,:);
           meanPostCk = mean(Ck(:,end));
           Sk = Seq((k-1)*M+1:k*M,:);
           mm = min(NIter,m);
           meanPostSk = mean(Sk(NIter-mm+1:NIter,end));
         
           [~,I] = sort(Ck(:,end),'descend');
           Ck = Ck(I,:);
                     
           worstPostCk = Ck(end,end);
           bestPostCk = Ck(1,end);
           muk = mean(Ck(:,1:end-1));
           Sigk = cov(Ck(:,1:end-1));

          theta_t = Sk(NIter,1:end-1);
          PostTheta = Sk(NIter,end); 
          alpha_k = exp(meanPostCk-meanPostSk);  
          try
               R = chol(cn^2*Sigk);
               if  alpha_k <= T
                   theta_tnew = theta_t  + randn(1,nopt)*R;
               else
                   theta_tnew = muk + randn(1,nopt)*R;
               end
          catch ME
%               Seq = -999;
%               return
               theta_t0 = theta_t;
               a = Ck(:,1:end-1);
               MU = zeros(1,nopt);
               STD = MU;
               for i = 1:nopt
                   MU(i) = mean(Ck(:,i));
                   STD(i) = std(Ck(:,i));
                   a(:,i) = (Ck(:,i) - MU(i))/STD(i);
                   theta_t0(i) = (theta_t0(i) - MU(i))/STD(i);
               end
              
                                    
               [coeffs,scores,latents] = pca(a);
               if isempty(latents)
                  flag = 1; 
                  return;
               end
               SumVar = sum(latents);
%                latents2 = latents;
               latents(latents/SumVar < 0.0001) = 0.0001*SumVar;  
               latents = diag(latents);
               muk = mean(scores);
               R = cn*sqrt(latents);
               
               rng(1230+k);
               if  alpha_k <= T
                   tmp = randn(1,nopt);
                   theta_tnew0 = theta_t0  + tmp*R;
                  
               else
                   theta_tnew0 = muk + randn(1,nopt)*R;
               end
               theta_tnew = (theta_tnew0 * coeffs') .* STD + MU;
               
%                a2 = coeffs*latents*coeffs';
%                DIFFpct = (a2-cov(a))./cov(a)
          end
                      
%            eig(Sigk)'
%            Sigk
       %%%%% STEP II: COMPUTE ALPHA_K AND THETA_tnew
%             alpha_k
       %%%%%
       %%%%% STEP III: CHECK VALIDITY AND COMPUTE THE POSTERIOR OF
       %%%%% THETA_TNEW
           tmp1 = theta_tnew - bl;
           tmp2 = theta_tnew - bu;
%            lhs5 = min(xdata) + exp(theta_tnew(2)) + exp(theta_tnew(5));
%             lhs7 = min(xdata)+theta_tnew(2) + theta_tnew(7);
           if min(tmp1) < 0 || max(tmp2) > 0 % || lhs7 > max(xdata) 
               PostThetaNew = -1e+12;
           else
%             theta_tnew(theta_tnew<bl) = bl(theta_tnew<bl) + (bu(theta_tnew<bl)-bl(theta_tnew<bl))/100;
%             theta_tnew(theta_tnew>bu) = bu(theta_tnew>bu) - (bu(theta_tnew>bu)-bl(theta_tnew>bu))/100;
              PostThetaNew = functn(nopt,theta_tnew,xdata,ydata);
           end
          %Omega = PostThetaNew/PostTheta;
          Omega = exp(PostThetaNew-PostTheta);
          Z = rand(1,1);
%           PostThetaNew-PostTheta
%           Ck
       %%%%%
       %%%%% STEP IV: UPDATING Ck AND Sk
          if Omega >= Z
             Ck(1,:) = [theta_tnew,PostThetaNew];   %%% replace Ck_1 with theta_tnew
             acc(k) = acc(k)+1;
          else
             theta_tnew = theta_t;
             PostThetaNew = PostTheta;
          end
%           acc
          Sk(NIter+1,:) = [theta_tnew,PostThetaNew];
       %%%%%
       %%%%% STEP V: DECIDE WHETHER THE WORST POINT IN Ck SHOULD BE UPDATED
%           Gamma_k = bestPostCk/(worstPostCk+eps);
          Gamma_k = exp(bestPostCk-worstPostCk);
          if Gamma_k > T || PostThetaNew > worstPostCk
              Ck(end,:) = [theta_tnew,PostThetaNew];
          end
       %%%%%
          C((k-1)*m+1:k*m,:) = Ck;
          Seq((k-1)*M+1:k*M,:) = Sk;
       end
    %%%%%
    end
end