function lnL = getPost(error)
%     pp = [2];
%     Lp = length(pp);
%     %qq = [0,1,2,3];
       n = length(error);
%     BIC = 1e+16;
%     for j = 1:Lp
%       qq = [0];
%       for k = 1:length(qq)   
%             mod = arima(pp(j),0,qq(k));
%             try        
%                 [fit,~,logL] = estimate(mod,error,'Display','off');
%                 LOGL = logL;
%                 PQ = pp(j)+qq(k);
%                 [~,bic] = aicbic(LOGL,PQ+1,length(error));
%                 if bic < BIC
%                     BIC = bic;
%                     FIT = fit;
%                 end
%             catch ME
% %                 if (strcmp(ME.identifier),'econ:arima:arima:UnstableAR')
% %                     
% %                 end
%             end
%        end
%     end
%     error = error(:);
%     
%     if BIC < 1e+16    
%         P = FIT.P;
%         phi = ones(1,P)*nan;
%         for i = 1:P
%             phi(i) = FIT.AR{i};
%         end
% 
%         %%%%%%%%% construct matrix L and M
%         L = zeros(n-P,n-P);
%         M = zeros(n-P,P);
%         M(P+1:end,:) = 0;
%        
%         for i = 1:n-P
%            for j = 1:i   %%% matrix L
%                if j == i
%                    L(i,j) = 1;
%                else
%                    if i-j <= P
%                       L(i,j) = -phi(i-j);
%                    else
%                       L(i,j) = 0;
%                    end
%                end
%            end
%            for j = i:P
%                if i <= P
%                    M(i,j) = -phi(P-abs(i-j));
%                end
%            end
%         end
% 
%         eta = error(P+1:end);
%         eta_star = error(1:P);
% 
%         a = L*eta + M*eta_star;
%     else
%         a = error;
%         P = 1;
%     end
    a = error(:);
%     lnL = -(n)/2*(1+log((a'*a)/(n)));
     lnL = -n/2*log((a'*a)/n)-12/2*log(n);
     if isnan(lnL) || lnL == -inf
         lnL = -1e+016;
     end
%    lnL = (a'*a)^(-1/2);

end