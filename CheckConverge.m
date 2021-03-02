function GR = CheckConverge(g,nopt,offset,q,M,Seq)
   Means = zeros(q,nopt);
   W = zeros(q,nopt);
   for i = 1:q
       s = Seq((i-1)*M+1:i*M,1:nopt);
       s = s(offset+1:g,:);
       Means(i,:) = mean(s);
       W(i,:) = var(s);
   end
   W = mean(W);
   B = var(Means)*g;
   G = sqrt((g-1)/g+(q+1)/(q*g).*B./W);
   GR = max(G);
end