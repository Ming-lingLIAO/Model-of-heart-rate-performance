  function f=functn(nopt,x,xdata,ydata,Case,Transform)
  
    if Case == 0
       OrigPar2 = x(2);
       OrigPar5 = x(5);
       OrigPar7 = x(7);
       X7_Baseline = min(xdata); 
       OrigPar10 = x(10);
       OrigPar11 = x(11);
    else                 
       OrigPar2 = exp(x(2));
       OrigPar5 = exp(x(5));
       OrigPar7 = exp(x(7));
       if Case == 1
           X7_Baseline = min(xdata);
       else
           X7_Baseline = min(xdata)+OrigPar2;
       end
	   
       OrigPar10 = exp(x(10));
       OrigPar11 = exp(x(11));
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

    error = MDL(:) - ydata(:);
    f = getPost(error);
return
     