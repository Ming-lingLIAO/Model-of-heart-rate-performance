  function f=functn(nopt,x,xdata,ydata)
  
    x1 = exp(x(1));      
    x2 = min(xdata)+x(2); 
    x3 = exp(x(3));
    x4 = exp(x(4));
    x5 = x2+x(5);
    x6 = x(6);

    F1 = x1 * xdata/x2 .* exp(x3*(1/x2-1./xdata))./(1+exp(x4*(1/x5-1./xdata))) + x6;
    MDL = F1;
    error = MDL(:) - ydata(:);
    f = getPost(error);
return
     
