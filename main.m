function main
    clc; 
    close all;
    clear all;
    format long;

    S = input('Enter stock price S=');
    X = input('Enter strike price X=');
    H = input('Enter high barrier H (must be smaller than S)=');
    t = input('Enter maturity t=');
    sigma = input('Enter annual volatility sigma [example: 0.25]=');
    r = input('Enter interest rate r [example: 0.15] =');
    q = input('Enter coninuous annualized dividend yield q [example: 0.05]=');
    n = input('Enter number of periods n=');
    
    %S=95;
    %X=100;
    %H=90;
    %t=1;
    %sigma=0.25;
    %r=0.15;
    %q=0.05;
    %n=192;
    
    %CLOSED FORM FORMULA (BLACK-SCHOLES)
        %Calculating c2
    x=((log((H^2)/(S*X))/log(exp(1)))+(r-q+(sigma^2)/2)*t)/(sigma*sqrt(t));
    lambda=(r-q+(sigma^2)/2)/(sigma^2);
    c2=S*exp(-q*t)*(H/S)^(2*lambda)*normcdf(x)-X*exp(-r*t)*(H/S)^(2*lambda-2)*normcdf(x-sigma*sqrt(t));
    
    %TREE CONSTRUCTION
        %Calculations
    deltaT = t/n;
    u = exp(sigma*sqrt(deltaT));
    d = 1/u;
    R = exp(r*deltaT);
    p = (exp((r-q)*deltaT)-d)/(u-d); %Taking into account dividends
    
    %Start creating the tree
    h=floor(log(H/S)/log(u));
    
    C=zeros(1,n+1);
    for i=0:n
        C(i+1)=max(0,S*(u^(n-i))*(d^i)-X);
    end
    
    if ((mod(n-h,2)==0)&((0<=((n-h)/2))&(((n-h)/2)<=n)))
        C(((n-h)/2)+1)=0;
    end
    
    for j=n-1:-1:0
        for i=0:j
            C(i+1)=(p*C(i+1)+(1-p)*C(i+2))/R;
        end
        if ((mod(j-h,2)==0)&((0<=((j-h)/2))&(((j-h)/2)<=j)))
            C(((j-h)/2)+1)=0;
        end
    end
     
    D=zeros(1,n+1);
    for i=0:n
        D(i+1)=max(0,S*(u^(n-i))*(d^i)-X);
    end
    for j=n-1:-1:0
        for i=0:j
            D(i+1)=(p*D(i+1)+(1-p)*D(i+2))/R;
        end
    end
    
    valueTree=D(1)-C(1);
    valueFormula=c2;
    relativeError=100*((valueTree-valueFormula)/valueFormula);
   
    %Showing values
    X=['The value from the tree is: ', num2str(sprintf('%0.7f',valueTree))];
    Y=['The value from the formula is: ', num2str(sprintf('%0.7f',valueFormula))];
    Z=['The value from the tree is: ', num2str(sprintf('%0.7f',relativeError)),'(%)'];
    disp('-------------------------')
    disp(X)
    disp(Y)
    disp(Z)
    disp('-------------------------')
end