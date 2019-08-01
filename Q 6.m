t = 3;
K = 140;
h= 1/t;
r= 0.05;
u = exp(0.3*sqrt(h));
d =  exp(-0.3*sqrt(h));
discount = exp(-0.05*h);
S = 150;
p= (exp(-r*h) -d)/(u-d);
N= t+1;


for j=t:-1:1 % Total periods
    for i=1:1:j % i times go up
       for k = max(j,2*i):1:(i+j) %k -j is the range of the maximum value          
            S_t(j,i)=S*((u^(i-1)*d^(j-i)));
            
          v(j,i,k) = S *u^(k-j) - S*u^(2*i-j);
          if j == t
              FLS(j,i,k) = v(j,i,k);
          else
              FLS(j,i,k) =exp(-r*h)* (p*v(j+1,i+1,max(k+1,2*i)) + (1-p)*v(j+1,i,k+1));
          end
          
        end
    end
end

idx = find(FLS);
last = max(idx);
Premi = FLS(last);
