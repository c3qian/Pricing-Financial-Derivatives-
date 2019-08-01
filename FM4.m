t = 25;
K = 140;
h= 1/t;
r= 0.05;
u = exp(0.3*sqrt(h));
d =  exp(-0.3*sqrt(h));
discount = exp(-r*h);
S = 150;
N= t+1;
% risk neural probablity
p= (exp(r*h) -d)/(u-d);
q= 1-p;


for j=N:-1:1
    for i=1:1:j
        S_t(j,i)=S*((u^(i-1)*d^(j-i)));
        if j==N 
            PEU(j,i)=max(K-S_t(j,i),0);
        else
           PEU(j,i)= discount*max((p*PEU(j+1,i+1)+q*PEU(j+1,i)),0);
           
        end
       
    end
end






