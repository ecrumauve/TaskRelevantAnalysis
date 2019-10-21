%%%%Correction of 1 or 0 Hit or Fa rate adapted from:

function [DP]=Dprime2(Phit,Pfa,N)

if Phit==1
    Phit=1-(1/(2*N));
elseif Phit==0
    Phit=1/(2*N);
end

if Pfa==1
    Pfa=1-(1/(2*N));
elseif Pfa==0
    Pfa=(1/(2*N));
end


ZHit=norminv(Phit,0,1);
ZFa=norminv(Pfa,0,1);
DP=abs(ZHit-ZFa);
