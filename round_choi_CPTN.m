clear
n = 2; % the number of copies
d = 2; % dimension of qudits

% Assume JN is obtained from SDP

J1 = round(JN, 7);
J2 = (J1 + J1')/2;
Min = PartialTrace(J2, 2, [d_n, d]);      
Min = (Min + Min')/2;
tni_ok = IsPSD( eye(d_n) - Min );
J3 = J2;  
if IsPSD(J3)
    J3_new = J3;
else
    [V,D] = eig(J3);
    lam = diag(D);

    pos = lam > 0;
    neg = lam < 0;

    Spos = sum(lam(pos));
    nneg = sum(lam(neg));          

    eta = (trace(J3) - nneg) / (Spos - nneg);

    lam(pos) = eta .* lam(pos);
    lam(neg) = 1 - eta;

    J3_new = V*diag(lam)*V';       
    J3_new = (J3_new + J3_new')/2;
end
J4 = J3_new/trace(J3_new)*trace(JN);

