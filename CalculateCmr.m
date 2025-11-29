function [Cmr] = CalculateCmr(r, m, h_uav, b) %Pr [m signals out of r are captured]
    Npts = 2^14;
    L = 40 + h_uav;
    D = flip(L * 5);
    eta = 2;
    success = 0;
    %--------------------
    for i = 1 : Npts
       radn = sqrt(rand(1, r)) * D;
       P = (L^2 + radn.^2).^-(eta/2); 
       P = sort(P,'desc');
       x = 0;
       %--------------------
       for j = 1 : r
           if (P(j) / (sum(P(1:r)) - P(j) )) > b
               x = x + 1;
           end          
       end
       %--------------------
       if(x == m)
        success = success + 1;
       end
    end
    %--------------------
    Cmr = success / Npts;
end