function [l_hat,w_hat] = decode_lsf( idx, lsfCb, l_hat_prev)


  l_hat       = lsfCb.L1(:,idx.i1+(1));
  l_hat(1:5)  = l_hat(1:5)  + lsfCb.L2(:,idx.i2+(1));
  l_hat(6:10) = l_hat(6:10) + lsfCb.L3(:,idx.i3+(1));
  
  % rearrange l_opt to guarantee min dist (twice with different values of J) 
  l_hat = rearrange(l_hat,.0012);
  l_hat = rearrange(l_hat,.0006);

  % gen best LSF vector, perform stability check
  
  L0 = lsfCb.L0{idx.i0+(1)};
  alpha = (1 - sum(L0,2));
  
  w_hat = alpha.*l_hat + sum(L0.*l_hat_prev,2);
  w_hat = stability_check(w_hat);
  
end