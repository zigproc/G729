function [l_hat,w_hat,idx] = quantize_lsf(w,l_hat_prev,lsfCb)

  % Section 3.2.4
  assert(length(w) == 10);

  % compute weights based on unquantized LSF 
  weights = adapt_weights(w);
  
  % compute prediction based on MA
  [E0,l1_0,l2_0,l3_0] = lsf_vq(w,...
                               weights,
                               l_hat_prev,...
                               lsfCb.L0{0+(1)},...
                               lsfCb.L1,...
                               lsfCb.L2,...
                               lsfCb.L3);
                               
  [E1,l1_1,l2_1,l3_1] = lsf_vq(w,...
                               weights,
                               l_hat_prev,...
                               lsfCb.L0{1+(1)},...
                               lsfCb.L1,...
                               lsfCb.L2,...
                               lsfCb.L3);
 
  if( E0 < E1)
    L0_idx = 1;
    L1_idx = l1_0;
    L2_idx = l2_0;
    L3_idx = l3_0;
  else
    L0_idx = 2;
    L1_idx = l1_1;
    L2_idx = l2_1;
    L3_idx = l3_1;  
  end
    
  % convert to indices (for transmission)  
  idx.i0 = L0_idx-(1);
  idx.i1 = L1_idx-(1);
  idx.i2 = L2_idx-(1);
  idx.i3 = L3_idx-(1); 
  
  % decode to give proper quantized w vector for synthesis filter
  [l_hat,w_hat] = decode_lsf(idx, lsfCb, l_hat_prev);
     
end

function [E_opt,i1,i2,i3] = lsf_vq(w,weights,l_hat_prev,L0,L1,L2,L3)
 
  assert(length(w) == 10);
  assert(all(size(l_hat_prev) == [10 4]));
 
  % assume l_hat_prev is columns of previously quantized 
  % vectors, 10x4. 
  % 4 is the history depth
  alpha = 1 - sum(L0,2);
  pred = sum(L0.*l_hat_prev,2);
  
  % find the L1 codebook index which minimizes the unweighted error
  l_hat   = (w' - pred)./alpha;
  tmp    = bsxfun(@minus,l_hat,L1);   %  error
  E      = sum(tmp .^ 2,1);
  [~,i1] = min(E);
  l_hat  = L1(:,i1);
  
  % find the L2 codebook which minimizes the weighted MSE  
  l_cand        = repmat(l_hat,1,32);
  l_cand(1:5,:) = l_cand(1:5,:) + L2;
  w_cand = bsxfun(@times,l_cand,alpha);  % generate the w vectors for each codebook .
  w_cand = bsxfun(@plus,w_cand,pred);    % candidate, and MA prediction
  w_cand = rearrange(w_cand,.0012);
  E      = weighted_mse(weights,w',w_cand);
  [~,i2] = min(E);  
  l_hat  = l_cand(:,i2);
    
  % find the L3 codebook which minimizes the weighted MSE
  l_cand         = repmat(l_hat,1,32);
  l_cand(6:10,:) = l_cand(6:10,:) + L3;
  w_cand = bsxfun(@times,l_cand,alpha);
  w_cand = bsxfun(@plus,w_cand,pred);
  w_cand = rearrange(w_cand,.0012);
  E      = weighted_mse(weights,w',w_cand);
  [~,i3] = min(E);  
  l_hat  = l_cand(:,i3);
  
  E_opt = E(i3);  
end

function weights = adapt_weights(w)
  % default weights
  weights = ones(1,10);
  
  % Equation (22)
  
  if( (w(2) - 0.04*pi - 1) <= 0)
    weights(2) = 10*(w(2) - 0.04*pi - 1)^2 + 1;
  end

  % set weights(i) = 10*(w(i+1) - w(i-1) - 1)^2 if
  % (w(i+1) - w(i-1) -1) <= 0
  diffs = w(3:10) - w(1:8) - 1;
  idx = find( diffs <= 0);
  weights(idx + 1) = 10*diffs(idx).^2;
  
  if( (-w(9) + 0.92*pi - 1) <= 0)
    weights(10) = 10*(-w(9) + 0.92*pi - 1)^2 + 1;
  end  

  % 
  weights(5:6) = 1.2*weights(5:6);

end

% Equation (21) 
% 
function E = weighted_mse(weights,w,w_cand)

  se = (bsxfun(@minus,w,w_cand)).^2;
  E  = weights * se;
end

