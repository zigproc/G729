function [l_hat,w_hat,idx] = quantize_lsf(w,l_hat_prev,lsp_codebook)

  % Section 3.2.4
  assert(length(w) == 10);

  l_hat = zeros(1,10);
  
  weights = adapt_weights(w);
  
  % test versions of codebooks
  L0_0 = randn(10,4);
  L0_1 = randn(10,4);    % selecting L0_0 vs. L0_1 implies 1 bit
  L1   = randn(10,128);  % 7 bit
  L2   = randn(10,32);   % 5 bit
  L3   = randn(10,32);   % 5 bit
  
  % compute prediction based on MA
  [E0,l1_0,l2_0,l3_0] = lsp_vq(w,l_hat_prev,L0_0,L1,L2,L3);
  [E1,l1_1,l2_1,l3_1] = lsp_vq(w,l_hat_prev,L0_1,L1,L2,L3);
 
  if( E0 < E1)
    L0_idx = 0;
    L0     = L0_0;
    L1_idx = l1_0;
    L2_idx = l2_0;
    L3_idx = l3_0;
  else
    L0_idx = 1;
    L0     = L0_1;
    L1_idx = l1_1;
    L2_idx = l2_1;
    L3_idx = l3_1;  
  end
    

    
  l_opt = L1(:,L1_idx);
  l_opt(1:5) = l_opt(1:5) + L2(:,L2_idx);
  l_opt(6:10) = l_opt(6:10) + L3(:,L3_idx);
  
  % rearrange l_opt to guarantee min dist of .0006
  l_opt = rearrange(l_opt,.0006);
  
  % gen best LSF vector, perform stability check
  w_opt = (1 - sum(L0,2))*l_opt + sum(L0.*l_hat_prev,2);
  w_opt = stability_check(w_opt);

  % return values  
  l_hat = l_opt;
  w_hat = w_opt;
  
  idx.i0 = L0_idx;
  idx.i1 = L1_idx;
  idx.i2 = L2_idx;
  idx.i3 = L3_idx;    
  
  
  % 
end

function [E_opt,L1_idx,L2_idx,L3_idx] = lsp_vq(w,l_hat_prev,L0,L1,L2,L3)
 
  % assume l_hat_prev is columns of previously quantized 
  % vectors, 10x4. 
  % 4 is the history depth
  scale = 1 - sum(L0,2);
  pred = sum(L0.*l_hat_prev,2);
  
  % TBD do below for both MA predictors, and select best
  
  % find the L1 codebook index which minimizes the unweighted error
  w_cand = repmat(scale,1,128).*L1 + repmat(pred,1,128);
  E      = weighted_mse(ones(1,10),repmat(w',1,128),w_cand);
  [~,L1_idx] = min(E);
  w_L1   = w_cand(:,L1_idx);
  
  % find the L2 codebook which minimizes the weighted MSE
  w_cand = repmat(w_L1,1,32);
  w_cand(1:5,:) = w_cand(1:5,:) + L2;
  % TBD rearrange to guarantee min dist of .0012
  E = weighted_mse(weights,repmat(w',1,32),w_cand);
  [~,L2_idx] = min(E);
  w_L2 = w_cand(:,L2_idx);
  
  % find the L3 codebook which minimizes the weighted MSE
  w_cand = repmat(w_L2,1,32);
  w_cand(6:10,:) = w_cand(6:10,:) + L3;
  % TBD rearrange to guarantee min dist of .0012
  E = weighted_mse(weights,repmat(w',1,32),w_cand);
  [~,L3_idx] = min(E);
  w_L3 = w_cand(:,L3_idx);  
  
  E_opt = E(L3_idx);  
end

function weights = adapt_weights(w)
  % default weights
  weights = ones(1,10);
  
  % Equation (22)
  
  if( (w(2) - 0.04*pi - 1) <= 0)
    weights(2) = 10*(w(2) - 0.04*pi - 1)^2 + 1;
  end

  % set weights(i) = 10*(w(i+1) - w(i-1) - 1)^2 if
  % (w(i+1) - w(i-1) - -1) <= 0
  diffs = w(3:10) - w(1:8) - 1;
  idx = find( diffs <= 0);
  weights(idx - 1) = 10*diffs(idx).^2;
  
  if( (-w(9) + 0.92*pi - 1) <= 0)
    weights(10) = 10*(-w(9) + 0.92*pi - 1)^2 + 1;
  end  

  % 
  weights(5:6) = 1.2*weights(5:6);

end

% Equation (21) 
% 
function E = weighted_mse(weights,w,w_hat)

  se = (w - w_hat).^2;
  
  E = weights * se';
end

% rearrange candidate vectors.
% x is DxN, where D is dimension, and N number of vectors
% J is min distance
function y = rearrange(x,J)
  [D,N] = size(x);
  
  y = zeros(D,N);
  
  for i=1:N
    v = x(:,i);
    for k=2:D
      if(v(k-1) > v(k) - J)
        v(k-1) = (v(k) + v(k-1) - J)/2;
        v(k)   = (v(k) + v(k-1) + J)/2;
      end
    end
    
    y(:,i) = v;
end


function w_check = stability_check(w)

  w_check = sort(w);
  if(w_check(1) < .005)
    w_check(1) = .005;
  end

  diff = w_check(2:10) - w_check(1:9);
  idx  = find(diff < 0);
  w_check(idx + 1) = w_check(idx) + .0391;
  
  if( w_check(10) > 3.135)
    w_check(10) = 3.135;
  end
  
end
