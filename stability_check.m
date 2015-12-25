function w_check = stability_check(w)
% LSF stability check
% section 3.2.4
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