function retval = calcThetaSigma(a_j_tilde)  
  # disp ("Running calcThetaSigma...");
  m = length(a_j_tilde);
  n = size(a_j_tilde(1),1);
  disp ("LALALA_m = ");
  disp(m);
  if (m == 1)
    retval  = [2* a_j_tilde(1), zeros(n); zeros(n), zeros(n)];
  elseif m == 2
    retval  = [2* a_j_tilde(1), a_j_tilde(2); a_j_tilde(2), zeros(n)];
  else
    retval  = [2* a_j_tilde(1), a_j_tilde(2); a_j_tilde(2), 2* a_j_tilde(3)];
  end
  
  for i1 = 4:2:m
    if i1 ~= m
      retval = [retval, [zeros((i1/2 -1)*n, n)]; zeros(n, (i1/2 -1)*n), a_j_tilde(i1), 2*a_j_tilde(i1+1)];
    else
      retval = [a_j_tilde, [zeros((i1/2 -1)*n, n); a_j_tilde(i1)]; zeros(n, (i1/2 -1)*n), a_j_tilde(i1), zeros(n)];
    endif
  endfor
  retval = -0.5*retval;
endfunction