function ret_sum = func_constraint464(i, Q, a, M, N, n, z)  
  ret_sum = 0;
  I = eye(n);
  
  for k = (0: i)
    ret_sum = ret_sum + nchoosek(i,k)*a'*func_P(0, i-k, I, M)*Q*N*func_P(0, k, I, M)*a-z'*N*func_P(n, i, I, M)*a;
  endfor
endfunction