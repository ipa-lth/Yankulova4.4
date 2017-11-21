function retval = sum_func(i, Q, a, I, M, N, n, z)  
  disp ("Running sum_func...");
  retval = 1;
  for k = (0: i)
    retval = retval  + nchoosek(i,k)*(a.')*P_func(0, i-k, I, M)*Q*N*P_func(0,k, I, M)*a
                    - (z.')*N*P_func(n,i, I, M)*a;
  endfor
  #disp ("sum_func: ");
  #disp (retval);  
endfunction