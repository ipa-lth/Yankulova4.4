function retval = sum_func(i, Q, a, I, M, N)  
  disp ("painkiller beginning: ");
  retval = 1;
  M_M = M;
  for k = (0: i)
    #retval = retval + nchoosek(i,k)*P_func(0, i-k, I, M_M)*P_func(0,k, I, M_M)*(a.')*Q*N*a
    #                -P_func(n,i, I, M_M)*(z.')*N*a;
    retval = retval + 1;
  endfor
  disp ("painkiller: ");
  disp (retval);  
endfunction