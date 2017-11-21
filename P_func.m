function retval = P_func(l, k, I, M)
  disp("Running P_func...")
  retval = I;
  if k != 0
    for q = (0: k-1)
      retval = retval * ((l-q)*I + M);
    endfor      
  endif
  #disp(retval);
endfunction