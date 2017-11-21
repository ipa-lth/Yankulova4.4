function retval2 = P_func(l, k2, I, M)
  retval2 = I;
  if k2 != 0
    for q = (0: k2-1)
      retval2 = retval2 * ((l-q)*I + M);
    endfor      
  endif
endfunction