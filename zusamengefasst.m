n = 2;
# diagonal matrix
M =  diag(1: n-1);
# l
# k: k positive integer
# I
# M
# P_(l,k)
P_lk = @(k) P(l, k) .* (k != 0) 
            + I .* (k == 0); 

% define Q
Q = sdpvar(n,n);                
% A, a, b, z
% Q is positivedifinited
% 4.59
F = [P >=0];

% 4.60
F = [F, Q*(A' + a*b') + (A+b*a')*Q - z*b' - b*z' <= 0];

% 4.61
F = [Q*N + N*Q <= 0];

% 4.62
F = [F, [1, x';x, Q] >=0];

% 4.63
F = [F, [u_max^^2 -a'*Q*a + 2*a'z, z';z, Q] >=0];
% 4.64
F = [F, sum]

function retval = P (l,k)
  retval = 1;
  for q = (0: k-1)
    retval = retval * (l-q)*I + M;
  endfor
  return
endfunction