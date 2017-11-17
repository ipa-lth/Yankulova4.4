n = 2;
# b: input vector
b = [0; 0;  1];
# k_v: position control vector

# I: I of n dimension
I = eye(n);
# M = diag(0,1 ..., n-1) 
M =  diag(1: n-1);

# N = diag(-n,..., -2,-1)
N = diag(-n,-1);

# P_(l,k)
#P_lk = @(k) P(l, k) .* (k != 0) + I .* (k == 0); 

% define Q
Q = sdpvar(n,n);                

# A uncontrolled matrix system
A = [0 ,  1,  0;
      0,  0,  1;
      0,  0,  -0.005];

# a_head: closed loop controled coefficients
a_head = a - b*k_v;

# a: route coefficients

# z
z = Q*a_head;

% Q is positivedifinited
% 4.59
F = [Q >0];

% 4.60
F = [F, Q*((A.') + a*(b.')) + (A+b*(a.'))*Q - z*(b.') - b*(z.') < 0];

% 4.61
F = [F, Q*N + N*Q < 0];

% 4.62
F = [F, [1, (x.'); x, Q] >0];

% 4.63
F = [F, [u_max^^2 -(a.')*Q*a + 2 .* (a.')*z, (z.') ; z, Q] >=0];
% 4.64
# i: belong to {0,1, ..., m}
for i = 0 : m
  F = [F, sum(i) >= 0];
end

function retval = P (l,k)
  retval = 1;
  for q = (0: k-1)
    retval = retval * ((l-q)*I + M);
  endfor
  return retval;
endfunction

function retval = sum(i)
  retval = 1;
  for k = (0: i)
    retval = retval + nchoosek(i,k)*(a.')*P(0, i-k)*Q*N*P(0,k)*a-(z.')*N*P(n,i)*a;
  endfor
  return retval;
endfunction
sdpvar