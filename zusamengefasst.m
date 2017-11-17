% ___________Define variables

n = 3;

# x
x = sdpvar(n,1);

# b: input vector
b = [0; 0;  1];

# I: I of n dimension
I = eye(n);
# M = diag(0,1 ..., n-1) 
M =  diag(1: n-1);

# N = diag(-n,..., -2,-1)
N = diag(-n:-1);

% define Q
Q = sdpvar(n,n);                

# A uncontrolled matrix system
A = [0 ,  1,  0;
      0,  0,  1;
      0,  0,  -0.005];
      
# u_max
u_max = 2.5 * 10^(-5)

# m: size of a_head
m = 3;

# a
a = sdpvar(n,1);
#a = [4.4469*10^-8;  2.3073*10^-5; 4.9148*10^-3];

# a_head: closed loop controled coefficients
a_head = [4.4469*10^-8;  2.3073*10^-5; 4.9148*10^-3];

# a: route coefficients

# z
z = Q*a_head;

% ___________Define constraints
% 4.59
F = [Q >0];

% 4.60
F = [F, Q*((A.') + a*(b.')) + (A+b*(a.'))*Q - z*(b.') - b*(z.') < 0];

% 4.61
F = [F, Q*N + N*Q < 0];

% 4.62
F = [F, [1, (x.'); x, Q] >0];

% 4.63
#disp(size(u_max^2 -(a.')*Q*a + 2 .* (a.')*z));
#disp(size(z.'));
#disp(size(z));
#disp(size(Q));
#disp(size([(u_max^2 -(a.')*Q*a + 2 * (a.')*z), (z.')]));
#disp(size( [ z, Q]));
F = [F,  [[(u_max^2 -(a.')*Q*a + 2 * (a.')*z), (z.')]; [ z, Q]] >=0];
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
  return
endfunction

function retval = sum(i)
  retval = 1;
  for k = (0: i)
    retval = retval + nchoosek(i,k)*(a.')*P(0, i-k)*Q*N*P(0,k)*a-(z.')*N*P(n,i)*a;
  endfor
  return
endfunction

% ___________Define an objective
objective = x;

% ___________optimize
optimize(F,objective)
