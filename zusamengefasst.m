% ___________Define variables
# number of dimention
n   = 3;

# define Q= R^(-1), ljapunov gebiet
Q       = sdpvar(n,n);        # Q = R^(-1)
# a_head: closed loop controled coefficients
# a_head  = sdpvar(n,1); #(a_head is the result!!!)
z       = sdpvar(n,1); 

# x: input matrix (could be e.g. initial position, velocity, ect.), right?
x       = sdpvar(n,1);

# b: input vector, what does it means?
b   = [0; 0;  1];

# I: I of n dimension
I   = eye(n);

# N = diag(-n,..., -2,-1)
N   = diag(-n:-1);

# M = diag(0,1 ..., n-1); 
M   =  diag(0: n-1);

            

# A current system
A = [0 ,  1,  0;
      0,  0,  1;
      0,  0,  -0.005];
      
# u_max
u_max = 2.5 * 10^(-5)

# m: size of a: dimention of A (nxn) ????
m = 3;

# a, A-lamda*I = 0; lamda = (a0, a1, ..., an). If we repleace A(current system)
# by A* (system we want), we have a* as a result wanted eigen vector
a = [40000; -200; 1];     #I calculate from A, is that ok?

# a: route coefficients


# z = Q*a_head;                                         #(a_head is the result!!!)
# z is propably is sth related to result, what is z means?

% ___________Define constraints
% 4.59
F = [Q >0];

% 4.60
F = [F, Q*((A.') + a*(b.')) + (A+b*(a.'))*Q - z*(b.') - b*(z.') < 0];

% 4.61
F = [F, Q*N + N*Q < 0];

% 4.62
F = [F, [[1, (x.')]; [x, Q]] > 0];

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
tmp = 0;
for i = 0 : m
  tmp = tmp + sum(i,Q,a, I, M);
end

F = [F, tmp >= 0];

function retval = P (l, k, I, M)
  retval = 1;
  for q = (0: k-1)
    retval = retval * ((l-q)*I + M);
  endfor
  return
endfunction

function retval = sum(i,Q,a, I, M)
  retval = 1;
  for k = (0: i)
    retval = retval + nchoosek(i,k)*(a.')*P(0, i-k, I, M)*Q*N*P(0,k, I, M)*a
                    -(z.')*N*P(n,i, I, M)*a;
  endfor
  return
endfunction

% ___________Define an objective
objective = 0;

% ___________optimize
optimize(F,objective)
