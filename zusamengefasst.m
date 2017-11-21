# loading path
addpath ("/home/ipa325/Downloads/YALMIP-master")
addpath ("/home/ipa325/Downloads/YALMIP-master/extras")
addpath ("/home/ipa325/Downloads/YALMIP-master/solvers")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules/parametric")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules/global")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules/sos")
addpath ("/home/ipa325/Downloads/YALMIP-master/operators")


# final aim is to find the Q and z. A, b should be in canonical form 
# ("Steuerung normal form")
# ___________Define variables
# number of dimention
n   = 3;

# define Q= R^(-1), ljapunov gebiet
Q       = sdpvar(n,n);        # Q = R^(-1)
# a_head: closed loop controled coefficients
# a_head  = sdpvar(n,1); #(a_head is the result!!!)
z       = sdpvar(n,1); 
#z       = [1; 1; 1]; 

# x: one line of Xi_o: initial state of the system
Xi_o    = [ 10,  0.05, 0.0046;
            10,  0.05, -0.0046;
            10,  -0.05, 0.0046;
            10,  -0.05, -0.0046;
            -10,  0.05, 0.0046;
            -10,  0.05, -0.0046;
            -10,  -0.05, 0.0046;
            -10,  -0.05, -0.0046;
            ]
            
x       = Xi_o(1,:).';
disp("x = ");
disp(x);
disp("Xi_o NoR = ");
Xi_o_num = size (Xi_o, 1);
disp(Xi_o_num);

# b: input vector, what does it means?
b   = [0; 0;  1];

# I: I of n dimension
I   = eye(n);

# N = diag(-n,..., -2,-1)
N   = diag(-n:-1);

# M = diag(0,1 ..., n-1); 
M   =  diag(0: n-1);

            

# A: current system matrix in canonical form ("Steuerung normal form")
A = [0 ,  1,  0;
      0,  0,  1;
      0,  0,  -0.005];
      
# a is the last line of canonical form of A
A_nol = size(A,1);      # number of line of A
disp ("a = ");
a = A (A_nol,:).';     #I calculate from A, is that ok? => sure, a = last line of A
disp (a);
      
# u_max
u_max = 2.5 * 10^(-5)

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
for i = 1: Xi_o_num
  # F = [F, [[1, (x.')]; [x, Q]] > 0];
  F = [F, [[1, Xi_o(i,:)]; [(Xi_o(i,:).'), Q]] > 0];
end


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
#disp ("Q = ");
#disp(Q);
#disp ("a = ");
#disp(a);
#disp ("I = ");
#disp(I);
#disp ("M = ");
#disp(M);
#disp ("N = ");
#disp(N);

# m: size of a: dimention of A (nxn) ????
m = 3;
for i = 0 : m
  tmp = tmp + sum_func(i, Q, a, I, M, N, n, z);
end

F = [F, tmp >= 0];

% ___________Define an objective
objective = 0;

% ___________optimize
optimize(F,objective)
Pfeasible = value(Q);
