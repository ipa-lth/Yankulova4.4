# loading YALMIP
addpath ("/home/ipa325/Downloads/YALMIP-master")
addpath ("/home/ipa325/Downloads/YALMIP-master/extras")
addpath ("/home/ipa325/Downloads/YALMIP-master/solvers")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules/parametric")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules/global")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules")
addpath ("/home/ipa325/Downloads/YALMIP-master/modules/sos")
addpath ("/home/ipa325/Downloads/YALMIP-master/operators")

%pkg load signal
%pkg load control


# ___________Define variables
# number of dimention, define Q= R^(-1), ljapunov gebiet
n       = 3;
Q       = sdpvar(n,n);        # Q = R^(-1)
z       = sdpvar(n,1); 

# x: one line of Xi_o
# Xi_o: initial state matrix of the system

Xi_o      = [ 20,  10, 10;
            20,  10, -10;
            20,  -10, 10;
            20,  -10, -10;
            -20,  10, 10;
            -20,  10, -10;
            -20,  -10, 10;
            -20,  -10, -10;
            ]
            

disp("Xi_o NoR = ");
Xi_o_num = size (Xi_o, 1);
disp(Xi_o_num);

# I: I of n dimension
I   = eye(n);

# N = diag(-n,..., -2,-1)
N   = diag(-n:-1);

# M = diag(0,1 ..., n-1); 
M   =  diag(0: n-1);

A1  = [    0,        1,       0;
         -10,   -1.167,      25;
           0,        0,    -0.8];

b1  = [0; 0; 2.4];

c1  = [1; 0; 0];     

# A, b,c  should be in canonical form ("Steuerung normal form")
[A,b,c,d,Ti,Qi] = get_Steuerungsnormalform(A1, b1, c1, 0);
%disp("dimension of T");
%disp( size(Ti,1));
%disp( size(Ti,2));

Xi_o_R = (Ti*(Xi_o.')).';
disp(Xi_o_R);
x_R         = Xi_o_R(1,:).';
disp("x_R = ");
disp(x_R);

# a is the last line of canonical form of A
A_nol           = size(A,1);      # number of line of A

a               = (-1)*A (A_nol,:).';

disp ("a = ");
disp (a);
      
# u_max
u_max           = 10.5


################################
#
# Define constraints
#
################################

%===========================4.59
F = [Q >= 0];

%===========================4.60
F = [F, Q*((A.') + a*(b.')) + (A+b*(a.'))*Q - z*(b.') - b*(z.') <= 0];

%===========================4.61
F = [F, Q*N + N*Q <= 0];

%===========================4.62
for i = 1: Xi_o_num
  F = [F, [[1, Xi_o_R(i,:)]; [(Xi_o_R(i,:).'), Q]] > 0];
end
%F = [F, [[1, x_R.']; [x_R, Q]] > 0];

%===========================4.63
%disp(size(u_max^2 -(a.')*Q*a + 2 .* (a.')*z));
%disp(size(z.'));
%disp(size(z));
%disp(size(Q));
%disp(size([(u_max^2 -(a.')*Q*a + 2 * (a.')*z), (z.')]));
%disp(size( [ z, Q]));
F = [F,  [[(u_max^2 -(a.')*Q*a + 2 * (a.')*z), (z.')]; [ z, Q]] >=0];

%===========================4.64
# i: belong to {0,1, ..., m}
tmp = 0;
%disp ("Q = ");
%disp(Q);
%disp ("a = ");
%disp(a);
disp ("I = ");
disp(I);
disp ("M = ");
disp(M);
disp ("N = ");
disp(N);

# m: size of a: dimention of A (nxn) ????
m = 5;
tmp = 0;
for i = 0 : m
  tmp = tmp + sum_func(i, Q, a, I, M, N, n, z);
  F = [F, tmp >= 0];
end

%===========================4.5.2: Bisection
F_preserved = F;
Beta_lower = 0.1;                   % this is supposed to be feasible. However, problem is I can't find the feasible Beta
Beta_upper = Beta_lower * 2;        % make sure Beta_lower > 0
Beta_work  = Beta_lower; 
%Beta = 100;
F = [F_preserved, Q*(A.' + a*b.') + (A +b*a.')*Q - z*b.' -b*z.' + 2*Beta_upper*Q < 0];
sol = optimize(F,[], sdpsettings('solver', 'SDPT3'))
while (sol.problem==0)
    Beta_upper = Beta_upper * 2;
    F = [F_preserved, Q*(A.' + a*b.') + (A +b*a.')*Q - z*b.' -b*z.' + 2*Beta_upper*Q < 0];
    sol = optimize(F,[], sdpsettings('solver', 'SDPT3'))
end


################################
#
# Define an objective
#
################################

#===========================4.5.1
%objective = -logdet(Q);
%objective = -geomean(Q);

################################
#
# Optimize
#
################################
%sol = optimize(F,objective, sdpsettings('solver', 'SDPT3'))

#===========================4.5.2



% Bisection code for finding B_work, which is max
tol = 0.01;
%Beta_work  = Beta_lower; 
while (Beta_upper - Beta_lower)>tol
  Beta_test = (Beta_upper + Beta_lower)/2;
  disp([Beta_lower Beta_upper Beta_test])
  F = [F_preserved, Q*(A.' + a*b.') + (A +b*a.')*Q - z*b.' -b*z.' + 2*Beta_test*Q < 0];
  sol = optimize(F,[], sdpsettings('solver', 'SDPT3'))
  if sol.problem==1
    Beta_upper = Beta_test;
  elseif sol.problem==0
    Beta_lower = Beta_test;
    Beta_work = Beta_test;
  else
    display('Something else happens!!!')
    break
 end
end

%diagnostics = bisection(F,objective, sdpsettings('solver','bisection','bisection.solver','mosek'))

Q_value = value(Q);
R1  = Q_value^-1;
z_value = value(z);

disp ("Q = ");
disp(Q_value);
disp ("R1 = ");
disp(R1);
disp ("a_head = ");
disp(R1*z_value);
disp ("Ti = ");
disp(Ti);