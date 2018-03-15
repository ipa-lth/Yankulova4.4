clear

# loading YALMIP
init_yalmip

# loading Fusions Reaktor
init_fusion_reactor

# Init values (may come from optimization loop)
pmin = 1/20
mu1 = 1.5
zeta0 = 2.5
zeta1 = 5



######################
# Löse A.15 mit mu=1 #
######################

# Define Variables
Q = sdpvar(n, n);  # Q = R^(-1)
#Q = sdpvar(n, n, 'full');  # Q = R^(-1)

z0 = sdpvar(n, 1);
z1 = sdpvar(n, 1);


# Bisection parameter
#g = cvxpy.Parameter(sign='positive')
sdpvar t;

mu = 1;
#mu = cvxpy.Parameter(sign='positive')
          
# Define Constraints
constraints_A15 = [Q >= 0]; # Q is Semidefinite

# (A.10)
for i = 1:size(X0,2)
  constraints_A15 = [constraints_A15,
                    [Q,       [X0(:,i)];
                    [X0(:,i)',       1]] >= 0];
end

# (A.11)
constraints_A15 = [constraints_A15,
                  [Q,   z0;
                   z0',  1 ] >= 0];

# (A.12
constraints_A15 = [constraints_A15,
                  [Q,   z1;
                   z1',  mu**2 ] >= 0];

# (A.13)
constraints_A15 = [constraints_A15,
                  Q*A+A'*Q-b*z0'-z0*b'<=0];

# (A.14)
constraints_A15 = [constraints_A15,
                  Q*A + A'*Q - b*z1' - z1*b' + 2*t*Q<=0];


ops = sdpsettings('solver','bisection','bisection.solver','sdpt3');
Objective = -t; #Maximize
diagnostics = optimize(constraints_A15, Objective, ops);
#test_it(value(Q), value(z0), value(z1), mu, value(t), A, b, X0)
#check(constraints_A15)

disp("l*(mu=1)=")
l1m0 = inv(value(Q))*value(z1)

disp("lambda^{hat}(p=1)=")
roots_m0_p1 = eig(A-b*l1m0')'



#######################
# Löse A.15 mit mu>=1 #
#######################

# Define Variables
Q = sdpvar(n, n);  # Q = R^(-1)
#Q = sdpvar(n, n, 'full');  # Q = R^(-1)

z0 = sdpvar(n, 1);
z1 = sdpvar(n, 1);


# Bisection parameter
#g = cvxpy.Parameter(sign='positive')
sdpvar t;

mu = mu1;
#mu = cvxpy.Parameter(sign='positive')
          
# Define Constraints
constraints_A15 = [Q >= 0]; # Q is Semidefinite

# (A.10)
for i = 1:size(X0,2)
  constraints_A15 = [constraints_A15,
                    [Q,       [X0(:,i)];
                    [X0(:,i)',       1]] >= 0];
end

# (A.11)
constraints_A15 = [constraints_A15,
                  [Q,   z0;
                   z0',  1 ] >= 0];

# (A.12
constraints_A15 = [constraints_A15,
                  [Q,   z1;
                   z1',  mu**2 ] >= 0];

# (A.13)
constraints_A15 = [constraints_A15,
                  Q*A+A'*Q-b*z0'-z0*b'<=0];

# (A.14)
constraints_A15 = [constraints_A15,
                  Q*A + A'*Q - b*z1' - z1*b' + 2*t*Q<=0];


ops = sdpsettings('solver','bisection','bisection.solver','sdpt3');
Objective = -t; #Maximize
diagnostics = optimize(constraints_A15, Objective, ops);
#test_it(value(Q), value(z0), value(z1), mu, value(t), A, b, X0)
#check(constraints_A15)

disp("l*(mu=1)=")
l1m1 = inv(value(Q))*value(z1)

disp("lambda^{hat}(p=1)=")
roots_m1_p1 = eig(A-b*l1m1')'

disp("lambda^{hat}(p=pmin)=")
roots_m0_pmin = zeta0 * real(roots_m0_p1) + 1i*imag(roots_m0_p1)
#roots_m0_pmin = zeta0 * roots_m0_p1  # Shifting on both real and imaginary axis in zeta

disp("lambda^{hat}_*(p=pmin)=")
roots_m1_pmin = zeta1 * real(roots_m1_p1) + 1i*imag(roots_m1_p1)
#roots_m1_pmin = zeta0 * roots_m1_p1  # Shifting on both real and imaginary axis in zeta

[k0_0, k0_1] = k_explizit_Ab2(roots_m0_p1, roots_m0_pmin, pmin, A, b)
figure;
plot_moving_poles(A, b, c, d, k0_0, k0_1, pmin)

[k1_0, k1_1] = k_explizit_Ab2(roots_m1_p1, roots_m1_pmin, pmin, A, b)
figure;
plot_moving_poles(A, b, c, d, k1_0, k1_1, pmin)

###############################
# Optimisation problem (4.40) #
###############################

# Define Variables
R0 = sdpvar(n, n); #semidef
R1 = sdpvar(n, n);

G0 = sdpvar(n, n, 'skew'); # G -> skew
G1 = sdpvar(n, n, 'skew'); # G_* -> skew

D0 = sdpvar(n, n, 'sym'); # D -> sym
D1 = sdpvar(n, n, 'sym'); # D_* -> sym

# Define Parameters
pm = pmin;

alpha = 1-pm;
beta = 1+pm;

# TODO: Check if Â_0 and Â_1 are correct (find proof)
A0_0 = A-b*k0_0';                                               # A^{hat}_0
A0_1 = b*k0_1';                                                  # A^{hat}_0
#A0_1 = A-b*k0_1'                                               # A^{hat}_1

T0_0 = A0_1'*R1 + R1*A0_1;                                       # Theta_0
T0_1 = A0_0'*R1 + R1*A0_0 + A0_1'*R0 + R0*A0_1;                 # Theta_1
T0_2 = A0_0'*R0 + R0*A0_0;                                       # Theta_2

# TODO: Check if Â_{*,0} and Â_{*,1} are correct (find proof)
A1_0 = A-b*k1_0';                                                # A^{hat}_{*, 0}
A1_1 = b*k1_1';                                                  # A^{hat}_{*, 0}
#A1_1 = A-b*k1_1'                                               # A^{hat}_{*, 1}

T1_0 = A1_1'*R1 + R1*A1_1;                                      # Theta_{*,0}
T1_1 = A1_0'*R1 + R1*A1_0 + A1_1'*R0 + R0*A1_1;                # Theta_{*,1}
T1_2 = A1_0'*R0 + R0*A1_0;                                      # Theta_{*,2}

constraints_439 = [R0 >= 0];

constraints_439 = [constraints_439,
                  [R0 + R1 >= 0]];
             
constraints_439 = [constraints_439,
                  [ [-D0,  G0;
                      G0', D0] - [T0_0 + 0.5*beta*T0_1 + 0.25*beta**2*T0_2, 0.25*alpha*(T0_1 + beta*T0_2);
                                  0.25*alpha*(T0_1 + beta*T0_2),            0.25*alpha**2*T0_2           ]] >= 0 
                  ];

for i = 1:size(X0,2)
  constraints_439 = [constraints_439,
                    [1.0 - X0(:,i)'*(R0 + R1)*X0(:,i)] >= 0 ];
end

constraints_439 = [constraints_439,
                  [ [[1,    k0_0'];
                     [k0_0,    R0]] + 1.0/pm * [[0,    k0_1'];
                                                [k0_1,    R1]]] >= 0 ]; 

#constraint_439e_ = cvxpy.bmat([[1,               k0_0'+1.0/pm*k0_1'],
#                              [k0_0+1.0/pm*k0_1,         R0+1.0/pm*R1]]) >> 0

constraints_439 = [constraints_439,
                  [ [[1,    k0_0'];
                     [k0_0,    R0]] + [[0,    k0_1'];
                                       [k0_1,    R1]]] >= 0 ]; 

#constraint_439f_ = cvxpy.bmat([[1,    k0_0'+k0_1'],
#                              [k0_0+k0_1,     R0+R1]]) >> 0                                                

constraints_439 = [constraints_439,
                  [ [-D1,  G1;
                      G1', D1] - [T1_0 + 0.5*beta*T1_1 + 0.25*beta**2*T1_2, 0.25*alpha*(T1_1 + beta*T1_2);
                                  0.25*alpha*(T1_1 + beta*T1_2),            0.25*alpha**2*T1_2           ]] >= 0 
                  ];

#check(constraints_439)

obj_439 = trace(R0+1.0/pm*R1); # objective

diagnostics = optimize(constraints_439, obj_439, sdpsettings('solver', 'sdpt3'))

disp("end")