###############################
# Optimisation problem (4.40) #
###############################

# Define Variables
R0 = sdpvar(n, n, 'full'); #semidef
R1 = sdpvar(n, n, 'full');

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

# Problematic constraint
constraints_439 = [constraints_439,
                  [ [[1,    k0_0'];
                     [k0_0,    R0]] + 1.0/pm * [[0,    k0_1'];
                                                [k0_1,    R1]]] >= 0 ]; 


#constraint_439e_ = cvxpy.bmat([[1,               k1_0'+1.0/pm*k1_1'],
#                              [k1_0+1.0/pm*k0_1,         R0+1.0/pm*R1]]) >> 0


# Problematic constraint
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

diagnostics = optimize(constraints_439, obj_439, sdpsettings('solver', sol))
#diagnostics = optimize(constraints_439, [], sdpsettings('solver', sol))