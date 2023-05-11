clc; close; clear all;

%%
syms iL1 vc L1 L2 RL1 RL2 C Rc Ro vd N1 N2 d s

%% ON

iL1_d = -vd/L1 + (RL1/L1)*iL1;
vc_d = -(1/(C*(Rc+Ro))) * vc;

vo = (Ro/Ro+Rc)*vc;

% Vectores
x = [iL1; vc];
f = [iL1_d; vc_d];
h = vo;

%% State space 

%Linealizacion del sistema (Jacobiano)
Aon = jacobian(f, x);  %df/dx
Bon = jacobian(f, vd);  %df/du
Con = jacobian(h, x);  %dh/dx
Don = jacobian(h, vd);  %dh/du

% Puntos de equlibrio REEMPLAZAR POR VALORES POSTA
% x1_eq = 0;
% x2_eq = 0;
% x3_eq = 0;
% x4_eq = 0;
% u_eq = 0;

% Xeq = [x1_eq x2_eq x3_eq x4_eq];
% delta_xeq = [0 0 0 0];
% Ueq = [u_eq];

%Calculo mmatrices en pto de equlibrio
% A = subs(A, {'x1','x2','x3', 'x4', 'u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});
% B = subs(B, {'x1','x2','x3', 'x4','u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});
% C = subs(C, {'x1','x2','x3', 'x4','u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});
% D = subs(D, {'x1','x2','x3', 'x4','u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});

% Gss : Sistema linealizado en espacio de estados
% A = double(A);
% B = double(B);
% C = double(C);
% D = double(D);

%% OFF

n = N1/N2;

vc_d = -(Ro/(C*(Ro+Rc)))*(n*iL1 + vc/Ro);

vo = Ro*(-n*iL1-C*vc_d);

%vo = vc + Rc * C * vc_d;

iL1_d = (n/L1)*(-n*iL1*RL2 + vo); %PREGUNTAR POR EL MENOS QUE DEBERIA IR AL PRINCIPIO DE TODO

%iL1_d = -n/L1*(1-(Rc*Ro)/(Ro*(Rc+Ro))) * vc + n*n/L1 * (RL2 + (Rc*Ro)/(Rc+Ro)) *iL1;

%  iL1_d = n/L1*(-n*(Ro+RL2)*iL1 + Ro*((vc + n * iL1*Ro)/(Rc+Ro)));
%  vc_d = -(1/(C*(Ro+Rc)))*(vc + n*Ro*iL1);
%  
%  vo = vc - Rc * C * vc_d;

% Vectores
x = [iL1; vc];
f = [iL1_d; vc_d];
h = vo;

%% State space 

%Linealizacion del sistema (Jacobiano)
Aoff = jacobian(f, x);  %df/dx
Boff = jacobian(f, vd);  %df/du
Coff = jacobian(h, x);  %dh/dx
Doff = jacobian(h, vd);  %dh/du

%% Promedios 

A_a = Aon*d+Aoff*(1-d);
B_a = Bon*d+Boff*(1-d);
C_a = Con*d+Coff*(1-d);
D_a = Don*d+Doff*(1-d);

%% Funcion transferencia

Phi=inv(s*eye(2)-A_a);

H=C_a*Phi*B_a+D_a;

X = Phi * B_a * vd; % los estados considerando D de flyback


%% Transferencia vo respecto a d

H_d = C_a * Phi * ((Aon - Aoff) * X + (Bon - Boff)*vd) + (Con-Coff)*X;


%% ceros y polos

[N, D] = numden(H_d)

%%

N_coeffs = coeffs(N, s)
D_coeffs = coeffs(D, s)

%% Reemplazo valores

Ro_eq = 100;
L1_eq = 10e-3;
L2_eq = 16e-6;
L3_eq = 100e-6;
C_eq = 4.7e-6;
RL1_eq = 1;
RL2_eq = 1; 
Rc_eq = 0.5;
d_eq = 0.5;
vd_eq = 300;
N1_eq = 50;
N2_eq = 2;

Hd_values = subs(H_d, {'L1', 'L2', 'RL1', 'RL2', 'C', 'Rc', 'Ro', 'vd', 'N1', 'N2', 'd'}, {L1_eq L2_eq RL1_eq RL2_eq C_eq Rc_eq Ro_eq vd_eq N1_eq N2_eq d_eq});

%%
[N, D] = numden(Hd_values)
N_coeffs = coeffs(N, s)
D_coeffs = coeffs(D, s)

%%
sys = tf(double(N_coeffs), double(D_coeffs))
step(sys)
figure('Name', 'Bode'); bode(sys)
margin(sys)
grid on