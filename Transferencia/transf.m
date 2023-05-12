clc; close; clear all;

%%
syms iL1 vc L1on L1off L2 RL1 RL2 C Rc Ro vd N1 N2 d s R1 R2 C1 C2


%% ON

iL1_d = -vd/L1on + (RL1/L1on)*iL1
vc_d = -(1/(C*(Rc+Ro))) * vc

vo = (Ro/(Ro+Rc))*vc

% Vectores
x = [iL1; vc];
f = [iL1_d; vc_d];
h = vo;

%% State space 

%Linealizacion del sistema (Jacobiano)
Aon = jacobian(f, x)  %df/dx
Bon = jacobian(f, vd)  %df/du
Con = jacobian(h, x)  %dh/dx
Don = jacobian(h, vd)  %dh/du

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

vc_d = -(Ro/(C*(Ro+Rc)))*(n*iL1 + vc/Ro)

%vo = Ro*(-n*iL1+C*vc_d); 

vo = vc + Rc * C * vc_d

iL1_d = (1/(n*L2))*(-n*iL1*RL2 + vo)

% iL1_d = -n/L1*(1-(Rc*Ro)/(Ro*(Rc+Ro))) * vc + n*n/L1 * (RL2 + (Rc*Ro)/(Rc+Ro)) *iL1;

% % %  iL1_d = n/L1off*(-n*(Ro+RL2)*iL1 + Ro*((vc + n * iL1*Ro)/(Rc+Ro)));
% % %  vc_d = -(1/(C*(Ro+Rc)))*(vc + n*Ro*iL1);
% % %  
% % %  vo = Ro*(-n*iL1-C*vc_d); 

% Vectores
x = [iL1; vc];
f = [iL1_d; vc_d];
h = vo;

%% State space 

%Linealizacion del sistema (Jacobiano)
Aoff = jacobian(f, x)  %df/dx
Boff = jacobian(f, vd)  %df/du
Coff = jacobian(h, x)  %dh/dx
Doff = jacobian(h, vd)  %dh/du

%% Promedios 

A_a = Aon*d+Aoff*(1-d);
B_a = Bon*d+Boff*(1-d);
C_a = Con*d+Coff*(1-d);
D_a = Don*d+Doff*(1-d);


%% Funcion transferencia

Phi=inv(s*eye(2)-A_a);

H=C_a*Phi*B_a;

X = Phi * B_a * vd; % los estados considerando D de flyback

X = subs(X, s, 0);

%% Transferencia vo respecto a d

H_d = C_a * Phi * ((Aon - Aoff) * X + (Bon - Boff)*vd) + (Con-Coff)*X;



%% ceros y polos

[N, D] = numden(H_d)

%%

%% Reemplazo valores

Ro_eq = 100;
L1on_eq = 11e-3;
L1off_eq = 10e-3;
L2_eq = 16e-6;
L3_eq = 100e-6;
C_eq = 1e-6;
RL1_eq = 80e-3;
RL2_eq = 1e-3;
RL3_eq = 1e-3; 
Rc_eq = 0.5;
d_eq = 0.35;
vd_eq = 300;
N1_eq = 50;
N2_eq = 2;
N3_eq = 5;

Hd_values = subs(H_d, {'L1on', 'L1off', 'L2', 'RL1', 'RL2', 'C', 'Rc', 'Ro', 'vd', 'N1', 'N2', 'd'}, {L1on_eq L1off_eq L3_eq RL1_eq RL3_eq C_eq Rc_eq Ro_eq vd_eq N1_eq N3_eq d_eq});

%%
[N, D] = numden(Hd_values)
[N_coeffs, a] = coeffs(N, s)
[D_coeffs, b] = coeffs(D, s)

%%
sys = tf(double(N_coeffs), double(D_coeffs))
step(sys)
figure('Name', 'Bode'); bode(sys)
margin(sys)
grid on

%%
step(sys)
pzplot(sys)
nyqlog(sys)
%%
EA= -(1+R2*C1*s)/(R1*(R2*C1*C2*s*s+(C1+C2)*s));
boost=54
fs=100e3
f0=fs/5
fz= f0/10 %tan(boost*pi/360 + pi/4)
fp= f0*10 %tan(boost*pi/360 + pi/4)*f0

%%
R1_eq = 2.2e3;
R2_eq = 15e3;
C1_eq =  4.7e-9 %1/(2*pi*R2_eq*fz);
C2_eq =  1e-12 %1/(2*pi*R2_eq*fp);
EA_values = subs(EA, {'R1', 'R2', 'C1', 'C2'}, {R1_eq R2_eq C1_eq C2_eq})
[N, D] = numden(EA_values)
[N_coeffs, a] = coeffs(N, s)
[D_coeffs, b] = coeffs(D, s)
D_coeffs = [D_coeffs(1), D_coeffs(2), 0];
sysEA = tf(double(N_coeffs), double(D_coeffs))
figure('Name', 'Bode2'); bode(sysEA)
margin(sysEA)
grid on

%%
nyqlog(sysEA)
%%
L = -sysEA*sys
figure('Name', 'Bode2'); bode(L)
margin(L)
grid on
%%
nyqlog(sysEA)