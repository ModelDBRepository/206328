%% Kinetic parameters of state transitions during AMPA receptor activation
%% Milstein-Nicoll Scheme
%%

k_RtoRG=Glu_conc*1*10^7;%in M^-1.s^-1.
k_RGtoR=5*10^4;%in s^-1.

k_RGtoC1=3.65*10^4;%in s^-1.
k_C1toRG=4.55*10^3;%in s^-1.

k_C1toC2=3*10^2;%in s^-1.
k_C2toC1=1*10^4;%in s^-1.

k_C1toO1=1*10^4;%in s^-1.
k_O1toC1=6*10^3;%in s^-1.

k_C2toO2=1*10^4;%in s^-1.
k_O2toC2=6*10^3;%in s^-1.

k_C1toD1=1.1*10^3;%in s^-1.
k_D1toC1=1;%in s^-1.

k_C2toD2=3*10^2;%in s^-1.
k_D2toC2=10;%in s^-1.

%Number of states=8;
%Sequence of variables
%[R RG C1 C2 O1 O2 D1 D2]

%Transition from R to rest
Row1=[0 k_RtoRG 0 0 0 0 0 0];
%Transition from RG to rest
Row2=[k_RGtoR 0 k_RGtoC1 0 0 0 0 0];
%Transition from C1 to rest
Row3=[0 k_C1toRG 0 k_C1toC2 k_C1toO1 0 k_C1toD1 0];
%Transition from C2 to rest
Row4=[0 0 k_C2toC1 0 0 k_C2toO2 0 k_C2toD2];
%Transition from O1 to rest
Row5=[0 0 k_O1toC1 0 0 0 0 0];
%Transition from O2 to rest
Row6=[0 0 0 k_O2toC2 0 0 0 0];
%Transition from D1 to rest
Row7=[0 0 k_D1toC1 0 0 0 0 0];
%Transition from D2 to rest
Row8=[0 0 0 k_D2toC2 0 0 0 0];


Transition_matrix=[Row1;Row2;Row3;Row4;Row5;Row6;Row7;Row8];
