%% Numerical simulation of the glutamate diffusion with non-zero sigma_Dglut.
%%
%Critical geometric dimensions involved in simplified synapse geometry.
h=20;%height of cleft in nm.
eps=5;%height of the local volume above the PSD in nm.
Rpsd=100;%radius of the PSD in nm.
Rcleft=200;%radius of the entire cleft space in nm.
Dist_absrptn=40;%distance of the absorbtion boundary from Rcleft in nm.
Rabs=Rcleft+Dist_absrptn;%resultant radius of absorption in nm.

%Defining time-dimension.
Ti=0;Tf=1000;%time scale is in microseconds. So the diffusion is observed for 1 millisecond interval.
Ntme=10001;%No. of temporal discretisations.
dt=(Tf-Ti)/(Ntme-1);%time-step of brownian simulation is 0.1 microsecond.
t=linspace(Ti,Tf,Ntme);%time vector on microseconds scale.

%Glutamate diffusion parameters and its distribution features.
Ntrns=2000;%No. of transmitter molecules released from one vesicle.
mu_D=200;%mean diffusion coefficient with units in nm^2.us^-1.
sigma_D=40;%standard deviation around mean diffusion coefficient with units in nm^2.us^-1.
theta=(sigma_D^2)/mu_D;%Scale parameter of gamma distribution.
kappa=mu_D/theta;%Shape parameter of gamma distribution.
D=gamrnd(kappa,theta,Ntrns,1);%Assigning distribution of Dglut to the glutamate molecules.
Sqrt_distance=sqrt(2*dt*D);

%Glutamate diffusion coordinates.
X=zeros(Ntrns,Ntme);%x-coordinate of glutamate molecules.
Y=zeros(Ntrns,Ntme);%y-coordinate of glutamate molecules.
Z=zeros(Ntrns,Ntme);%z-coordinate of glutamate molecules.
R=zeros(Ntrns,Ntme);%radial coordinate of glutamate molecules.

%Release site in terms of radial distance from center.
Dist_release=0;%nm away from the center of cleft disc along the radius.
R(:,1)=Dist_release;
X(:,1)=Dist_release*cos(0);
Y(:,1)=Dist_release*sin(0); 

for i=1:Ntme-1
    
    %Normal random number generation.
    Xrandom=randn(Ntrns,1);
    Yrandom=randn(Ntrns,1);
    Zrandom=randn(Ntrns,1);
    
    %Diffusion step.
    X(:,i+1)=X(:,i)+(Sqrt_distance.*Xrandom);
    Y(:,i+1)=Y(:,i)+(Sqrt_distance.*Yrandom);
    Z(:,i+1)=Z(:,i)+(Sqrt_distance.*Yrandom);
    R(:,i+1)=sqrt((X(:,i+1).^2)+(Y(:,i+1).^2));
    
    %Checking for absorption boundary.
    Index1=find(R(:,i+1)>=Rabs);
    X(Index1,i+1)=10000;Y(Index1,i+1)=10000;Z(Index1,i+1)=0;
    R(Index1,i+1)=sqrt((X(Index1,i+1).^2)+(Y(Index1,i+1).^2));
    
    %Checking for upper reflecting boundary.
    Index1=find(R(:,i+1)<Rabs);
    Index2=find(Z(:,i+1)>0);
    Index3=intersect(Index1,Index2);
    Z(Index3,i+1)=-Z(Index3,i+1);
    
    %Checking for lower reflecting boundary.
    Index2=find(Z(:,i+1)<-h);
    Index3=intersect(Index1,Index2);
    Z(Index3,i+1)=-(2*h)-Z(Index3,i+1);
end


Glu_count=sum(double(R<=Rpsd).*double(Z<=eps-h),1);%recording the temporal profile of glutamate population within the prescribed
                          %local volume of radius Rpsd=100nm and height eps=5nm
                          %above the PSD.
%Result:                          
Glu_conc_Molar=10*(Glu_count/(pi*5*Rpsd^2))/6.023;%converting the temporal profile of glutamate population
                                                 %into the molar concentration.
