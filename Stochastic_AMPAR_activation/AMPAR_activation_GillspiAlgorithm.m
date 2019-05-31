%% Stochastic activation of AMPA receptors based on the Gillespie algorithm. 
%%
%Defining time in seconds unit.
t_refer=linspace(0,0.001,10001);%so the simulation runs for 0.001 seconds= 1 microseconds. 
dt_min=(0.001-0)/(10001-1);%time-step in units of seconds.
    
%Now, load the Glu_conc_Molar; the temporal profile of glutamate concentration in molar
%obtained from the numerical simulation of glutamate diffusion. This temporal profile
%would serve as a template for the stochastic AMPAR activation.
Template=Glu_conc_Molar;

Nrec=30;%No. of AMPA receptors located in the PSD.

State_variables=zeros(8,1);%Populations of AMPA receptors in the different states of the 
                           %Milstein-Nicoll scheme of AMPAR activation.
State_variables(1)=Nrec;%Initailly, the entire AMPAR population is in the closed-state. 

t=0;%Initializing t for Gillespie algorithm.

while t(end)<=0.001
    Glu_conc=Template(find(t_refer<=t(end),1,'last'));
    AMPAR_activation_MNscheme;%A seperate script containing the kinetic parameters of transitions
                              %among the various states of the Milstein-Nicoll scheme. 
        
    Transition_matrix_temp=(repmat(State_variables(:,end),1,8)).*Transition_matrix;
    Index1=find(Transition_matrix_temp~=0);
    Vector1=Transition_matrix_temp(Index1);
    Vector1_sum=sum(Vector1);
    Propensity_vector=zeros(1,length(Vector1)+1);
    
    for kk=1:length(Vector1)
        Propensity_vector(kk+1)=sum(Vector1(1:kk));
    end
    
    Propensity_vector=Propensity_vector./Vector1_sum;
    
    random1=rand(1);
    tau=(1/Vector1_sum)*log(1/random1);
    if tau>dt_min
        tau=dt_min;
        t=[t t(end)+tau];
        
        State_variables(:,end+1)=State_variables(:,end);
    else
        t=[t t(end)+tau];
    
        random2=rand(1);
        Index2=min(find(Propensity_vector>random2));
        Index3=Index1((Index2)-1);
        Transition_matrix_logic=zeros(8,8);
        Transition_matrix_logic(Index3)=1;
    
        State_variables(:,end+1)=State_variables(:,end)-sum(Transition_matrix_logic,2)+sum(Transition_matrix_logic,1)';
    end
        
end

%Result:    
C_Openstate=State_variables(5,:)+State_variables(6,:);%the temporal profile of total open-state population.
C_time=t;%the time-vector.
