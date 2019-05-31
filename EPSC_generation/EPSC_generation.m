%% Random to Uniform data extraction:
%This extracts the data over uniform time vector from the data over non-uniform
%time vector obtained from Gillespie algorithm under the stochastic activation of
%AMPA receptors.

Nonuniform_t=C_time;%Feed the non-uniform time vector obtained from Gillespie algorithm.
Uniform_t=t_refer;%Feed the t_refer mentioned in the Gillespie algorithm.

Raw_data=C_Openstate;%Feed the non-uniform data vector obtained from Gillespie algorithm.
Extracted_data=zeros(1,length(Uniform_t));
Extracted_data(1)=Raw_data(1);

for ii=2:length(Uniform_t)
    Extracted_data(ii)=Raw_data1(find(Nonuniform_t<=Uniform_t(ii),1,'last'));
end

Open_state=Extracted_data;%This will be used as a template for the EPSC generation.

%% EPSC generation

Open_prob=Open_state/Nrec;%Nrec is the no. of AMPARs in the PSD i.e. 30.
s_ampa=Open_prob/0.2141;%a constant calculated on the basis of nS conductance attributed to the individual 
                        %receptors in the open-state and Nrec.
                            
%Membrane potential dynamics 
V=zeros(1,Ntme);
V(1)=-70;%mV

Vl=-70;%Resting membrane potential in mV.
Cm=0.9;%Membrane capacitance in nanofarad.
gm=25;%Membrane leak conductance in nanoSieman.
Vampa=0;%Ampa reversal potential in mV.
gampa=0.0650;%Conductance of AMPA in nanoSieman.

for i=1:Ntme-1
    V(i+1)=V(i)+dt*(1/Cm)*((-1*gm*(V(i)-Vl))+(-1*w_ampa*gampa*s_ampa(i)*(V(i)-Vampa)));%Only includes AMPA mediated activation.    
end

%Result:
Iampa=-1*(w_ampa*gampa*s_ampa).*(V-Vampa);%Ampa excitatory current in picoAmpere.
