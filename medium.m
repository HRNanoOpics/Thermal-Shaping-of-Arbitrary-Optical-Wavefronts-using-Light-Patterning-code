%for glass refractive index value, cf https://www.schott.com/d/advanced_optics/3794eded-edd2-461d-aec5-0a1d2dc9c523/1.0/schott_tie-19_temperature_coefficient_of_refractive_index_eng.pdf
% and schott-optical-glass-catalogue-excel-june-2018
function [kappa,cp,b_0,beta,n,Diffusivity]=medium(material,T_0)
     
    beta=zeros(1,4);
    b_0=0;
    b_1=0;
    b_2=0;
    b_3=0;
    b_4=0;
   %!!! only water has 4 order for the beta!!!
    if strcmp(material,'BK7')   
        b_0=1.5167;
        beta(1)=3e-6; 
        b_1=beta(1);
        kappa=1.13;
        cp=840;
        rho=2.42e3;
        Diffusivity=kappa/(cp*rho);
    elseif strcmp (material, 'air') %https://www.oharacorp.com/o2.html
        b_0=1;
        beta(1)=-1.36e-6; 
        b_1=beta(1);
        kappa=0.0262;
        cp=1000;
        rho=1;
        Diffusivity=kappa/(cp*rho);
    elseif strcmp(material,'diamond') %http://www.iiviinfrared.com/Optical-Materials/cvd-diamond_substrates.html
        b_0=1;
        beta(1)=6.5e-6; 
        b_1=beta(1);
        kappa=2200;
        cp=0.536;
        rho=3.51e3;
        Diffusivity=kappa/(cp*rho);
    elseif strcmp(material,'SF57')
        b_0=1.8464;
        beta(1)=3e-6;
        b_1=beta(1);
        kappa=1.38; 
        cp=0;
        Diffusivity=0;
    elseif strcmp(material,'P-SF68')
        b_0= 2.0048;
        beta(1)=24.1e-6;
        b_1=beta(1);
        kappa=0.65; 
        rho=2.43;
        cp=0.37;
        Diffusivity=kappa/(cp*rho);
    elseif strcmp(material,'water')
        
        %Based on polynomial fit taken from ref “Refractive index of water and its dependence on wavelength, temperature, and density,” 
        %Journal of physical and chemical reference data, vol. 14, no. 4, pp. 933–945, 1985.
        %n(T)=sum(b_j*T^j)
        
        b_0=1.3325;      
        b_1=-5.3864e-6;
        b_2=-2.0985e-6;
        b_3=6.8405e-9;
        b_4=-1.25e-11;        
       
        %beta(m)=1/m!*dn^m/dn^m
        
        beta(1)=b_1+2*b_2*T_0+3*b_3*T_0^2+4*b_4*T_0^3;
        beta(2)=b_2+3*b_3*T_0+6*b_4*T_0^2;
        beta(3)=b_3*T_0+4*b_4;
        beta(4)=b_4;
        T=0:200;
        n=b_0+b_1*T_0+b_2*T_0.^2+b_3*T_0.^3+b_4*T_0.^4;

        kappa=0.6;
        cp=4180;
        rho=1e3;
        Diffusivity=kappa/cp/rho;        
    elseif strcmp(material,'glycerol')
        b_0=1.473;
        beta(1) = -2.7e-4;
        b_1=beta(1);
        kappa=0.285;
        cp=2430;
        rho=1.261*10^3;
        Diffusivity=kappa/cp/rho;
        
        elseif strcmp(material,'PDMS') %Markos, C., Vlachos, K. & Kakarantzas, G. Thermo-optic effect of an index guiding photonic crystal fiber with elastomer inclusions. Proc.  SPIE7753,
        b_0=1.4000;
        beta(1) = -4.5e-4;
        b_1=beta(1);
        kappa=0.15;
        cp=1460;
        rho=0.97*10^3;
        Diffusivity=kappa/cp/rho;
     
        
    elseif strcmp(material,'sapphire')
        %source https://www.schott.com/d/advanced_optics/1bc1afd0-a532-4d31-a9b0-7d70324b6ba1/1.3/schott-sapphire-may-2013-eng.pdf
        b_0=1.764;
        beta(1) = 13e-6;
        b_1=beta(1);
        kappa=27.21;
        cp=77.87;
        rho=3.98e3;
        Diffusivity=kappa/(cp*rho);    
    elseif strcmp(material,'high_kappa100000')
        b_0=1.764; 
        beta(1) = 13e-6;
        b_1=beta(1);
        kappa=100000;
        cp=77.87;
        rho=3.98e3;
        Diffusivity=kappa/(cp*rho); 
   
    elseif strcmp(material,'none')
        beta(1) = 0;
        
        kappa=1e7;
        cp=0;
        Diffusivity=0;
    else
        error(['the material ''' material ''' is not in the database.\n'])
    end
    
       T=0:200;
       n=b_0+b_1*T+b_2*T.^2+b_3*T.^3+b_4*T.^4;

 
end