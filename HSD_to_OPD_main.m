
% 3D mapping of temperature, refractive index variation calculation and
% phase-shift induced by a 2D heat source. Effect of 2 or three medium layers
% Two configuration are possible: 
%   - 2 media, 1 and 2, the heat source lays at the interface of the two media
%   - 3 media, 1,2 and 3, the heat source lays in medium 2 and is usually a t the interfacte between media 1 and 2 (z_source=0) 
%
% We define a 3D coordinate system x,y,z with x,y the plan parallele to the
% interface layer. z=0 represents the interface between medium 1 and 2.
% Code used in article entitled "Thermal Shaping of Arbitrary Optical Wavefronts using Light Patterning" 
% from Hadrien M.L. Robert, Martin Cicala and Marek Piliarik*, 
%Institute of Photonics and Electronics of the Czech Academy of Sciences, Chaberská 1014/57, 18251 Prague, Czech Republic. *piliarik@ufe.cz  

% This code is available under GNU General Public Licence, version 3.0.


%% Simulation parameters

%%% Optical beam parameter and absorption
Incoming_power=500e-3;                  %Heating laser power (W)
lambda=488e-9;                          %Reference beam wavelength (m) (not the heating beam!)
eff_abs=0.55;                          %Gold nanoparticle absorption ratio at lambda heating

%%%  Heat source parameters
diameter=500e-6;                        %Heat source diameter (m)
radius=diameter/2;
heat_source_thickness=20e-6;            %Heat source thickness (if 3D heat source)        
 
shape_option='Chessboard';               % Shape can be: 'Disk', 'Square','1D grating', '2D grating', 'Chessboard', 'Vortex' or 'Gauss-Laguerre', 
                                        %to add other shapes, the function HSD_generation need to be uptated
period=200e-6;                        %period of the grating or chessboard in meter
l_vortex=1;                             %Vortex parameter             
dist_laguerre=0;                        %Gauss-Laguerre parameter 
maskcustom_path='custom_mask.tif';      %Image need to be in greyscale with the same meshing defined by the parameters

if heat_source_thickness>0
            heatsource3D=1;
        else
            heatsource3D=0;
end
        
%%%  Meshing
pix_size=5e-6;                          %Pixel size (m)

%%% Room temperature
T_0=20;                                 %Temperature (°C)

%%% Layers parameters
three_layers=1;                         %If =0, 2 medium, if=1, 3 medium

% materials can be: BK7, sapphire, PDMS, glycerol, air, diamond or water.
% to add more materials, the function material needs to be updated.

material_1='sapphire';                       %medium 1
material_2='PDMS';                  %medium 2
material_3=material_1;                  %medium 3 

t_1=340e-6;                             %Thickness medium 1
t_2=20e-6;                             %Thickness medium 2
t_3=340e-6;                             %Thickness medium 3=0 if only 2 medium

z_source=0e-6;                          %Heat source height from the interface medium 1/2


%% Total field size calulation and number of pixel along x,y,z

if three_layers==0

    t_3=0;

end

field_size=t_1+t_2+t_3;                 %Field_size en m

Nx=round(field_size/pix_size);          %Nb pixel along x   

if mod(Nx,2)>0

    Nx=Nx+1
end

Ny=Nx;                                  %Nb pixel along y 
Nz=Nx;                                  %Nb pixel along z

%% Thermal coefficient extraction

[k1,cp_1,b0_1,beta_1,n_1,Diffusivity_1]=medium(material_1,T_0);
[k2,cp_2,b0_2,beta_2,n_2,Diffusivity_2]=medium(material_2,T_0);
[k3,cp_3,b0_3,beta_3,n_3,Diffusivity_3]=medium(material_3,T_0);

%% Heating delivered power calculation

s_eff=pi*radius^2;                      %Illuminated surface (m^2)
Irradiance=Incoming_power/s_eff;        %Irradiance (W/m^2)
Q_uW_um2=Irradiance*10^-6*eff_abs       %delivered heat(µW/µm^2)
Q_Ir=Irradiance*eff_abs;                %delivered heat W/m2


%% HSD (Heat source density) map generation

[HSD,Q_tot]=HSD_generation(Nx,Ny,radius,pix_size,Q_Ir,shape_option,period,l_vortex,dist_laguerre,maskcustom_path); 
HSD_3D=zeros(Nx,Ny,round(t_2/pix_size));
figure(25)
imagesc(HSD)

%% 2D heat source case

if heatsource3D==0      
    
    % Green function calculation
    green_T=Green_generation(Nx,Ny,Nz,z_source,t_1,t_2,pix_size,k1,k2,k3,three_layers);

    % 3D map temperature calculation
    delta_T_map=HSD_to_deltaT(HSD,green_T);
    
else                    
    %% 3D HSD calculation
    for i=1:round(heat_source_thickness/pix_size)

        HSD_3D(:,:,i)=HSD/round(heat_source_thickness/pix_size);               

    end
    
    z_s=pix_size:pix_size:heat_source_thickness;   
    green_T_big=zeros(Nx,Ny,Nz,length(z_s));

    delta_T_map_big=zeros(Nx,Ny,Nz,length(z_s));
    TOM_layer=reshape(HSD_3D(:,round(Ny/2),:),[Nx,round(t_2/pix_size)]);   
    sym_HSD=sum(sum(abs(TOM_layer-fliplr(TOM_layer)))); %Check if it has symetry

    if sym_HSD>0 %Symetrize the meshing if the meshing is symetric to earn speed

            for i=1:length(z_s)
                
                % subGreen functions calculation
                green_T_big(:,:,:,i)=Green_generation(Nx,Ny,Nz,z_s(i),t_1,t_2,pix_size,k1,k2,k3,three_layers);
                % sub temperature map calculation
                delta_T_map_big(:,:,:,i)=HSD_to_deltaT(reshape(HSD_3D(:,:,i),[Nx,Ny]),reshape(green_T_big(:,:,:,i),[Nx,Ny,Nz]));
                figure(100)
                imagesc(imrotate(reshape(delta_T_map_big(Nx/2,:,:,i),[Ny,Nz]),90))
                colormap hot 

            end
            
    else

        if mod(length(z_s),2)==0

            for i=1:length(z_s)/2
                
                % subGreen functions calculation
                green_T_big(:,:,:,i)=Green_generation(Nx,Ny,Nz,z_s(i),t_1,t_2,pix_size,k1,k2,k3,three_layers);                       
                green_T_big(:,:,:,length(z_s)-i+1)=flip(green_T_big(:,:,:,i),3);
                
                % sub temperature map calculation
                delta_T_map_big(:,:,:,i)=HSD_to_deltaT(reshape(HSD_3D(:,:,i),[Nx,Ny]),reshape(green_T_big(:,:,:,i),[Nx,Ny,Nz]));
                delta_T_map_big(:,:,:,length(z_s)-i+1)=flip(delta_T_map_big(:,:,:,i),3);
                figure(100)
                imagesc(imrotate(reshape(delta_T_map_big(Nx/2,:,:,i),[Ny,Nz]),90))
                colormap hot
                figure(101)
                imagesc(imrotate(reshape(delta_T_map_big(Nx/2,:,:,length(z_s)-i+1),[Ny,Nz]),90))
                colormap hot

            end

        else

            for i=1:(length(z_s)-1)/2
                % subGreen functions calculation
                green_T_big(:,:,:,i)=Green_generation(Nx,Ny,Nz,z_s(i),t_1,t_2,pix_size,k1,k2,k3,three_layers);                       
                green_T_big(:,:,:,length(z_s)-i+1)=flip(green_T_big(:,:,:,i),3);
                
                % sub temperature map calculation
                delta_T_map_big(:,:,:,i)=HSD_to_deltaT(reshape(HSD_3D(:,:,i),[Nx,Ny]),reshape(green_T_big(:,:,:,i),[Nx,Ny,Nz]));
                delta_T_map_big(:,:,:,length(z_s)-i+1)=flip(delta_T_map_big(:,:,:,i),3);
                figure(100)
                imagesc(imrotate(reshape(delta_T_map_big(Nx/2,:,:,i),[Nx,Ny]),90))
                colormap hot
                figure(101)
                imagesc(imrotate(reshape(delta_T_map_big(Nx/2,:,:,length(z_s)-i+1),[Nx,Ny]),90))
                colormap hot

            end
                central_index=(length(z_s)+1)/2;
                green_T_big(:,:,:,central_index)=Green_generation(Nx,Ny,Nz,z_s(central_index),t_1,t_2,pix_size,k1,k2,k3,three_layers);
                delta_T_map_big(:,:,:,central_index)=HSD_to_deltaT(reshape(HSD_3D(:,:,i),[Nx,Ny]),reshape(green_T_big(:,:,:,central_index),[Nx,Ny,Nz]));
                figure(103)
                imagesc(imrotate(reshape(delta_T_map_big(Nx/2,:,:,central_index),[Nx,Ny]),90))
                colormap hot   


        end

    end
        %sum of the temperature map
        delta_T_map=sum(delta_T_map_big,4);
        %sum of the Green function
        green_T=sum(green_T_big,4);

end    



T_map=delta_T_map+T_0;                  %Room temperature addition

%% Refractive index and OPD variation calculation OPD=Optical path difference=sum(refractive index time thickness of the medium)

delta_n=zeros(Nx,Ny,Nz);

if three_layers==0

    t_1=field_size/2;
    t_2=t_1;

    delta_n_1=beta_1(1)*T_map(:,:,1:Nz/2)+beta_1(2)*T_map(:,:,1:Nz/2).^2+beta_1(3)*T_map(:,:,1:Nz/2).^3+beta_1(4)*T_map(:,:,1:Nz/2).^4;             % Delta n map in the medium 1
    delta_n_2=beta_2(1)*T_map(:,:,Nz/2+1:Nz)+beta_2(2)*T_map(:,:,Nz/2+1:Nz).^2+beta_2(3)*T_map(:,:,Nz/2+1:Nz).^3+beta_2(4)*T_map(:,:,Nz/2+1:Nz).^4; % Delta n map in the medium 2

    delta_n(:,:,1:Nz/2)=delta_n_1;
    delta_n(:,:,Nz/2+1:Nz)=delta_n_2;

    OPD_1=sum(delta_n(:,:,1:round(Nz/2)),3)*pix_size-(n_1(round(T_0))-b0_1)*t_1;        %OPD in the medium 1
    OPD_2=sum(delta_n(:,:,round(Nz/2)+1:Nz),3)*pix_size-(n_2(round(T_0))-b0_2)*t_2;     %OPD in the medium 2
    OPD_3=0;                                                                            %OPD in the medium 3=0 for a 2 layers system

else

    N_1=round(t_1/pix_size);        %Height of the medium 1/2 interface in pixel
    N_2=round(t_2/pix_size)+N_1;    %Height of the medium 2/2 interface in pixel  

    delta_n_1=beta_1(1)*T_map(:,:,1:N_1)+beta_1(2)*T_map(:,:,1:N_1).^2+beta_1(3)*T_map(:,:,1:N_1).^3+beta_1(4)*T_map(:,:,1:N_1).^4;                 % Delta n map in the medium 1
    delta_n_2=beta_2(1)*T_map(:,:,N_1+1:N_2)+beta_2(2)*T_map(:,:,N_1+1:N_2).^2+beta_2(3)*T_map(:,:,N_1+1:N_2).^3+beta_2(4)*T_map(:,:,N_1+1:N_2).^4; % Delta n map in the medium 2
    delta_n_3=beta_3(1)*T_map(:,:,N_2+1:Nz)+beta_3(2)*T_map(:,:,N_2+1:Nz).^2+beta_3(3)*T_map(:,:,N_2+1:Nz).^3+beta_3(4)*T_map(:,:,N_2+1:Nz).^4;     % Delta n map in the medium 3

    delta_n(:,:,1:N_1)=delta_n_1;    
    delta_n(:,:,N_1+1:N_2)=delta_n_2;
    delta_n(:,:,N_2+1:Nz)=delta_n_3;

    OPD_1=sum(delta_n(:,:,1:N_1),3)*pix_size-(n_1(round(T_0))-b0_1)*t_1;        %OPD in the medium 1
    OPD_2=sum(delta_n(:,:,N_1+1:N_2),3)*pix_size-(n_2(round(T_0))-b0_2)*t_2;    %OPD in the medium 2
    OPD_3=sum(delta_n(:,:,N_2+1:Nz),3)*pix_size-(n_3(round(T_0))-b0_3)*t_3;     %OPD in the medium 3

end

%% Phase-shift calculation
OPD_total=OPD_1+OPD_2+OPD_3;
phase=2*pi/lambda*OPD_total; 


%% Display result 

figure(1)
imagesc(T_map(:,:,round(t_1/pix_size)))
title 'Temperature profile along xy at the medium 1 and 2 interface'
colorbar
colormap hot
figure(2)
imagesc(imrotate(reshape(T_map(:,round(Ny/2),:),[Nx,Nz]),90))
title 'Temperature profile along xz'
colorbar
colormap hot
figure(3)
imagesc(OPD_total)
title 'ODP map'
colorbar
colormap hot
figure(4)
imagesc(phase)
title 'Phase map'
colorbar
colormap hot
figure(5)
x=(round(-Nx/2)+1:round(Nx/2))*pix_size*10^6;
plot(x,phase(:,round(Ny/2)))
title 'Phase profile'
xlabel 'x position (µm)'
ylabel 'Phase-shift (radian)'









