function [green_T]=Green_three_layers(z_source,t_1,t_2,Nx,Nz,pix_size,k1,k2,k3)
   %Green function for three layers z is the orhogonal axis to the layer, r the radial distance from the z axis 

   
    Nx_big=round(sqrt(2)*Nx/2);  %Border enhancement with sqrt(2) factor because of the circular reconstuction
    z_1=0;
    z_2=t_2;  
    
    %%% z meshing, to avoid singularity, we replace with pix_size/(4*log(1+sqrt(2)))
    big_z=((1:Nz)*pix_size-pix_size/2)-t_1; 

    
    %%% r meshing, to avoid singularity, we replace with pix_size/(4*log(1+sqrt(2)))
    big_r=(1:round(Nx_big))*pix_size-pix_size/2;
    dh=100;
    
    
    delta=z_2-z_1;    
    z_s=(round(z_source/pix_size))*pix_size;           %Heatsource height
    half_green_T=zeros(round(Nx_big),Nz);
    green_T=zeros(2*Nx_big,Nz);                          
    BIG_R=zeros(Nx_big,Nz);
    
   %%%%%% Calculation of half of the 2D Green function
    for i=1:Nz 
        
        z=big_z(i);
        
        if z<=z_1      
            
            rexp = min([2*delta;2*z_2-z_s-z;z_s-z]);
            
           if abs(rexp)<pix_size  
               
               rexp=pix_size; 
               
           end
            
%          dh = 0.05/(20*rexp);
           
            N = 10*round(3/(rexp*dh));             
            big_h=(1:N)*dh;            
            G=zeros(round(Nx_big),N);            
            bigN(i)=N;            
            bigdh(i)=dh;            
            bigrexp(i)=rexp;
            
            for j=1:N
                
                h=big_h(j); 
                
                detA=(k2+k3)*(k2+k1)+(k2-k1)*(k3-k2)*exp(-2*h*delta);
                
                c1= (k2+k3)*exp(-h*(z_s-z));
                c2= (k2-k3)*exp(-h*(2*z_2-z_s-z));
                G(:,j)=(c1+c2)*besselj(0,h.*big_r)/(2*pi*detA)*dh;
                
            end
            
            half_green_T(:,i)=sum(G,2);
            
        end
        
        if (z>z_1 && z<=z_2)
            
            rexp = min([2*delta;2*delta+z_s-z;2*z_2-z_s-z;-2*z_1+z_s+z;2*delta-z_s+z]);
            
            if  abs(rexp)<pix_size, 
                rexp=pix_size;               
            end
            
            %dh = 0.05*1/(20*rexp);
            
            N = 10*round(3/(rexp*dh));
            big_h=(1:N)*dh;
            G=zeros(round(Nx_big),N);
            bigdh(i)=dh;
            bigN(i)=N;
            bigrexp(i)=rexp;
                      
            for j=1:N
                
                h=big_h(j);
                
                detA=(k2+k3)*(k2+k1)+(k2-k1)*(k3-k2)*exp(-2*h*delta);
                
                c1= (k2-k3)*(k2-k1)*exp(-h*(2*delta+z_s-z));
                c2= (k2-k3)*(k2+k1)*exp(-h*(2*z_2-z_s-z));
                G1=(c1+c2)*besselj(0,h.*big_r)/(4*pi*k2*detA)*dh;

                c3= (k2-k1)*(k2+k3)*exp(-h*(z_s-2*z_1+z));
                c4= (k2-k3)*(k2-k1)*exp(-h*(2*delta+z-z_s));
                G2=(c3+c4)*besselj(0,h.*big_r)/(4*pi*k2*detA)*dh;
                G(:,j)=G1+G2;
                
            end
            
            R=sqrt(big_r.^2+(z-z_s).^2);
            half_green_T(:,i)=1./(4*pi*k2*R)+sum(G,2)';
            
        end
        
        if (z>z_2)
            
            rexp = min([2*delta;z_s+z-2*z_1;z-z_s]);
            
            if  abs(rexp)<pix_size, 
                rexp=pix_size;               
            end
            
            %dh=0.05*1/(20*rexp);
            
            N=10*round(3/(rexp*dh));
            G=zeros(round(Nx_big),N);
            bigdh(i)=dh;
            bigN(i)=N;
            bigrexp(i)=rexp;
            
            
            for j=1:N
                
                h=big_h(j);               

                detA=(k2+k3)*(k2+k1)+(k2-k1)*(k3-k2)*exp(-2*h*delta);
                
                c1= (k2-k1)*exp(-h*(z_s+z-2*z_1));
                c2= (k2+k1)*exp(-h*(z-z_s));
                G(:,j)=(c1+c2)*besselj(0,h.*big_r)/(2*pi*detA)*dh;
                
            end
            
            half_green_T(:,i)=sum(G,2);
           
            
        end
        
        % Full 2D Green function calculation (flip symetric along z) 
        green_T(1:round(Nx_big),i)=flip(half_green_T(:,i),1);
        green_T(round(Nx_big)+1:2*Nx_big,i)=half_green_T(:,i);
        
        figure(1156)
        imagesc(imrotate(log(green_T),90))
        title('Green function construction (log scale)')
        colormap hot

        
    end
     
    % Full 3D Green function calculation (rotation symetric along z) 
    green_T_3D_big=zeros(Nx+1,Nx+1,Nz);
 
    for i=1:Nz
        
        img=rot_sym_1D_to_2D(green_T(:,i));
        img(Nx_big+1,Nx_big+1)=max(img(:));
        green_T_3D_big(:,:,i)=img(round(2*Nx_big-Nx)/2+1:2*Nx_big-round(2*Nx_big-Nx)/2+1,round(2*Nx_big-Nx)/2+1:2*Nx_big-round(2*Nx_big-Nx)/2+1);
        
     
    end  
    
    green_T=imresize(green_T_3D_big,Nx/(Nx+1));
    
     
end

