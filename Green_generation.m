function [green_T]=Green_generation(Nx,Ny,Nz,z_source,t_1,t_2,pix_size,k1,k2,k3,three_layers)
        
    %Green function calculation
    if three_layers==1
       
        green_T=Green_three_layers(z_source,t_1,t_2,Nx,Nz,pix_size,k1,k2,k3);    
      
         
    else

        %Mesh generation
        xx=((1:Nx)-round(Nx/2))*pix_size-pix_size/2;
        yy=((1:Ny)-round(Nx/2))*pix_size-pix_size/2;
        zz=((1:Nz)-round(Nz/2))*pix_size-pix_size/2;
        [XX,YY,ZZ]=meshgrid(xx,yy,zz);
        rho=sqrt(XX.^2+YY.^2+ZZ.^2);
        
       
        k0=(k1+k2)/2;
        green_T=1./(rho*4*pi*k0);  
       
    end

    
end