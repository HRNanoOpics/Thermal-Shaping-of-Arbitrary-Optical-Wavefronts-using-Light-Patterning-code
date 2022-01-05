%Generation d'un profile HSD en forme de disque de rayon "rayon", Nx,Ny le
%nombre de pixel selon x,y de l'image, pix_size taille du pixel en mètre,
%Q_Ir chaleur emise par la source (on suppose distribution de source
%homogène) en W/m^2
%HSD= profil HSD
% Q_tot= chaleur totale emis par la source en W

function [HSD,Q_tot]=HSD_generation(Nx,Ny,radius,pix_size,Q_Ir,shape_option,period,l_vortex,dist_laguerre,maskcustom_path)

    power_per_pix=Q_Ir*pix_size^2;              %Puissance émise par pixel
    diameter=2*radius;
    

    HSD=zeros(Nx,Ny);
    
    if strcmp(shape_option,'Disk')
        for i=1:Nx

            for j=1:Ny

                    if (sqrt((i-Nx/2)^2+(j-Ny/2)^2)<radius/pix_size) 

                        HSD(i,j)=1;
                        

                    end            

            end

        end
        
    elseif strcmp(shape_option,'Square') 
        
            HSD(Nx/2-radius/pix_size:Nx/2+radius/pix_size-1,Ny/2-radius/pix_size:Ny/2+radius/pix_size-1)=1;
            
    elseif(strcmp(shape_option,'1D grating')) 
        
        pattern_sin_x=sin(2*pi*(1:round(diameter/pix_size))/(period/pix_size));
        pattern_square=(pattern_sin_x>0==1);          
        HSD(round((Nx-diameter/pix_size)/2):round((Nx-diameter/pix_size)/2)+round(diameter/pix_size)-1,round((Ny-diameter/pix_size)/2):round((Ny-diameter/pix_size)/2)+round(diameter/pix_size)-1)=meshgrid(pattern_square);
     
    elseif (strcmp(shape_option,'2D grating'))
        
        pattern_sin_x=sin(2*pi*(1:round(diameter/pix_size))/(period/pix_size));
        pattern_square=(pattern_sin_x>0==1);
        HSD(round((Nx-diameter/pix_size)/2):round((Nx-diameter/pix_size)/2)+round(diameter/pix_size)-1,round((Ny-diameter/pix_size)/2):round((Ny-diameter/pix_size)/2)+round(diameter/pix_size)-1)=pattern_square.*pattern_square';
        
        
    elseif (strcmp(shape_option,'Chessboard'))        

        check=checkerboard(round(period/pix_size/2),round(diameter/pix_size),round(diameter/pix_size))>0.5;
        HSD(round((Nx-diameter/pix_size)/2):round((Nx-diameter/pix_size)/2)+round(diameter/pix_size)-1,round((Ny-diameter/pix_size)/2):round((Ny-diameter/pix_size)/2)+round(diameter/pix_size)-1)=check(1:round(diameter/pix_size),1:round(diameter/pix_size));
       
        
    elseif (strcmp(shape_option,'Vortex'))
       
      
        x=-round(diameter/pix_size/2):1:round(diameter/pix_size/2)-1;
        y=x;
        [X,Y]=meshgrid(x,y);
        [TH,~]=cart2pol(X,Y);        
        vortex=wrapTo2Pi(TH*l_vortex)/2/pi;
        HSD(round((Nx-diameter/pix_size)/2):round((Nx-diameter/pix_size)/2)+round(diameter/pix_size)-1,round((Ny-diameter/pix_size)/2):round((Ny-diameter/pix_size)/2)+round(diameter/pix_size)-1)=vortex;
        
        
    elseif (strcmp(shape_option,'Gauss-Laguerre'))
        
        x=linspace(-dist_laguerre,dist_laguerre,round(diameter/pix_size));
        y=x;
        Z=0.6;
        [X,Y]=meshgrid(x,y);
        [TH,R]=cart2pol(X,Y); 

         n = abs(l_vortex)+0;	% radial index; n=|l|,|l|+2,|l|+4 ...
         D = sqrt(2);   % is a constant for normalization;

        % Analytical functions
        G = @(r,z) D./sqrt(1+z.^2).*exp(-r.^2./(1+z.^2)).*exp(-1i/4*(z.*r.^2)./(1+z.^2));
        A = @(r,z) (sqrt(2)*r./sqrt(1+z.^2)).^abs(l_vortex).*laguerreL((n-abs(l_vortex))/2,abs(l_vortex),2*r.^2./(1+z.^2));
        PHI = @(th) exp(1i*l_vortex*th);
        PSI = @(z) exp(-1i*(n+1)*atan(z));
        P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*PSI(z).*exp(-1i*t);


        % Compute profile for a seleted time 't':
        p1=P(TH,R,Z,0);
        Lvortex=angle(p1);
        Lvortex=(Lvortex-min(Lvortex(:)))/(max(Lvortex(:))-min(Lvortex(:)));
        HSD(round((Nx-diameter/pix_size)/2):round((Nx-diameter/pix_size)/2)+round(diameter/pix_size)-1,round((Ny-diameter/pix_size)/2):round((Ny-diameter/pix_size)/2)+round(diameter/pix_size)-1)=Lvortex;
        
         
    elseif (strcmp(shape_option,'Custom 2D mask'))
        
        filename=maskcustom_path;           
        im_HSD_RGB=(imread(filename));
        im_HSD=cast(reshape(im_HSD_RGB(:,:,1),[size(im_HSD_RGB,1),size(im_HSD_RGB,2)]),'double');
        maxim=max(im_HSD(:));
        minim=min(im_HSD(:));
        imnorm=(im_HSD-minim)/(maxim-minim);
        HSD=imresize(imbinarize(imnorm,0.5),[Nx,Ny]);            
        
    end
    
    HSD=HSD*power_per_pix; %Power normalization   
    Q_tot=sum(HSD(:));
end