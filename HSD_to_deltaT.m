function [delta_T_map]=HSD_to_deltaT(HSD,green_T)

   %Convolution 3D calculation with Fourier transform
   
    HSD_F=fftn(fftshift(HSD));
    Green_T_F=fftn(fftshift(green_T));
    delta_T_map=fftshift(ifftn(Green_T_F.*HSD_F));
    

end