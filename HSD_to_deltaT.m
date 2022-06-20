function [delta_T_map]=HSD_to_deltaT(HSD,green_T)

   %Convolution 3D calculation with Fourier transform
   % Code used in article entitled "Thermal Shaping of Arbitrary Optical Wavefronts using Light Patterning" 
   % from Hadrien M.L. Robert, Martin Cicala and Marek Piliarik*, 
   %Institute of Photonics and Electronics of the Czech Academy of Sciences, Chaberská 1014/57, 18251 Prague, Czech Republic. *piliarik@ufe.cz  

   
    HSD_F=fftn(fftshift(HSD));
    Green_T_F=fftn(fftshift(green_T));
    delta_T_map=fftshift(ifftn(Green_T_F.*HSD_F));
    

end
