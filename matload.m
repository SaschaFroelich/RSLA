function [file] = matload()
   [filename,filepath]=uigetfile(['' pwd filesep '*.mat']);
   
   strcat(filepath,filename)
   
   loadfile=strcat(filepath,filename);
   
   file=load(['' loadfile '']);
   
end