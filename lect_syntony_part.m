function [data01,datatmp,premms]=lect_syntony_part(fid,nbr_ech,apasser,datatmp,premms)
% fid : identifiant du fichier � lire
% nbr_ech : nombre d'�chantillons � lire
% � passer : nombre d'�chantillons � ne pas utiliser en d�but de fichier
% datatmp : vide au 1er passage dans la fonction
% premms : fixer � 1 pour le 1er passage dans la fonction

if premms==1
    ms1int8=fread(fid,nbr_ech/8,'uint8');
    ms1bits=reshape((fliplr(dec2bin(ms1int8,8)-'0')).',1,[]).';
    ms2int8=fread(fid,nbr_ech/8,'uint8');
    ms2bits=reshape((fliplr(dec2bin(ms2int8,8)-'0')).',1,[]).';
    datams1=ms1bits(apasser+1:end); % les bits qu'on veut de la 1ere ms
    datams2=ms2bits(1:apasser); % les bits qu'on veut de la 2eme ms
    data01=[datams1;datams2]; % une ms de data d�cal�e du nombre d'�chantillons � passer
    datatmp=ms2bits(apasser+1:end); % les bits � conserver pour le prochain appel
    premms=0;
else
    % datatmp contient d�j� les bits qu'on veut de la 1ere ms (d�j� lue)
    ms2int8=fread(fid,nbr_ech/8,'uint8');
    ms2bits=reshape((fliplr(dec2bin(ms2int8,8)-'0')).',1,[]).';
    datams2=ms2bits(1:apasser); % les bits qu'on veut de la 2eme ms
    data01=[datatmp;datams2]; % une ms de data d�cal�e du nombre d'�chantillons � passer
    datatmp=ms2bits(apasser+1:end); % les bits � conserver pour le prochain appel
end