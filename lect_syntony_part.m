function [data01,datatmp,premms]=lect_syntony_part(fid,nbr_ech,apasser,datatmp,premms)
% fid : identifiant du fichier à lire
% nbr_ech : nombre d'échantillons à lire
% à passer : nombre d'échantillons à ne pas utiliser en début de fichier
% datatmp : vide au 1er passage dans la fonction
% premms : fixer à 1 pour le 1er passage dans la fonction

if premms==1
    ms1int8=fread(fid,nbr_ech/8,'uint8');
    ms1bits=reshape((fliplr(dec2bin(ms1int8,8)-'0')).',1,[]).';
    ms2int8=fread(fid,nbr_ech/8,'uint8');
    ms2bits=reshape((fliplr(dec2bin(ms2int8,8)-'0')).',1,[]).';
    datams1=ms1bits(apasser+1:end); % les bits qu'on veut de la 1ere ms
    datams2=ms2bits(1:apasser); % les bits qu'on veut de la 2eme ms
    data01=[datams1;datams2]; % une ms de data décalée du nombre d'échantillons à passer
    datatmp=ms2bits(apasser+1:end); % les bits à conserver pour le prochain appel
    premms=0;
else
    % datatmp contient déjà les bits qu'on veut de la 1ere ms (déjà lue)
    ms2int8=fread(fid,nbr_ech/8,'uint8');
    ms2bits=reshape((fliplr(dec2bin(ms2int8,8)-'0')).',1,[]).';
    datams2=ms2bits(1:apasser); % les bits qu'on veut de la 2eme ms
    data01=[datatmp;datams2]; % une ms de data décalée du nombre d'échantillons à passer
    datatmp=ms2bits(apasser+1:end); % les bits à conserver pour le prochain appel
end