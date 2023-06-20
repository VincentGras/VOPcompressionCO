function VOP = VOPcompressionCLU(filename, Q10g, compress)
% writes the Q-matrices used for SAR calculation inside a Siemens-specified file as an input to the VOP-compression.exe program 
% AA, May 2014, NB March 2015

Nch = size(Q10g, 1) ; % number of TX-channels
assert(ndims(Q10g) <= 3, 'Bad input Q10g');
assert(size(Q10g, 2) == Nch, 'Bad input Q10g');
nb_Q = size(Q10g, 3) ; % number of Q-matrices

fid=fopen(filename, 'w');
CoilStr = '8Tx8Rx_Head_Rapid'; % antenne Rapid
CoilID = '8Tx8Rx_Head_RAPID';
patientPosition = 'head first supine';

fwrite(fid, hex2dec('53313047') , 'uint32'); % magic number identifying file type
fwrite(fid, nb_Q , 'uint32');  % number of Q-matrices
fwrite(fid, Nch, 'uint32'); % number of TX-channels
%fwrite(fid, mass, 'float'); % exposed mass in kg
fprintf(fid, CoilStr); % coil name
fwrite(fid,zeros(128-length(CoilStr),1),'uint8');
fprintf(fid, CoilID);  % coil ID
fwrite(fid,zeros(128-length(CoilID),1),'uint8');
fprintf(fid,patientPosition); 
fwrite(fid,zeros(128-length(patientPosition),1),'uint8');

% regenerate SAR matrices in proper format
for iq = 1:nb_Q
    Smat = transpose(Q10g(:,:,iq));
    S2 = zeros(Nch, 2*Nch);
    S2(:,1:2:end) = real(Smat);
    S2(:,2:2:end) = imag(Smat);
    % Write VOP data
    for c=1:Nch
       fwrite(fid,S2(c,:),'float32');
    end
end

% Write scaling factors
scale = ones(nb_Q,1); % The VOPs matrices should already be in W/kg/V^2
fwrite(fid,scale,'float32');

fclose(fid);

% if (exist('VOPCompression') ~= 2)
%     addpath C:\Users\Public\VOP_management  
%     % prefer the code on \\sel\lrmn but should be also there : C:\Users\Public\MATLAB\SiemensResources\pTx\Step_2\VOP_management    
% end

VOPCompression (filename, filename, sprintf('%.3f', compress));
VOP = VOPFileAccess('readVOPs',filename);

end
