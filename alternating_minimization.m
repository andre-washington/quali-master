% mimo multicast
% base stations: 3
% user stations: 6
% scheme: BS1 transmit to US1 and US2. BS2 transmit to US3 and US4. BS3 transmit to US5 and US6.
clc;clear; close; 
format short;
Tx = 3; Rx = 6; % number of Tx and Rx
Nt =3; Nr = Nt; %num antenas at Tx and Rx
Ch = 1; % num of channels to simulate

% prealocation for speed
H = zeros(Nt,Nr,Tx,Rx);
% Channel Matrix H(Tx,Rx)
for i=1:Tx
    for j=1:Rx
    H(:,:,i,j) = sqrt(1/2)*(randn(Nr,Nt)+1j*randn(Nr,Nt)); % rayleight channel matrix 
    end
end
% --- Alternating Minimization Algorithm ----%
%% Inicialize V
V =complex( rand(3,1,3), rand(3,1,3))
%% Find interference subspace for each user
Sub = zeros(3,3,Rx);
phi = zeros(3,1,Rx);
for k =1:Rx
    for j=1:Tx
%         transmissor=k
%         receptor = j
        if( j~= k)
            mtrx1 = H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)';  %fazer somatorio
           [eigvec,eigval] = eig( mtrx1);
%            [row,col]=(find(eigval == max(max(eigval)))); % eigenvalue dominant
%            phi(:,:,k) = eigvec(:,col); % eigenvector dominant
         %retirar do laço for
          smallest = sort(diag(real(eigval)));
          [row,col]=find(real(eigval) == smallest(1,1)); % smallest eigenvalue
          phi(:,:,k) = eigvec(:,col); % subspace for each k user
          
        end
    end
end

%% Find the precoders (S dominant eigenvectors)
for k =1:Rx
    for j=1:Tx
%         transmissor=k
%         receptor = j
        if( j~= k)
            mtrx2 = H(:,:,j,k)*phi(:,:,j)*phi(:,:,j)'*H(:,:,j,k)'; %fazer somatorio
           %retirar do laço for
            [eigvec,eigval] = eig( mtrx2);
            [row,col]=(find(eigval == max(max(eigval)))); % eigenvalue dominant
            V(:,:,j) = eigvec(:,col); % eigenvector dominant
         
           
        end
    end
end


