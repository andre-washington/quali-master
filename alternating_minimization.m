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
%% Inicialize the set of precoders V
V = complex( rand(3,1,3), rand(3,1,3));
%% Find interference subspace for each user
Sub = zeros(3,3,Rx);
phi = zeros(3,1,Rx);
 
 
for rep =0:5
%     kl = phi(:,:,1)
%     kz = V(:,:,1)
    matrix1 = zeros(3,3);
    for k =1:Rx
        
        for j=1:Tx
    %         transmissor=k
    %         receptor = j
            if(k == 1 || k == 2)
                 if( j~= 1)
                   matrix1 = matrix1 + H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)';  %subspace for user k
                 end
            end
            if(k == 3 || k == 4)
                 if( j~= 2)
                    matrix1 = matrix1 + H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)';  %subspace for user k
                 end
            end
            if(k == 5 || k == 6)
                 if( j~= 3)
                    matrix1 = matrix1 + H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)';  %subspace for user k
                 end
            end

        end
         
%         [eigvec,eigval] = eig( matrix1);
%         [dominant,idx] = sort(diag(real(eigvec)));      
%[row,col]=find(real(eigval) == smallest(1,1)); % least dominant eigenvalue
%[row,col]=(find(eigval == max(max(eigval)))); %dominant eigenvalue     
%         phi(:,:,k) = eigvec(:,col); % subspace for each Rx k 
        [U,S,W] = svd(matrix1);
        phi(:,:,k) =  U(:,1);
    end
   

    %% Find the precoders (S dominant eigenvectors)
%-    from dissertation
%     for j =1:Tx
%         matrix2 = zeros(3,3);
%         for k=1:Rx
%     %         transmissor=k
%     %         receptor = j
%             if(j==1)
%                 if( k ~= 1 && k ~=2)
% 
%                     matrix2 = matrix2 + H(:,:,j,k)*phi(:,:,k)*phi(:,:,k)'*H(:,:,j,k)'; %
%                 end
%             end
%              if(j==2)
%                 if( k ~= 3 && k ~=4)
% 
%                     matrix2 = matrix2 + H(:,:,j,k)*phi(:,:,k)*phi(:,:,k)'*H(:,:,j,k)'; %
%                 end
%              end
%              if(j==3)
%                 if( k ~= 5 && k ~=6)
% 
%                     matrix2 = matrix2 + H(:,:,j,k)*phi(:,:,k)*phi(:,:,k)'*H(:,:,j,k)'; %
%                 end
%             end
%         end
% 
%         [eigvec,eigval] = eig( matrix2);
%         [row,col]=(find(eigval == max(max(eigval)))); %dominant eigenvalue 
%         V(:,:,j) = eigvec(:,col); %dominant eigenvector for Tx j
%     end

%--   from paper
    matrix2 = zeros(3,3);
    for j =1:Tx
       
        for k=1:Rx
        %         transmissor=k
        %         receptor = j
            if(j==1)
                if( k ~= 1 && k ~=2)

                    matrix2 = matrix2 + H(:,:,j,k)'*(eye(3,3) - phi(:,:,k)*phi(:,:,k)')*H(:,:,j,k); %
                end
            end
            if(j==2)
                if( k ~= 3 && k ~=4)

                    matrix2 = matrix2 + H(:,:,j,k)'*(eye(3,3) - phi(:,:,k)*phi(:,:,k)')*H(:,:,j,k); %
                end
            end
            if(j==3)
                if( k ~= 5 && k ~=6)

                    matrix2 = matrix2 + H(:,:,j,k)'*(eye(3,3) - phi(:,:,k)*phi(:,:,k)')*H(:,:,j,k); %
                end
            end
        end

%         [eigvec,eigval] = eig( matrix2);
%         [row,col]=(find(eigval == max(max(eigval)))); %dominant eigenvalue 
%         V(:,:,j) = eigvec(:,col) %dominant eigenvector for Tx j
        [U,S,W] = svd(matrix2);
        V(:,:,j) = U(:,1);
    end
end
%% objective function
J = 0;
for k=1:Rx
    for j=1:Tx
        if(k == 1 || k == 2)
             if( j~= 1)
                J= J + trace(phi(:,:,k)'*H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)'*phi(:,:,k));
             end
        end
        if(k == 3 || k == 4)
             if( j~= 2)
                J= J + trace( phi(:,:,k)'*H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)'*phi(:,:,k));
             end
        end
        if(k == 5 || k == 6)
             if( j~= 3)
                J= J + trace( phi(:,:,k)'*H(:,:,j,k)*V(:,:,j)*V(:,:,j)'*H(:,:,j,k)'*phi(:,:,k));
             end
        end
    end
end
J

%% from paper
J = 0;
for k=1:Rx
    for j=1:Tx
        if(k == 1 || k == 2)
             if( j~= 1)
                J= J + norm(H(:,:,j,k)*V(:,:,j) - phi(:,:,k)*phi(:,:,k)'*H(:,:,j,k)* V(:,:,j), 'fro')^2;
             end
        end
        if(k == 3 || k == 4)
             if( j~= 2)
                J= J + norm(H(:,:,j,k)*V(:,:,j) - phi(:,:,k)*phi(:,:,k)'*H(:,:,j,k)* V(:,:,j), 'fro')^2;
             end
        end
        if(k == 5 || k == 6)
             if( j~= 3)
                J= J + norm(H(:,:,j,k)*V(:,:,j) - phi(:,:,k)*phi(:,:,k)'*H(:,:,j,k)* V(:,:,j), 'fro')^2;
             end
        end
    end
end
J

