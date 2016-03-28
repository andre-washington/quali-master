% mimo multicast
% base stations: 3
% user stations: 6
% scheme: BS1 transmit to US1 and US2. BS2 transmit to US3 and US4. BS3 transmit to US5 and US6.
clc;clear; close; 
Tx = 3; Rx = 6; % number of Tx and Rx
Nt =3; Nr = Nt; %num antenas at Tx and Rx
Ch = 1; % num of channels to simulate

% prealocation for speed
H = zeros(Nt,Nr,Tx,Rx);
for i=1:Tx
    for j=1:Rx
    H(:,:,i,j) = sqrt(1/2)*(randn(Nr,Nt)+1j*randn(Nr,Nt)); % rayleight channel matrix 
    end
end
% --- IA Zero Forcing ----%
%% Inicialize v3
% v3 = eye(3,1); %eigen vectors
[vec, val]= eig( (H(:,:,3,4)\H(:,:,1,4)) * (H(:,:,1,6) \ H(:,:,2,6))  * (H(:,:,2,1)\ H(:,:,3,1)) );
ymin=zeros(1,6)
for i=1:3
    v3 = vec(:,i);
    v1 = (H(:,:,1,4)\H(:,:,3,4))*v3;
    v2 = (H(:,:,2,1)\H(:,:,3,1))*v3;

    %% find u1 so that: u1*H(:,:,2,1)*v2 = 0, and so on...
    u1 = rand(1)*(H(:,:,2,1)*v2)/norm((H(:,:,2,1)*v2)); 
    % cross( H(:,:,2,1)*v2,u1)
    u2 = rand(1)*(H(:,:,2,2)*v2)/norm((H(:,:,2,2)*v2)); 
    u3 = rand(1)*(H(:,:,1,3)*v1)/norm((H(:,:,1,3)*v1)); 
    u4 = rand(1)*(H(:,:,1,4)*v1)/norm((H(:,:,1,4)*v1)); 
    u5 = rand(1)*(H(:,:,1,5)*v1)/norm((H(:,:,1,5)*v1)); 
    u6 = rand(1)*(H(:,:,1,6)*v1)/norm((H(:,:,1,6)*v1)); 

    %% calculate the minimum
    y=zeros(3,6);
    y(:,1)=cross(u1,(H(:,:,1,1)*v1));
    y(:,2)=cross(u2,(H(:,:,1,2)*v1));
    y(:,3)=cross(u3,(H(:,:,2,3)*v2));
    y(:,4)=cross(u4,(H(:,:,2,4)*v2));
    y(:,5)=cross(u5,(H(:,:,3,5)*v3));
    y(:,6)=cross(u6,(H(:,:,3,6)*v3));
    
    for j=1:2
        temp= norm(y(:,j))^2
        if (ymin(:,j)==0)
            ymin(:,j) = temp;
        end
        if(ymin(:,j) >  temp)
            ymin(:,j) =  temp
        end
    end
end





