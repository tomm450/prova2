function [ matout ] = ellipse_corr( matin )
% - matin -> matrice coordinata in ingresso ([1 0 1]',[up LE dwn]')

[~,iLE] = min(matin(:,1));

up  = flipud(matin(1:iLE,:)); 

if min(up(:,1)) ~= 0
    error(''); 
end
if max(up(:,1)) ~= 1
    error(''); 
end
dwn = matin(iLE:end,:);
if min(dwn(:,1)) ~= 0
    error(''); 
end
if max(dwn(:,1)) ~= 1
    error(''); 
end


%figure(199)
%plot(up(:,1),up(:,2),'rx--')
%hold on
%plot(dwn(:,1),dwn(:,2),'rx--')


dwn(:,2) = -dwn(:,2);
% posso fare la stessa procedura per entrambi
[ up ]  = a_b_finder( up );
[ dwn ] = a_b_finder( dwn );

dwn(:,2) = -dwn(:,2);

%plot(up(:,1),up(:,2),'*')
%hold on
%plot(dwn(:,1),dwn(:,2),'*')

matout = [flipud(up(2:end,:)); dwn(1:end,:)];

end

