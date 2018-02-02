%--------------------------------------------------------------------------
function Y = slopefilter(X, nfilt)

if nfilt==0, return; end
% try
X=X';
Y = zeros(size(X), class(X));
switch nfilt
    case 1
        Y(1:end-1,:) = diff(X);
    case 2
        Y(2:end-2,:) = 2*diff(X(2:end-1,:)) + (X(4:end,:)-X(1:end-3,:));
    case 3
        Y(3:end-3,:) = 4*diff(X(3:end-2,:)) + 2*(X(5:end-1,:)-X(2:end-4,:)) + (X(6:end,:)-X(1:end-5,:));
    case 4
        Y(4:end-4,:) = 8*diff(X(4:end-3,:)) + 4*(X(6:end-2,:)-X(3:end-5,:)) + 2*(X(7:end-1,:)-X(2:end-6,:)) + (X(8:end,:)-X(1:end-7,:));
    otherwise
        error('ndiff_: nDiff_filt is valid for 1-4');
end
Y=single(Y');
end %func