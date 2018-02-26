function op=hthresh(X,T);

%A function to perform hard thresholding on a 
%given an input vector X with a given threshold T
% H=hthresh(X,T);

    ind=find(abs(X)<=T);
    X(ind)=0;
    op=X;
    