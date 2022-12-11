%% Belief Propagation-Soft Decision
clc;clear;
c=[1 1 0 1 1 1];%LDPC code received by the receiver
p=0.8;%Set the initial prior probability
H=[1 1 0 1 0 0; 1 0 1 0 1 0; 0 1 1 0 0 1];%Parity check matrix
[x,y]=size(H);
pro1=log((1-p)/p);%The log likelihood ratio of a code point with the value 1
pro0=log(p/(1-p));%The log likelihood ratio of a code point with the value 0
clog=c;
clog(find(clog==1))=pro1;
clog(find(clog==0))=pro0;
M=zeros(x,y);
z=1;
n=1;%Number of iterations
%The H matrix is converted to the M variable matrix according to the log likelihood ratio
for m=1:length(c)
    if c(m)==1
        M(:,m)=pro1*H(:,m);
    else
        M(:,m)=pro0*H(:,m);
    end
end
E=zeros(x,y);
%For a code point, it refers to the "opinions" of other relevant code points in the check equation
for i=1:x
    for j=1:y
        prod=1;
        for m=1:y
            if m~=j && M(i,m)~=0
                prod=prod*atan(M(i,m)/2);

            end
            if M(i,j)~=0
                E(i,j)=log((1+prod)/(1-prod));
            end
        end

    end
end%Each line of E represents the electrical level probability of each code point
%Adding up the elements of each column and adding the "initial logarithmic natural ratio" gives an exactelectrical level probability
Emax=sum(E,1)+clog;
message=Emax;
%A bit decision (negative decision 1, positive decision 0) can be made based on the "posterior probability vector".
message(find(Emax<=0))=1;
message(find(Emax>=0))=0;
s=mod(H*message',2);%The output is evaluated against the check matrix H

if s==0
    z=0;
else
    n=n+1;%Proceed to the next iteration after the error judgment
    while(z)
        Estor=E;
        M=zeros(x,y);
        %The new M matrix is the log likelihood ratio combined with the information of the probability matrix E
        for i=1:x
            for j=1:y
                if(E(i,j)~=0)
                    M(i,j)=sum(E(:,j),1)-E(i,j)+clog(j);
                else
                    M(i,j)=0;
                end
            end
        end
        E=zeros(x,y);
        %Repeat to produce a new probability matrix E
        for i=1:x
            for j=1:y
                prod=1;
                for m=1:y
                    if m~=j && M(i,m)~=0
                        prod=prod*atan(M(i,m)/2);
                    end
                    if M(i,j)~=0
                        E(i,j)=log((1+prod)/(1-prod));
                    end
                end

            end
        end
        Emax=sum(E,1)+clog;
        message=Emax;
        message(find(Emax<=0))=1; 
        message(find(Emax>=0))=0;
        s=mod(H*message',2);%The same as the previous judgment
        if s==0
            z=0;
        else
            z=1;n=n+1;
        end
    end
end
message



                

    