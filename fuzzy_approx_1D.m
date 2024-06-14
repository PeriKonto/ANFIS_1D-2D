% ====== Fuzzy Approximator ======
% Given an initial mesh of points it can use fuzzy rules 
% to approximate any surface (3-D) 
%
% OUTPUT = fuzzy_approx(X_X, Y_Y, Fp, N, infeng, mftype);
% X_X      Parameters on X-axis
% Y_Y      Parameters on Y-axis
% Fp       Approximated nonlinearity
% N        Number of MFs on both axis [3 4 5 7 13]
% infeng   Inference engine (1:Product 2:Minimum 3:Zadeh)
% mftype   MF Type (1:Triangular 2:Gauss 3:Trapezoidal 4:Bell 5:Sigmoid)
%
% OUTPUT   Approximated fuzzy surface
% BETA version 20/04/04
% MF number must change
% Peri Kontoroupis
load jun10_nonlinear
plot(states(:,1),par(:,1),'.')

%function OUTPUT = fuzzy_approx(X_X, Y_Y, Fp, N, infeng, mftype);
X_X=states(:,1);
Fp=par(:,1);
N=5;
infeng=1;
mftype=3;

%Start of the modification

n=1; %default value number of inputs 

%r=1:1:n;

%%N=[5 5];   % number of Mfs on both axis %3 4, 5, 7, 13 possible number of MFs to use

%M=prod(N); % No of fuzzy rules
M=N; % this is the modified  number of rules 

%X_X=XXp(1,:); Y_Y=YYp(:,1)';

szx=length(X_X); %szy=length(Y_Y);

%Place MFs on optional places
MF_x=N(1); %MF_y=N(2);

spacing_x=round((szx-2)/(MF_x-1)); %decide the spacing you use to put on the MFs 
%spacing_y=round((szy-2)/(MF_y-1));

%predeclare
%%Y_Y_mf=zeros(MF_y,1);
X_X_mf=zeros(MF_x,1);

%take the first and last node to contain a fuzzy rule
X_X_mf(1) = X_X(1);
%Y_Y_mf(1) = Y_Y(1);

for i = 2:MF_x  
    X_X_mf(i) = X_X(1)+(i-1)*(X_X(end)-X_X(1))/(MF_x-1);    
%    Y_Y_mf(i) = Y_Y(1)+(i-1)*(Y_Y(end)-Y_Y(1))/(MF_y-1);     
end    

% %this needs checking.......... 
% I(1,1)=X_X(1); I(1,2)=X_X(szx); % akrotata in X (original)
% I(2,1)=Y_Y(1); I(2,2)=Y_Y(szy); % akrotata in Y (original)
I(1,1)=X_X_mf(1); I(1,2)=X_X_mf(end); % akrotata in X_mf
%I(2,1)=Y_Y_mf(1); I(2,2)=Y_Y_mf(end); % akrotata in Y_mf (reduced to compensate the MFs)
xstep=X_X_mf(2)-X_X_mf(1);            % here step is left unoptimised... just step used as input for x
%ystep=Y_Y_mf(2)-Y_Y_mf(1);            % similarly for y
%%Steps(r)=[xstep ystep];         % save them as steps ...


%for i=1:2,                      % for two inputs something is missing
%    Plithos=[szx szy];          % possible fuzzy sets to be used... initial matrix
%end;
Plithos=szx;

% variable Stepmber needs checking...
for i=1%:2,
    Stepmber(i)=(I(i,2)-I(i,1))/(N(i)-1);  % <---- this guides where to place tbe MF
    Iindex(i,1)=1;
    Iindex(i,2)=N(i);
end;

for i=1:length(X_X_mf)
    for k=1:length(X_X_mf)
        Y_Ys((i-1)*length(Y_Y_mf)+k) = Y_Y_mf(k);
        X_Xs((i-1)*length(X_X_mf)+k) = X_X_mf(i); 
    end
end

PLEGMA=[X_Xs' Y_Ys'];

for i=1:MF_x
    for j=1:MF_y
        Y1(i, j)=Fp((i-1)*spacing_x+1, (j-1)*spacing_y+1);
    end
end

clear Y

Y=Y1(:);

% figure(7); plot(PLEGMA(:,1),PLEGMA(:,2),'r.'); 
% 
% figure(8); surf(reshape(PLEGMA(:,1),length(X_X_mf), length(X_X_mf)), ...
%     reshape(PLEGMA(:,2),length(X_X_mf), length(X_X_mf)) , reshape(Y,length(X_X_mf), length(X_X_mf)) ); 
% hold; plot3(PLEGMA(:,1),PLEGMA(:,2),Y','k.');legend('Nodes'); title('Initial selected mesh');

H = [(I(1,2)-I(1,1))/(szx-1) (I(2,2)-I(2,1))/(szy-1)]; %this is the default spacing for the original mesh

handle = waitbar(0,'Progress 1/3 ...');
INPUT = dimenscomb(I,H,handle);
K=prod(Plithos);  % INPUT AND K must have the same size

%MF parameters (set locations for all MFs) 
if mftype==1   %triangular MF
    
    for i=1:n,
        J=1:1:N(i);
        bvar(i,J)=I(i,1)+(J-1)*Stepmber(i);
        avar(i,J)=bvar(i,J)-Stepmber(i);
        dvar(i,J)=bvar(i,J)+Stepmber(i);
        avar(i,1)=I(i,1);
        bvar(i,1)=avar(i,1);
        dvar(i,N(i))=I(i,2);  
        bvar(i,N(i))=dvar(i,N(i));  
    end; 
    
elseif mftype==2  %gaussian MF 
    sigma=0.1; %gaussian variance
    for i=1:n,
        J=1:1:N(i);
        centervar(i,J)=I(i,1)+(J-1)*Stepmber(i); 
        centervar(i,1)=I(i,1);
        centervar(i,N(i))=I(i,2);  
    end;
    
elseif mftype==3
    %div adjusts trapezoidal's upper base, for example for div=2 then c-b=Stepmber 
    div=4;
    for i=1:n,
        J=1:1:N(i);
        centervar(i,J)=I(i,1)+(J-1)*Stepmber(i);
        avar(i,J)=centervar(i,J)-Stepmber(i);
        bvar(i,J)=centervar(i,J)-Stepmber(i)/div;
        cvar(i,J)=centervar(i,J)+Stepmber(i)/div;
        dvar(i,J)=centervar(i,J)+Stepmber(i);
        avar(i,1)=I(i,1);
        bvar(i,1)=avar(i,1);
        cvar(i,N(i))=I(i,2);  
        dvar(i,N(i))=cvar(i,N(i));  
    end;
    
elseif mftype==4 %gbell mf
    beta=2.5;  % shape
    alpha=0.1; % variance
    for i=1:n,
        J=1:1:N(i);
        centervar(i,J)=I(i,1)+(J-1)*Stepmber(i); 
        centervar(i,1)=I(i,1);
        centervar(i,N(i))=I(i,2);  
    end;
    
elseif mftype==5 %sigmoid mf
    alpha=10; % shape
    for i=1:n,
        J=1:1:N(i);
        centervar(i,J)=I(i,1)+(J-1)*Stepmber(i); 
        centervar(i,1)=I(i,1);
        centervar(i,N(i))=I(i,2);  
    end;
        
end;

%until here it should ok
% check this....
H=Steps;
H(r)=1;
handle = waitbar(0,'Progress 2/3 ');
odigos=dimenscomb(Iindex,H,handle); %GENERATION OF MESH TO BE USED

%there is a mistake here somewhere ...spacing

handle = waitbar(0,'Progress 3/3 ');
for k=1:K,  % Node possible combinations
    for m=1:M; %FUZZY RULES MF_X * MF_Y
        w(k,m)=1;
        wmin=1;
        for i=1:n, %number of inputs
            
            if mftype==1 %triangular mf
                a=avar(i,odigos(m,i)); %left side
                b=bvar(i,odigos(m,i)); %top
                d=dvar(i,odigos(m,i)); %right side
                if a==b; array(1)=NaN;
                else array(1)=(INPUT(k,i)-a)/(b-a);
                end;
                array(2)=1;
                if b==d; array(3)=NaN;
                else array(3)=(d-INPUT(k,i))/(d-b);
                end;
                minar=min(array);  MF(i)=max(minar,0);
                
            elseif mftype==2, %gaussian mf
                center(i)=centervar(i,odigos(m,i));
                MF(i)=exp(-(INPUT(k,i)- center(i))^2/(2*sigma^2)); %mf type formula+
                
            elseif mftype==3, %trapezoidal mf
                a=avar(i,odigos(m,i));
                b=bvar(i,odigos(m,i));
                c=cvar(i,odigos(m,i));
                d=dvar(i,odigos(m,i));
                if a==b 
                    array(1)=Inf;
                else 
                    array(1)=(INPUT(k,i)-a)/(b-a);
                end;
                array(2)=1;
                if c==d
                    array(3)=Inf;
                else
                    array(3)=(d-INPUT(k,i))/(d-c);
                end;
                minar=min(array);      
                MF(i)=max(minar,0);         
                
            elseif mftype==4, %gbell mf
                center(i)=centervar(i,odigos(m,i));
                MF(i)= 1./((1+abs((INPUT(k,i)-center(i))/alpha))^(2*beta));
                
                
            elseif mftype==5, %sigmoid mf
                center(i)=centervar(i,odigos(m,i));
                MF(i)=1./(1 + exp(-alpha*(INPUT(k,i)-center(i))));     
            end;        
            
            %calculation of weights 
            if infeng==1     % product inference engine 
                w(k,m)=w(k,m)*MF(i);
            elseif infeng==2 % minimum inference engine
                w(k,m)=min(w(k,m),MF(i));
            elseif infeng==3 % Zadeh inference engine
                wmin=min(wmin,MF(i));
                w(k,m)=min(wmin,1-wmin);
            end;
        end;
    end;
    
    ar(k)=0;
    par(k)=0;
    
    for m=1:M,  %note nodes in sequence ... 
        ar(k)=ar(k)+Y(m)*w(k,m);  %algorithm _1: sum of product between weights and firing strength 
        par(k)=par(k)+w(k,m);     %algorithm _2: sum of weights i.e. sum(w,2);
    end;
    
    if par(k)==0
        par(k)=1e-50;
    end;
    
    %defuzzifier used: center of gravity 
    OUTPUT(k)=ar(k)/par(k);       %algorithm _3: defuzzification 
    waitbar(k/K,handle)
    
end;

close(handle);

OUTPUT=reshape(OUTPUT,Plithos(2),Plithos(1));














% 
% % ======= General Display =======
% 
% %Illustration of MFs***************
% spaceplot=0.001;
% stepmemberx=I(1,1):spaceplot:I(1,2);
% stepmembery=I(2,1):spaceplot:I(2,2);
% if mftype==1
%     for i=1:n,
%         for m=1:N(i);
%             a=avar(i,m);
%             b=bvar(i,m);
%             d=dvar(i,m);
%             params=[a b d];
%             if i==1,
%                 memberfunctx(m,:)=trimf(stepmemberx,params);
%             elseif i==2,
%                 memberfuncty(m,:)=trimf(stepmembery,params);
%             end;
%         end;
%     end;
% elseif mftype==2
%     for i=1:n,
%         for m=1:N(i);
%             center=centervar(i,m);
%             params=[sigma center];
%             if i==1,
%                 memberfunctx(m,:)=gaussmf(stepmemberx,params);
%             elseif i==2,
%                 memberfuncty(m,:)=gaussmf(stepmembery,params);
%             end;
%         end;
%     end;
% elseif mftype==3
%     for i=1:n,
%         for m=1:N(i);
%             a=avar(i,m);
%             b=bvar(i,m);
%             c=cvar(i,m);
%             d=dvar(i,m);
%             params=[a b c d];
%             if i==1,
%                 memberfunctx(m,:)=trapmf(stepmemberx,params);
%             elseif i==2,
%                 memberfuncty(m,:)=trapmf(stepmembery,params);
%             end;
%         end;
%     end;
%     
% elseif mftype==4
%     for i=1:n,
%         for m=1:N(i);
%             center=centervar(i,m);
%             params=[alpha beta center];
%             if i==1,
%                 memberfunctx(m,:)=gbellmf(stepmemberx,params);
%             elseif i==2,
%                 memberfuncty(m,:)=gbellmf(stepmembery,params);
%             end;
%         end;
%     end;
%     
% elseif mftype==5
%     for i=1:n,
%         for m=1:N(i);
%             center=centervar(i,m);
%             params=[alpha center];
%             if i==1,
%                 memberfunctx(m,:)=sigmf(stepmemberx,params);
%             elseif i==2,
%                 memberfuncty(m,:)=sigmf(stepmembery,params);
%             end;
%         end;
%     end;
%      
% end;
% 
% memberfunctx=memberfunctx';
% memberfuncty=memberfuncty';
% 
% figure(9);    
% subplot(1,2,1); plot(stepmemberx,memberfunctx);
% title('MFs on x'); axis tight; xlabel('x1');
% subplot(1,2,2); plot(stepmembery,memberfuncty);
% title('MFs on y'); axis tight; xlabel('x2');
% %*********************************
% 
% domain=[I(1,1) I(1,2) I(2,1) I(2,2)];
% 
% axis(domain);
% Xp=linspace(X_X(1), X_X(end), 13);
% Yp=linspace(Y_Y(1), Y_Y(end), 13);
% figure(10);surf(Xp,Yp,OUTPUT);
% xlabel('x1'); ylabel('x2'); title('Approximated parametic surface');
% ABSERROR=abs(Fp-OUTPUT);
% figure(11); mesh(Xp,Yp,ABSERROR); title('Absolute error between par & non-par surface')
% 
% for i=1:MF_x
%    for j=1:MF_y
%        Y1(i, j)=OUTPUT((i-1)*spacing_x+1, (j-1)*spacing_y+1);
%    end
% end
% 
% Y_final=Y1(:);
% figure(12); surf(reshape(PLEGMA(:,1),length(X_X_mf), length(X_X_mf)), ...
%    reshape(PLEGMA(:,2),length(X_X_mf), length(X_X_mf)) , reshape(Y,length(X_X_mf), length(X_X_mf)) ); 
% hold; plot3(PLEGMA(:,1),PLEGMA(:,2),Y_final','k.');title('Selected mesh after Fuzzy rule estimation');  
% 
% ABSERROR=abs( reshape(Y,length(X_X_mf), length(X_X_mf)) -  reshape(Y_final,length(X_X_mf), length(X_X_mf)) );
% 
% figure(13); mesh(reshape(PLEGMA(:,1),length(X_X_mf), length(X_X_mf)), ...
%    reshape(PLEGMA(:,2),length(X_X_mf), length(X_X_mf)) , ABSERROR ); 
% title('ABS error in selected mesh');


