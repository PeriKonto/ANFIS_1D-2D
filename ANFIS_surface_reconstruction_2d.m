%Peri 

clear all; close all;
epoch_n = 10;
mf_n = [6 6];

%its working ...

% ====== collect training data
point_n = 12; x = linspace(-10, 10, point_n); y = linspace(-10, 10, point_n); [xx, yy] = meshgrid(x, y);
tmp1 = sin(xx)./(xx); index = find(isnan(tmp1)==1); tmp1(index) = ones(size(index));
tmp2 = sin(xx.*yy)./(yy); index = find(isnan(tmp2)==1); tmp2(index) = ones(size(index));
zz = tmp1.*tmp2; trn_data = [xx(:) yy(:) zz(:)];

% ====== training options
ss = 0.1; ss_dec_rate = 0.9; ss_inc_rate = 1.1; mf_type = 'gbellmf';

% ====== generate the initial FIS 
in_fismat = genfis1(trn_data, mf_n, mf_type);

% ====== start training
[trn_out_fismat trn_error step_size] = ...
   anfis(trn_data, in_fismat, [epoch_n nan ss ss_dec_rate ss_inc_rate], ...
   [1,1,1,1]);

   % ====== compute ANFIS output 
z_hat = evalfis([xx(:) yy(:)], trn_out_fismat);

% ====== plot of training data
% genfig('training data'); blackbg;
subplot(221); surfl(xx, yy, zz); limit = [min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:)) ...
      min(zz(:)) max(zz(:))];
axis(limit); set(gca, 'box', 'on'); xlabel('X'); ylabel('Y'); title('Training data');
zz_hat = evalfis([xx(:) yy(:)], trn_out_fismat);
subplot(222); surfl(xx, yy, reshape(zz_hat, point_n, point_n));
axis(limit); set(gca, 'box', 'on'); xlabel('X'); ylabel('Y'); title('ANFIS Output');
subplot(223); plot(1:epoch_n, trn_error); xlabel('epoch number'); ylabel('root mean squared error');
title('error curve');
subplot(224); plot(1:epoch_n, step_size); xlabel('epoch number'); ylabel('step size'); title('step size curve');

% ====== plot MFs
% figH = genfig('MFs for SINC function training');
% blackbg;

% plot initial and final MFs on x and y
figure(2)
subplot(221); plotmf(in_fismat, 'input', 1);
subplot(222); plotmf(in_fismat, 'input', 2);
subplot(223); plotmf(trn_out_fismat, 'input', 1);
subplot(224); plotmf(trn_out_fismat, 'input', 2);
% delete(findobj(figH, 'type', 'text'));
subplot(221); title('Initial MFs on X');
subplot(222); title('Initial MFs on Y');
subplot(223); title('Final MFs on X');
subplot(224); title('Final MFs on Y');

figure(3);%blackbg;
subplot(221);subplot(221); surfl(xx, yy, zz); limit = [min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:)) ...
      min(zz(:)) max(zz(:))];
axis(limit); set(gca, 'box', 'on'); xlabel('X'); ylabel('Y'); title('Training data');
subplot(222)
pcolor(xx, yy, reshape(zz, point_n, point_n));
axis square; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');title('Plane view for training data');
shading interp;grid
subplot(223); surfl(xx, yy, reshape(zz_hat, point_n, point_n));
axis(limit); set(gca, 'box', 'on'); xlabel('X'); ylabel('Y'); title('ANFIS Output');
subplot(224)
pcolor(xx, yy, reshape(zz_hat, point_n, point_n));
axis square; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');title('Plane view for ANFIS data');
shading interp;grid

nu=2;
[in,out,rules]=fismat2InOut(trn_out_fismat,mf_n,nu);

%lets do it element by elemnt
for stiles=1:length(xx(:)),

    %calculation of membership functions
    for ii=1:nu 
        for i=1:mf_n(ii) 
            if ii==1
                m_A(:,i)=gbellmf(xx(stiles),in(i,:,ii));   			             %Calculation of mfs for A
            else
                m_B(:,i)=gbellmf(yy(stiles),in(i,:,ii)); 			             %Calculation of mfs for B      
            end
        end
    end
    
    m_A_index(stiles,:) = m_A;
    m_B_index(stiles,:) = m_B;

    index=0;
        for ik = 1 : mf_n(1)
            for ijk = 1: mf_n(2)
                index=index+1;
                weights(:,index) = m_A(:,ik).*m_B(:,ijk);
            end
        end
        
  weights_index(stiles,:) = weights;

         for i = 1 : rules
             f_x_y(i,:) = xx(stiles)*out(i,1) + yy(stiles)*out(i,2) + out(i,3);
         end
    
    f_x_y_index(stiles,:) = f_x_y';
    
    %Multiplication: weights with firing strengths
    w_f(stiles,:)=weights_index(stiles,:).*f_x_y_index(stiles,:);  
    
end

%sum of:   
weight_sum=sum(weights_index,2); %denominator
w_f_sum=sum(w_f,2); %numerator

%algorithm for reconstruction
recons_z=w_f_sum ./ weight_sum;

figure(4);
subplot(121); surfl(xx, yy, reshape(zz_hat, point_n, point_n));
axis(limit); set(gca, 'box', 'on'); xlabel('X'); ylabel('Y'); title('ANFIS Output');
subplot(122); surfl(xx, yy, reshape(recons_z, point_n, point_n));
axis(limit); set(gca, 'box', 'on'); xlabel('X'); ylabel('Y'); title('Reconstructed ANFIS Output');

