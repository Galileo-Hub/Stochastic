% Function 'Stochastic_evolution_of_correlations' models a stochastic
% evolution of a chain of electrons while its particle number is
% continuously monitored (ie continuously measured). The evolution of
% two-point correlation functions and entanglement entropy is computed and
% analysed at each time step.

% Physical description of the system
% --------------------------------------
% We consider a chain of length L that is occupied by N electrons. An
% electron can occupy a site with probability p. The electrons move
% from an odd site to an even site to the right and backwards with
% the rate sh*(1-dh), and from an even site to an odd site to the right
% and backwards with rate sh*(1+dh). The system is inially prepare in a
% state in wich each even site is occupied with probability 1, and odd site
% is unoccupied. The system freely evolves in time, and while the electrons
% can rearrange themselves an observer (eg a measurement facility)
% monitores the number of electrons at each site at at each time step.

% Computational method and logic
% ------------------------------
% The chain of our model of electrons can be described by a Gaussian state, 
% which can fully represented by a two-point correlation function D. D is
% LxL matrix. The diagonal elements correspond to occupation propabilities
% p_j at each site j=1..L. The off-diagonal elements D_ij represent the
% correlations of particles at different sites i and j.
% In the absence of monitoring the transitioning rates of electrons from site to site
% can be represented by an LxL matrix h: h_ij =
% sh*(1+dh*(-1)^j)(delta_{i,j+1}+delta_{i,j-1}).

%   A continuous monitoring of the particle number with a rate gamma alters the time
% evolution of the two-point correlation function, which is described by
% equation:
%  dD = -i(h*D-D*h)*dt - gamma*(D-D_diag)*dt +(D*dW+dW*D-2*D*dW*D) 
% where the diagonal elements of the LxL matrix dW represent independent Wiener processes
% (Brownian motions). At each time step we draw dW_jj from a normal
% distribution with zero mean and 'gamma*dt' variance. 'D_diag' denotes
% diagonal LxL matrix of diagonal element of D.
%   We represent the correlation function D = U*(U') by a LxN matrix U that
% obeys isometry relation: (U')*U = NxN identity matrix.
% We imply the Trotterization method to compute the evolution of dD more
% effectively.

% Entry parameters:
% 'Navg' is the number of quantum trajectories we average over
% 'dimL' stands for L
% 'dimLA' is the length of a subregion A, for wich quantum entanglement
% entropy is evaluated
% 'dt' represents the time step
% 'tfin' stands for the final time, for which the modelling stops
% 'gamma' is the monitoring rate
% 'dh' , 'sh' are hopping rates introduced in the describtion of the model

% Output parameters:
% 'StavgA' is vector of the averaged entanglement entropy of subregion A evaluated at each time step
% 'StavgB', similarly to StavgA, represents entanglement entropy of a subregion B of length dimLA+1
% 'mDavg' gives two-point correlation function at time 'tfin' averaged over 'Navg' trajectories
% 'mD' gives two-point correlation function at time 'tfin' of the last computed single trajectories 
% 'mU' is the value of matrix U (see the description) at time 'tfin' of the last computed single trajectories 
% 't' is a vector of time steps of length 'dt' from 0 until 'tfin'

function [StavgA,StavgB,mDavg,mD,mU,t]=Stochastic_evolution_of_correlations(Navg,dimL,dimLA,dt,tfin,gamma,dh,sh)


dimN=fix(dimL/2); % for a half filled system N = L/2
if dimLA ~= dimL
    dimLB=dimLA+1;
else
    dimLB=dimL; %if dimLA = dimL then we consider dimLB = dimL too
end

% Functions Dan and Dan2 are for check of validity in the absence of
% stochasticity (ie gamma=0). These are exact results for correclations D.
%--------------------------------------------------------------------------
    % 'Dan(time)' returns correlation function D at time 'time'
    function mD1=Dan(time)                            
        mC=mV*diag(exp(-1i*vE*time))*((mV')*mU0);
        mD1=mC*(mC');
    end
    % 'Dan2(time)' returns correlation function D at time 'time' projected 
    % on the eigenmodes of 'h'
    function mDan2=Dan2(time)
        mCmu=diag(exp(-1i*vE*time))*((mV')*mU0);
        mDan2=mCmu*(mCmu');
    end

%clear all
tic

t = 0:dt:tfin;      % define the time vector 't' with time steps 'dt'
lt = length(t);

% The Hamiltonian matrix h is denoted by 'mh':
vh=1+dh*((-1).^[1:dimL-1]);mh=diag(vh,-1)+diag(vh,1);
mh=sh*mh;

% Eigevectors and eigenvalues of h are:
[mV,vE]=eig(mh);vE=diag(vE);

%Initial Conditions (BC):
mU0=zeros(dimL,dimN);
for i = 1:dimN
    mU0(2*i,i)=1;    % only even sites up to dimN are occupied with probablity 1
    %mU0(i,i)=1;     % alternative BC: all sites up to dimN are occupied with probablity 1
end
    %mU0=mV(:,1:dimN);   % alternative BC: only eigensome alternative initial condition

mD=mU0*(mU0');
mD1=Dan(0);dmD=mD1-mD;
mDA(:,:)=mD(1:dimLA,1:dimLA);mEigtA =eig(mDA);
mDB(:,:)=mD(1:dimLB,1:dimLB);mEigtB =eig(mDB);
UH=expm(-1i*mh*dt);

% figure(1000) shows D(0), D(tfin) and averaged D(tfin) in the 1st row
%                    and entanglement entropy S(t) in the 2nd row
figure(1000)
clf;
subplot(2,3,1);
hold off
imagesc(real(mD));
axis square;
colormap(parula);
colorbar;
hold on;
for i = 1:dimL
plot([0.5,dimL+0.5],[i-0.5,i-0.5],'k-');
plot([i-0.5,i-0.5],[0.5,dimL+0.5],'k-');
end
hold on
plot([0.5,dimLA+0.5],[dimLA+0.5,dimLA+0.5],'w-','LineWidth', 2);
plot([0.5,dimLA+0.5],[dimLA+0.5,dimLA+0.5],'r--','LineWidth', 2);
plot([dimLA+0.5,dimLA+0.5],[0.5,dimLA+0.5],'w-','LineWidth', 2);
plot([dimLA+0.5,dimLA+0.5],[0.5,dimLA+0.5],'r--','LineWidth', 2);
title('$|D_{ij}|$ at $t = 0$','interpreter','latex');

mDavg=0;StavgA=0;StavgB=0;ndavg=0;mEigHDavg=0;
for ii=1:Navg    % Generate and average over 'Navg' quantum trajectories
dW =sqrt(gamma)*sqrt(dt)*randn(dimL,lt-1); % Wiener processes dW
mD=mU0*(mU0');mU=mU0;

nd=diag(mD); % 'nd' is a vector of electron occupation probabilities at each site
%nd1=diag(Dan(0)); %% Checkpoint: in the limit of 'gamma=0' nd has to be
%equal to nd1, which is given by function 'Dan'.


mEH=sum(mh.*mD,'all'); % 'mEH' is the evarage Energy. We plot it below.

nd=diag(mD);
mDA(:,:)=mD(1:dimLA,1:dimLA);mEigtA =eig(mDA);
mDB(:,:)=mD(1:dimLB,1:dimLB);mEigtB =eig(mDB);
mEigHD = diag((mV')*mD*mV);
for i = 1:lt-1 % evolution of a single trajectory through each time step 'i*dt'
    vn=sum(mU.*conj(mU),2);
    UW=diag(exp(dW(:,i)+dt*gamma*(2*vn-1)));
    mU=UW*UH*mU;
    [mQ,mR]=qr(mU,0);mU=mQ;
mD=mU*(mU');
%mD1=Dan(i*dt); dmD=(mD1-mD);  %%Checkpoint:  at 'gamma=0' mD has to be
%equal to mD1, which is given by function 'Dan' at each time step 'i*dt'
mDA(:,:)=mD(1:dimLA,1:dimLA);mEigtA = [mEigtA eig(mDA)];
mDB(:,:)=mD(1:dimLB,1:dimLB);mEigtB = [mEigtB eig(mDB)];
    nd = [nd,diag(mD)];
%     nd1 = [nd1,diag(Dan(i*dt))];
%    mMD1=[mMD1,diag(Dan2(i*dt))];
    mEH = [mEH,sum(mh.*mD,'all')];
    mEigHD = [mEigHD,diag((mV')*mD*mV)];
%    sdt=round((mh.*mD)*10^10)*10^(-10);
end
mEigtA=abs(mEigtA)+1i*1e-22; mEigtB=abs(mEigtB)+1i*1e-22;
% The imaginary addition above is to help Matlab compute x*log(x) when x is zero,
% otherwise it returns NaN and does not plot the corresponding values

% 'StA' and 'StB' denote entanglement entropies of subregions A and B at
% each time step. These are vectors of length 'lt+1'
 StA=-sum(mEigtA(:,:).*log2(mEigtA(:,:))+(1-mEigtA(:,:)).*log2(1-mEigtA(:,:)),1);
 StB=-sum(mEigtB(:,:).*log2(mEigtB(:,:))+(1-mEigtB(:,:)).*log2(1-mEigtB(:,:)),1);
mDavg = mDavg+mD; StavgA = StavgA+StA; StavgB = StavgB+StB; ndavg = ndavg+nd; mEigHDavg = mEigHDavg+mEigHD;
end
mDavg=mDavg/Navg; StavgA=StavgA/Navg; StavgB=StavgB/Navg; ndavg=ndavg/Navg; mEigHDavg=mEigHDavg/Navg;

toc  % Return the total computational time.

%---------------------------------------------
% Below we plot figures to illustrate the results of the computation,
% including spectral analysis of the occupation probability.

    mEigtA=real(mEigtA);mEigtB=real(mEigtB);
    figure(1);clf;
    subplot(2,2,1);hold off;
        plot(t(2:end),mEigtA(:,2:end));
        set(gca,'XScale', 'log', 'YScale', 'linear','FontSize',12);
        xlabel('$t$','Interpreter','latex');ylabel('$\lambda^A_\alpha$','Interpreter','latex');
    subplot(2,2,2);hold off;
        plot(t,mEigtA);
        set(gca,'XScale', 'linear', 'YScale', 'linear','FontSize',12);
        xlabel('$t$','Interpreter','latex');ylabel('$\lambda^A_\alpha$','Interpreter','latex');
    subplot(2,2,3);hold off;
        plot(t(2:end),mEigtB(:,2:end));
        set(gca,'XScale', 'log', 'YScale', 'linear','FontSize',12);
        xlabel('$t$','Interpreter','latex');ylabel('$\lambda^B_\alpha$','Interpreter','latex');
    subplot(2,2,4);hold off;
        plot(t,mEigtB);
        set(gca,'XScale', 'linear', 'YScale', 'linear','FontSize',12);
        xlabel('$t$','Interpreter','latex');ylabel('$\lambda^B_\alpha$','Interpreter','latex');
        
    sticks=[1:(1+fix(dimL/10)):fix(dimL/2),flip(dimL:-(1+fix(dimL/10)):(fix(dimL/2)+1+fix(dimL/10)))];
    figure(10);clf;
    subplot(2,2,1);hold off;
         uimagesc(log10(t(2:end)),1:dimL,abs(nd(:,2:end)));
         title('(a)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'log', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$j$','Interpreter','latex');
    subplot(2,2,2);hold off;
         imagesc((t(2:end)),1:dimL,abs(nd(:,2:end)));
         title('(b)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'linear', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$j$','Interpreter','latex');
    subplot(2,2,3);hold off;
         uimagesc(log10(t(2:end)),1:dimL,abs(ndavg(:,2:end)));
         title('(c)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'log', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$j$','Interpreter','latex');
    subplot(2,2,4);hold off;
         imagesc((t(2:end)),1:dimL,abs(ndavg(:,2:end)));
         title('(d)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'linear', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$j$','Interpreter','latex');
         
    figure(20);clf;
    subplot(2,2,1);hold off;
         uimagesc(log10(t(2:end)),1:dimL,real(mEigHD(:,2:end)));
         title('(a)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'log', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$\nu$','Interpreter','latex');
    subplot(2,2,2);hold off;
         imagesc(t,1:dimL,real(mEigHD));
         title('(b)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'linear', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$\nu$','Interpreter','latex');
    subplot(2,2,3);hold off;
         uimagesc(log10(t(2:end)),1:dimL,real(mEigHDavg(:,2:end)));
         title('(c)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'log', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$\nu$','Interpreter','latex');
    subplot(2,2,4);hold off;
         imagesc(t,1:dimL,real(mEigHDavg));
         title('(d)');
         axis square;colormap(parula);colorbar;
         set(gca,'YTick',sticks);
         set(gca,'XScale', 'linear', 'YScale', 'linear','FontSize',12);
         xlabel('$t$','Interpreter','latex');ylabel('$\nu$','Interpreter','latex');
         
    figure(30);clf;
        subplot(3,1,1);hold off;plot(t,mEH);set(gca, 'XScale', 'log', 'YScale', 'linear');
        subplot(3,1,2);hold off;plot(t(2:end),mEigHD(:,2:end));set(gca, 'XScale', 'log', 'YScale', 'linear');
        subplot(3,1,3);hold off;plot(t(2:end),mEigHDavg(:,2:end));set(gca, 'XScale', 'log', 'YScale', 'linear');

    nnd2=diff(nd,1,1);     
    figure(40);clf;subplot(2,1,2);hold off;plot(t,nd);
                  subplot(2,1,1);plot(t(2:end),nd(:,2:end));set(gca, 'XScale', 'log', 'YScale', 'linear');
    figure(41);clf;scolor='-';
                  subplot(2,1,2);hold on;plot(t,nnd2,scolor);
                  xlabel('$t$','Interpreter','latex');ylabel('$\langle n\rangle_{j+1}-\langle n\rangle_{j}$','Interpreter','latex');
                  subplot(2,1,1);hold on;plot(t(2:end),nnd2(:,2:end),scolor);set(gca, 'XScale', 'log', 'YScale', 'linear');
                  xlabel('$t$','Interpreter','latex');ylabel('$\langle n\rangle_{j+1}-\langle n\rangle_{j}$','Interpreter','latex');

                  
                  
    Fs=1;%lt;%1/tfin;
    Fn = 2^nextpow2(lt);Fdim = 2;Fvar=length(nnd2(:,1));
    FY = fft(nnd2,Fn,Fdim);
    figure(80);scolor='-';clf;
    %for i=1:Fvar
        %subplot(Fvar,2,2*i-1);hold on;
        subplot(2,1,1);hold on;
        plot(0:(Fs/Fn):(Fs-Fs/Fn),abs(FY(:,1:end)),scolor)
        title(['Row ',num2str(i),' in the Frequency Domain'])
        set(gca, 'XScale', 'log', 'YScale', 'log');
        %subplot(Fvar,2,2*i);hold on;
        subplot(2,1,2);hold on;
        plot(0:(Fs/Fn):(Fs-Fs/Fn),abs(FY(:,1:end)),scolor)
        title(['Row ',num2str(i),' in the Frequency Domain'])
        set(gca, 'XScale', 'linear', 'YScale', 'linear');
    %end
                  
    figure(50);clf;subplot(3,1,1);hold off;loglog(t,sum(nd(dimLA:end,:),1));
                    subplot(3,1,2);hold off;plot(t,sum(nd(dimLA:end,:),1));
                    subplot(3,1,3);hold off;loglog(gamma*t,gamma*sum(nd(dimLA:end,:),1));
        

%------------------------------
figure(1000) % Upadating figure(1000), initiated earlier
subplot(2,3,[4,6]);
hold on
loglog(t,real(StavgA),'k.');
set(gca, 'XScale', 'log', 'YScale', 'log');
% hold on
% loglog(t,real(StavgB),'r.');
% set(gca, 'XScale', 'log', 'YScale', 'log');
subplot(2,3,2);
hold off
imagesc(abs(mD));
axis square;
colormap(parula);
colorbar;
hold on;
for i = 1:dimL
plot([0.5,dimL+0.5],[i-0.5,i-0.5],'k-');
plot([i-0.5,i-0.5],[0.5,dimL+0.5],'k-');
end
hold on
plot([0.5,dimLA+0.5],[dimLA+0.5,dimLA+0.5],'w-','LineWidth', 2);
plot([0.5,dimLA+0.5],[dimLA+0.5,dimLA+0.5],'r--','LineWidth', 2);
plot([dimLA+0.5,dimLA+0.5],[0.5,dimLA+0.5],'w-','LineWidth', 2);
plot([dimLA+0.5,dimLA+0.5],[0.5,dimLA+0.5],'r--','LineWidth', 2);
title('$|D_{ij}|$ at $t = t_{\rm fin}$','interpreter','latex');
subplot(2,3,3);
hold off
imagesc(abs(mDavg));
axis square;
colormap(parula);
colorbar;
hold on;
for i = 1:dimL
plot([0.5,dimL+0.5],[i-0.5,i-0.5],'k-');
plot([i-0.5,i-0.5],[0.5,dimL+0.5],'k-');
end
hold on
plot([0.5,dimLA+0.5],[dimLA+0.5,dimLA+0.5],'w-','LineWidth', 2);
plot([0.5,dimLA+0.5],[dimLA+0.5,dimLA+0.5],'r--','LineWidth', 2);
plot([dimLA+0.5,dimLA+0.5],[0.5,dimLA+0.5],'w-','LineWidth', 2);
plot([dimLA+0.5,dimLA+0.5],[0.5,dimLA+0.5],'r--','LineWidth', 2);
title('$|\bar{D}_{ij}|$ at $t = t_{\rm fin}$','interpreter','latex');


end
% End of the function Main function 'Stochastic_evolution_of_correlations'
% and an orderly return of [StavgA,StavgB,mDavg,mD,mU,t] variables