clear

num_threads_SL(10); % to use multiple cores

% Hamiltonian parameter
NC = 1;
J = 0.0;
T = 1e-8;

% NRG parameter
Lambda = 2.5;
N = ceil(-2*log(T/10)/log(Lambda)) + 10; %% Lambda^(-N/2) ~ E_N ~ T / 10 ==> N ~ -2 log(T/10) / log(Lambda)
Nkeep = 1000;
RhoV2in = [0.5;0.5]; % values outside of the 'ozin' grid are assumed to be zero

% % % % U(1) charge * SU(2) spin
% Define operators
[FF,ZF,SF,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',NC);
% request a spinful-fermion(FermionS) local Hilbert space with abelian U(1) and non-Abelian SU(2) formed from NC channels
% FF: annihilator (spinor), ZF: fermionic parity (for anticommutation),
% SF: spin operator, f^dagger_alpha sigma_alpha,beta f_beta., IF: identity and other metadata


[Fs,Zs,Ss,Is] = setItag('s00','op',FF,ZF,SF,IF.E);
% operators Fs, Zs etc are now set to act on the `s00` bath zeroth site space.
[FL,ZL,SL,IL] = setItag('L00','op',FF,ZF,SF,IF.E);
% these operators instead are now set to act on the `L00` or impurity space. Only the spin SF 
% is retagged into SL, because the full spinor is not needed for the Kondo model.

% project onto half-filled, spin-1/2 subspace
% IL = getsub(IL,find(IL.Q{1}(:,1) == 0)); % IL.Q{1} gets the first leg (of the tensor, (:, 1) selects the first column which is the U(1) charge
Z_L00 = ZL; % getsub(ZL,find(ZL.Q{1}(:,1) == 0));
A0 = getIdentity(IL,2,Is,2,'K00*',[1 3 2]);
FL_K00 = contract(A0, '!2*', {FL, '!1', {Zs, '!1', A0}}, [1 3 2]);
EK = getIdentity(A0,2);
Fs_K00 = contract(A0, '!2*', {Fs, '!1', A0}, [1 3 2]);
HSS = (2*J)*contract(SL,'!2*',Ss,[2 1 3 4]); % S*S interaction

ff = doZLD([-1 1],RhoV2in,Lambda,N,1);
V0 = ff{1}(1);

NF_K00 = quadOp(FL_K00,FL_K00,[]); % particle number operators
Ns_K00 = quadOp(Fs_K00,Fs_K00,[]); % particle number operators
HV = V0 * (contract(FL_K00, '!2*', Fs_K00));
HV = HV + HV'; % hybridisation term (V c^dagger_dsigma c_0sigma + hc)
specFuncOp = FL_K00; % Op = [(contract(Fs,'!1',HSS,[3 4 1 5 2]) - contract(HSS,'!13',Fs)); SL];
% cflag = [1;-1]; % brackets in the correlation functions are fermionic (+1) or bosonic (-1)
% zflag = [1;0]; % operators are fermionic (1) or not (0)

wmin = 1e-8;

fig = figure( ...
    'Color','white', ...
    'Units','pixels', ...
    'Position',[100 100 760 520]);

layout = tiledlayout(fig,1,1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax = nexttile(layout);
hold(ax,'on');

% for Ub = [-0, -20, -40]
for U = [1, 3, 5]
	wmax = U / 2 + 2;
	Ub = -0 * U;
	Himp = (U / 2) * NF_K00 * (NF_K00 - EK) - (U/2) * NF_K00;
	Hbath = (Ub / 2) * Ns_K00 * (Ns_K00 - EK) - (Ub/2) * Ns_K00;
	H0 = contract(A0,'!2*',{HSS,'!13',A0}) + Himp + Hbath + HV + 1e-40*getIdentity(A0,2); % impurity Hamiltonian

	nrgdata = NRG_SL([],H0,A0,Lambda,ff{1}(2:end),FF,ZF,'Nkeep',Nkeep);
	[Es,Qs] = plotE(nrgdata);
	nrgdata = getRhoFDM(nrgdata,T,'-v');
	[odisc,specFuncDisc,sigmak] = getAdisc(nrgdata,specFuncOp,specFuncOp,ZF,'Z_L00',Z_L00);

	[ocont,specFuncCont] = getAcont(odisc, specFuncDisc, sigmak, T/5);

	%% Plot figures
	ipos = ocont >= wmin & ocont <= wmax;
	plot(ax, ...
		ocont(ipos),specFuncCont(ipos), ...
		'-', ...
		'LineWidth',1.8, ...
		'DisplayName',sprintf('U=%g,Ub=%g', Ub/U));
end

hold(ax,'off');

set(ax, ...
    'FontSize',15, ...
    'LineWidth',1.1, ...
    'TickDir','out', ...
    'TickLabelInterpreter','latex', ...
    'XMinorTick','on', ...
    'YMinorTick','on', ...
    'Layer','top');

grid(ax,'on');
box(ax,'on');

ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;

xlabel(ax,'$|\omega|$', ...
    'Interpreter','latex', ...
    'FontSize',21);

ylabel(ax,'$A_f(\omega)$', ...
    'Interpreter','latex', ...
    'FontSize',21);

title(ax,sprintf( ...
    '$U=%g,\\quad V=%g,J=%g,\\quad T=%.1e$', ...
    U,V0,J,T), ...
    'Interpreter','latex', ...
    'FontSize',17);

legend(ax, ...
    'Location','best', ...
    'Interpreter','latex', ...
    'Box','off');

xlim(ax,[wmin wmax]);

drawnow;
ax.LooseInset = max(ax.TightInset,0.01);

exportgraphics( ...
    ax, ...
    sprintf('Aimp_log_U=%g_T=%.0e.pdf',U,T), ...
    'ContentType','vector', ...
    'BackgroundColor','white');
