Codici vari per il progetto di ANEDP2 2011/2012

Autori: Elena Bulgarello, Ivan Fumagalli

Contiene solo files modificati da noi.

our_swesimodeladapt_gio_bc_predcorr.m è il main

griglie_varie.mat:
	p,e,t : points,edges,triangles della mesh square6x6.msh
	pp,ee,tt:		||		  square12x12.msh
	ppp,eee,ttt:		||		  square24x24.msh
	tspan,tstart,theta

quiverscale.m: scaricato da internet, serve per dare la scala delle frecce quando si usa quiver
temp.m:		file di appoggio per scriptini o prove

Modificati ma ancora in dropbox e non qui:
	assem_mi.m : introdotto qualche commento e provato a mettere il residuo a 0 o a 1 per agire sulla parte dello shock-capturing che dipende da esso
	assem_GLS_... : ~idem

NB le cartelle che finiscono con rescut0.1, a parte Ritter_rectang40x5_N192x24_1t5_100_rescut0.1, ce l'hanno perché non abbiamo più rimesso a posto assem_mi.m dopo averlo modificato per avere rescut0.1. Tuttavia, non è il caso di rifarle, dato che Ritter_rectang40x5_N192x24_1t5_100_rescut0.1 non mostra differenze significative dalla simulazione di riferimento

Da fine maggio: script salvati nelle cartelle delle simulazioni che funzionano

Per Ritter_val=f(c)_ogni5Corr: in assem_...adapt c'è un riassunto delle prove fatte dalla assem_...adapt_NEW fino alla 5bis, che è la "madre" della ogni5Corr	[07/06/13]

Per beach_parallel:
	1) a sinistra si discosta da soluzione esatta MA VA BENE, perché la soluzione esatta suppone di avere lunghezza infinita della diga in quella direzione, mentre noi abbiamo messo una parete (Dirichlet u.n=0)	[07/06/13]
	2) default: Discesa->slope = pi/36	Salita->slope = pi/36	[07/06/13]

cestino (since 07/06/13):
..._060613: son quelle che ho tentato di sistemare dopo che mi son dimenticato di come ho fatto la simulazione val=f(c)ogni2
..._NEW: quelle da cui sono state generate le ..._060613 e i file del commit 07/06/13
..._OLD: dovrebbero essere praticamente uguali a quelle del commit precedente al 07/06/13
..._Elena_Commit070613: hanno fatto funzionare tutto quello che al 07/06/13 ha la u=max imposta a mano (sostituite in tale commit da quelle che hanno anche val=f(c))