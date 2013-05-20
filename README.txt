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