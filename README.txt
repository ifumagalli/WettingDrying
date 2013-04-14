Codici vari per il progetto di ANEDP2 2011/2012

Autori: Elena Bulgarello, Ivan Fumagalli

Contiene solo i files modificati da noi.

our_swesimodeladapt_gio_bc_predcorr.m è il main

Ivan_our_... non serve a niente

griglie_varie.mat:
	p,e,t : points,edges,triangles della mesh square6x6.msh
	pp,ee,tt:		||		  square12x12.msh
	ppp,eee,ttt:		||		  square24x24.msh
	tspan,tstart,theta

quiverscale.m: scaricato da internet, serve per dare la scala delle frecce quando si usa quiver

Modificati ma ancora in dropbox e non qui:
	assem_mi.m : introdotto qualche commento e provato a mettere il residuo a 0 o a 1 per agire sulla parte dello shock-capturing che dipende da esso
	assem_mat_vect_gio_stab_i e ..._stab_adapt : idem, per spegnere tutta la stabilizzazione