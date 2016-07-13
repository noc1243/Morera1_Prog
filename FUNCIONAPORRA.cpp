else if (tipo=='M')
     {
          g =  FATORDC;
          netlist[i].gm = 0.0;
          netlist[i].gds = 0.0;
          netlist[i].gmb = 0.0;
          io = 0.0;

          double vds = vAtual[netlist[i].td]-vAtual[netlist[i].ts];

          if ((vds < 0 && netlist[i].pnmos == nmos) || (vds > 0 && netlist[i].pnmos == pmos))
          {
           	  int aux  = netlist[i].td;
           	  netlist[i].td = netlist[i].ts;
           	  netlist[i].ts = aux;
          }

		  #ifdef DEBUG
          	  //printf ("vd %f vs %f vb %f vg %f\n", vAtual[netlist[i].td],vAtual[netlist[i].ts], vAtual[netlist[i].tb], vAtual[netlist[i].tg]);
		  //getch();
		  #endif

          double vgs = vAtual[netlist[i].tg]-vAtual[netlist[i].ts];
          vds = vAtual[netlist[i].td]-vAtual[netlist[i].ts];
          double vbs = vAtual[netlist[i].tb]-vAtual[netlist[i].ts];

          vgs *= (netlist[i].pnmos==nmos?1.0:-1.0);
          vds *= (netlist[i].pnmos==nmos?1.0:-1.0);
          vbs *= (netlist[i].pnmos==nmos?1.0:-1.0);

          vbs = (vbs>netlist[i].phi/2.0?netlist[i].phi/2.0:vbs);



          double vt = netlist[i].vt + netlist[i].gama * (sqrt((netlist[i].phi - vbs)) - sqrt((netlist[i].phi)));
          //double vt = netlist[i].vt + netlist[i].gama * (sqrt((netlist[i].phi - vbs)) - sqrt((netlist[i].phi)));

          #ifdef DEBUG
          	  printf ("vt: %f  vgs: %f vds: %f vs: %f ts %d\n", vt,vgs,vds, vAtual[netlist[i].ts], netlist[i].ts);
		  #endif


          if (vgs < vt)
          {
               netlist[i].transistorOp = corte;
			   #ifdef DEBUG
               	   printf ("corte\n");
			   #endif
          }
          else if (vds < vgs - vt)
          {
        	   netlist[i].transistorOp = ohmica;
			   #ifdef DEBUG
               	   printf ("ohmica\n");
			   #endif
          }
          else
          {
        	  netlist[i].transistorOp = saturacao;
			  #ifdef DEBUG
              	  printf ("saturacao\n");
			  #endif
          }

          if (netlist[i].transistorOp == saturacao)
          {
        	  netlist[i].gm = netlist[i].k * (netlist[i].w/netlist[i].l)* (2.0*(vgs - vt)) * (1 + netlist[i].lambda* vds);
        	  netlist[i].gds = netlist[i].k * (netlist[i].w/netlist[i].l)* pow ((vgs - vt),2.0) * netlist[i].lambda;
            io = netlist[i].k * (netlist[i].w/netlist[i].l) * pow((vgs - vt),2.0) * (1.0 + netlist[i].lambda * vds) - (netlist[i].gm * vgs) - (netlist[i].gds * vds);
          }
          else if (netlist[i].transistorOp == ohmica)
          {
        	  netlist[i].gm = netlist[i].k * (netlist[i].w/netlist[i].l)*(2.0* vds)*(1+ netlist[i].lambda* vds);
        	  netlist[i].gds = netlist[i].k * (netlist[i].w/netlist[i].l) * (2.0*(vgs - vt) - 2 * vds + 4.0* netlist[i].lambda * (vgs - vt) * vds - 3.0*netlist[i].lambda * pow (vds,2.0));
            io = netlist[i].k * (netlist[i].w/netlist[i].l) * (2.0* (vgs - vt)*vds - pow (vds,2.0)) - (netlist[i].gm * vgs) - (netlist[i].gds * vds);
          }

          if (gm != 0.0)
            cout << "gm: " << gm << endl;
          //Q?
          //netlist[i].gmb = (netlist[i].gm*netlist[i].gama)/(sqrt(fabs(netlist[i].phi - (vAtual[netlist[i].tb] -vAtual[netlist[i].ts]))));
          //netlist[i].gmb = (netlist[i].gm*netlist[i].gama)/(2*sqrt(fabs(netlist[i].phi - (vAtual[netlist[i].tb] -vAtual[netlist[i].ts]))));
          //netlist[i].gmb = (netlist[i].gm*netlist[i].gama)/(2*sqrt(fabs(netlist[i].phi - (vbs) ) ));

          //GM ESTA IGUAL A 0.0 E FUNCIONA WTF?
          netlist[i].gmb = (netlist[i].gm*netlist[i].gama)/(2*sqrt(fabs(netlist[i].phi - (vbs) ) ));

          //io-= 0.0;
          io-= (netlist[i].gmb * vbs);
          io*=(netlist[i].pnmos == pmos?-1.0:1.0);

          tvgs = vgs;
          tvds = vds;
          tvbs = vbs;
          //netlist[i].vds = vds;
          //netlist[i].vbs = vAtual[netlist[i].tb]-vAtual[netlist[i].ts];//vbs;
          //ti0 = netlist[i].k * (netlist[i].w/netlist[i].l) * (2.0* (vgs - vt)*vds - pow (vds,2.0)) * (1.0 + netlist[i].lambda * vds);
          ti0 = netlist[i].k * (netlist[i].w/netlist[i].l) * pow((vgs - vt),2.0) * (1.0 + netlist[i].lambda * vds);
	   	 #ifdef DEBUG
		 	   //printf ("gm %e gmb %e  gds %e io %e \n\n", gm, gmb, gds, io);
		 #endif

          Yn[netlist[i].td][netlist[i].tb]+=netlist[i].gmb;
          Yn[netlist[i].ts][netlist[i].ts]+=netlist[i].gmb;
          Yn[netlist[i].td][netlist[i].ts]-=netlist[i].gmb;
          Yn[netlist[i].ts][netlist[i].tb]-=netlist[i].gmb;

          Yn[netlist[i].td][netlist[i].tg]+=netlist[i].gm;
          Yn[netlist[i].ts][netlist[i].ts]+=netlist[i].gm;
          Yn[netlist[i].td][netlist[i].ts]-=netlist[i].gm;
          Yn[netlist[i].ts][netlist[i].tg]-=netlist[i].gm;

          Yn[netlist[i].td][netlist[i].td]+=netlist[i].gds;
          Yn[netlist[i].ts][netlist[i].ts]+=netlist[i].gds;
          Yn[netlist[i].td][netlist[i].ts]-=netlist[i].gds;
          Yn[netlist[i].ts][netlist[i].td]-=netlist[i].gds;


         Yn[netlist[i].td][netlist[i].td]+=g;
         Yn[netlist[i].tg][netlist[i].tg]+=g;
         Yn[netlist[i].td][netlist[i].tg]-=g;
         Yn[netlist[i].tg][netlist[i].td]-=g;

         Yn[netlist[i].ts][netlist[i].ts]+=g;
         Yn[netlist[i].tg][netlist[i].tg]+=g;
         Yn[netlist[i].ts][netlist[i].tg]-=g;
         Yn[netlist[i].tg][netlist[i].ts]-=g;

         Yn[netlist[i].tb][netlist[i].tb]+=g;
         Yn[netlist[i].tg][netlist[i].tg]+=g;
         Yn[netlist[i].tb][netlist[i].tg]-=g;
         Yn[netlist[i].tg][netlist[i].tb]-=g;

          Yn[netlist[i].td][nv+1]-=io;
          Yn[netlist[i].ts][nv+1]+=io;

     }