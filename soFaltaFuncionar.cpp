/*
Programa de demonstracao de analise nodal modificada
Por Antonio Carlos M. de Queiroz acmq@coe.ufrj.br
Versao 1.0 - 6/9/2000
Versao 1.0a - 8/9/2000 Ignora comentarios
Versao 1.0b - 15/9/2000 Aumentado Yn, retirado Js
Versao 1.0c - 19/2/2001 Mais comentarios
Versao 1.0d - 16/2/2003 Tratamento correto de nomes em minusculas
Versao 1.0e - 21/8/2008 Estampas sempre somam. Ignora a primeira linha
Versao 1.0f - 21/6/2009 Corrigidos limites de alocacao de Yn
Versao 1.0g - 15/10/2009 Le as linhas inteiras
Versao 1.0h - 18/6/2011 Estampas correspondendo a modelos
Versao 1.0i - 03/11/2013 Correcoes em *p e saida com sistema singular.
Versao 1.0j - 26/11/2015 Evita operacoes com zero.
*/

/*
Elementos aceitos e linhas do netlist:
Resistor:  R<nome> <no+> <no-> <resistencia>
VCCS:      G<nome> <io+> <io-> <vi+> <vi-> <transcondutancia>
VCVC:      E<nome> <vo+> <vo-> <vi+> <vi-> <ganho de tensao>
CCCS:      F<nome> <io+> <io-> <ii+> <ii-> <ganho de corrente>
CCVS:      H<nome> <vo+> <vo-> <ii+> <ii-> <transresistencia>
Fonte I:   I<nome> <io+> <io-> <corrente>
Fonte V:   V<nome> <vo+> <vo-> <tensao>
Amp. op.:  O<nome> <vo1> <vo2> <vi1> <vi2>
As fontes F e H tem o ramo de entrada em curto
O amplificador operacional ideal tem a saida suspensa
Os nos podem ser nomes
*/

#define versao "1.0j - 26/11/2015"
#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <complex>
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_TIPO 5
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-15
#define DEBUG
#define FATORDC 10e9 //fator multiplicativo para capacitores e indutores DC
#define MAX_IT 50 //maximo de iteracoes
#define REFVAL 1 //valor de referencia utilizado nos calculos de convergencias
#define MINCONV 2 //minimo de iteracoes para considerar como estavel a solucao
#define NCONV 5 // o numero de vezes que o algoritmo pode tentar calcular uma solucao
                //mesmo que esta nao esteja convergindo

#define MAX_ERRO 1e-9

enum ptoOperacao { corte, ohmica, saturacao};

enum mosType { nmos, pmos};

typedef struct elemento { /* Elemento do netlist */
  char nome[MAX_NOME];
  double valor, modulo, fase;
  double l,w,k,vt,lambda,gama,phi,ld, cOx;
  double cgb, cgs, cgd;
  mosType pnmos;
  int a,b,c,d,x,y,td,tg,ts,tb; // nos dos elementos, incluindo os do MOSFET
  ptoOperacao transistorOp;
} elemento;

elemento netlist[MAX_ELEM]; /* Netlist */

int
  ne, /* Elementos */
  nv, /* Variaveis */
  nn, /* Nos */
  i,j,k;

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],ntd [MAX_NOME],ntg [MAX_NOME],nts [MAX_NOME],ntb [MAX_NOME], tTipo[MAX_TIPO], nL[MAX_NOME], nW[MAX_NOME], nK[MAX_NOME], nVt[MAX_NOME], nLambda[MAX_NOME], nGama[MAX_NOME], nPhi[MAX_NOME], nLd[MAX_NOME],
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;
FILE *arquivo;

double
  g,
  Yn[MAX_NOS+1][MAX_NOS+2],
  //erroAtual = 0,
  vAtual[MAX_NOS+1],
  vProximo[MAX_NOS+1],
  gm = 0,
  gds = 0,
  io = 0;

ptoOperacao transistorOpAtual [MAX_NOS +1]; //TIRANDO ISSO

bool ptOperacao = true; //comeca o programa na analise do ponto de operacao
bool mantemModelo = true; //manter modelo do transistor
bool convergiu = false; //indica que a solucao convergiu
int contadorConv = 0; //conta quantas vezes os calculos deram erros menores
int iteracoes = 0; //conta o numero de iteracoes do alogoritmo
int vezNConvergiu = 0; //Conta a quantidade de vezes que o algoritmo nao convergiu

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
int resolversistema(void)
{
 // int i,j,l, a;
  double t, p;

  int l;
  for (int i=1; i<=nv; i++) {
    double t=0.0;
    int a=i;
    for (l=i; l<=nv; l++) {
      if (fabs(Yn[l][i])>fabs(t)) {
     a=l;
	t=Yn[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=nv+1; l++) {
	p=Yn[i][l];
	Yn[i][l]=Yn[a][l];
	Yn[a][l]=p;
      }
    }
    if (fabs(t)<TOLG) {
      printf("Sistema singular\n");
      printf ("%lf\n", (fabs(t) - TOLG)*1000000000);
      return 1;
    }
    for (j=nv+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
      Yn[i][j]/= t;
      p=Yn[i][j];
      if (p!=0)  /* Evita operacoes com zero */
        for (l=1; l<=nv; l++) {
	  if (l!=i)
	    Yn[l][j]-=Yn[l][i]*p;
        }
    }
  }
  return 0;
}

/* Rotina que conta os nos e atribui numeros a eles */
int numero(char *nome)
{
  int i,achou;

  i=0; achou=0;
  while (!achou && i<=nv)
    if (!(achou=!strcmp(nome,lista[i]))) i++;
  if (!achou) {
    if (nv==MAX_NOS) {
      printf("O programa so aceita ate %d nos\n",nv);
      exit(1);
    }
    nv++;
    strcpy(lista[nv],nome);
    return nv; /* novo no */
  }
  else {
    return i; /* no ja conhecido */
  }
  printf("Dentro da Convergencia: %d", convergiu );
}

//void calculaCAPs (double *Cgb, double *Cgs, double *Cgd, int modoOperacao);

bool controleConvergencia ( double vAtual[], double vProximo[], int iteracoes )
{
  double maxVal = 0;
  double tempVar;
  if (iteracoes < MAX_IT && (!iteracoes))
  {
     for (int cont = 0; cont <= nv; cont++)
     {
           //define o modo como o erro sera tratado
           if (fabs(vProximo[cont]) > REFVAL)
             tempVar = fabs( ( vProximo[cont] - vAtual[cont] ) / vProximo[cont]);

           if (vProximo[cont] < REFVAL)
             tempVar = fabs(vProximo[cont] - vAtual[cont]);

           //pega sempre o mesmo valor
           if (tempVar > maxVal)
             maxVal = tempVar;
          }

          //faz com o que as tensoes futuras sejam as atuais
          for (int cont = 0; cont <= nv; cont++)
          {
             vAtual[cont] = vProximo[cont];
          }

          //printf("Erro: %f", MAX_ERRO);
          printf("\nErro atual: %.10f\n", maxVal);
          if (maxVal <= MAX_ERRO)
          {
              mantemModelo = true;
              contadorConv++;
          }
          else if (maxVal >= MAX_ERRO)
          {
              vezNConvergiu++;
          }
          else if ((maxVal >= MAX_ERRO) && (vezNConvergiu >= NCONV))
          {
               mantemModelo = false;
               contadorConv = 0;
               //troca o modelo de transistor
               //e reinicia a contagem
               //se nao for pra manter o modelo
               iteracoes = 0;
          }

          if (contadorConv >= MINCONV)
          {
               convergiu = true;
               printf("Convergiu \n");
          }
               //Mudanca 21/06
              //erroAtual = maxVal;

                //parte onde voce coloca o transistor
          if (mantemModelo)
          {
             iteracoes++;
          }

           else if (iteracoes >= MAX_IT)
           {
             //troca o modelo do transistor
             //e reinicia a contagem
             //se estourar o limite de iteracoes
             iteracoes = 0;
             mantemModelo = false;
             contadorConv = 0;
           }
     }
}

void CalculaCapacitancias (elemento *netlist)
{
	if (netlist->transistorOp == corte)
	{
		netlist->cgb = netlist->cOx * netlist->w * netlist->l;
		netlist->cgs = 0;
		netlist->cgd = 0;
	}
	else if (netlist->transistorOp == ohmica)
	{
		netlist->cgb = 0;
		netlist->cgs = (1.0/2 * netlist->cOx * netlist->w * netlist->l) + (netlist->cOx * netlist->w * netlist->ld);
		netlist->cgd = (1.0/2 * netlist->cOx * netlist->w * netlist->l) + (netlist->cOx * netlist->w * netlist->ld);
	}
	else
	{
		netlist->cgb = 0;
		netlist->cgs = (2.0/3 * netlist->cOx * netlist->w * netlist->l) + (netlist->cOx * netlist->w * netlist->ld);
		netlist->cgd = netlist->cOx * netlist->w * netlist->ld;
	}
}

//lembrar de limitar a tensao no substrato para que nao seja maior que a no Gate


int main(void)
{
  printf("Programa demonstrativo de analise nodal modificada\n");
  printf("Por Antonio Carlos M. de Queiroz - acmq@coe.ufrj.br\n");
  printf("Versao %s\n",versao);
  denovo:
  /* Leitura do netlist */
  ne=0; nv=0; strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomearquivo);
  arquivo=fopen(nomearquivo,"r");
  if (arquivo==0) {
    printf("Arquivo %s inexistente\n",nomearquivo);
    goto denovo;
  }

  for (int cont = 0; cont <= MAX_NOS; cont++ )
  {
     vAtual[cont] = 0.1;
     vProximo[cont]=0;
     //transistorOpAtual [cont] = corte;
     //transistorOpProximo [cont] = saturacao;
  }
  vAtual [0] = 0;

  bool linear = true;
  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo)) {
    ne++; /* Nao usa o netlist[0] */
    if (ne>MAX_ELEM) {
      printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
      exit(1);
    }
    txt[0]=toupper(txt[0]);
    tipo=txt[0];
    sscanf(txt,"%10s",netlist[ne].nome);
    p=txt+strlen(netlist[ne].nome); /* Inicio dos parametros */
    /* O que e lido depende do tipo */
    if (tipo=='R' || tipo == 'C' || tipo == 'L') {
      if(ptOperacao)
      {
        sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
        printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
        netlist[ne].a=numero(na);
        netlist[ne].b=numero(nb);
      }
    }

    else if (tipo=='I' || tipo=='V') {
          if(ptOperacao)
          {
            sscanf(p,"%10s%10s%lg%lg%lg",na,nb,&netlist[ne].modulo,&netlist[ne].fase,&netlist[ne].valor);
            printf("%s %s %s %g %g %g\n",netlist[ne].nome,na,nb,netlist[ne].modulo, netlist[ne].fase, netlist[ne].valor);
            netlist[ne].a=numero(na);
            netlist[ne].b=numero(nb);
          }
    }

    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      if (ptOperacao)
      {
        sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[ne].valor);
        printf("%s %s %s %s %s %g\n",netlist[ne].nome,na,nb,nc,nd,netlist[ne].valor);
        netlist[ne].a=numero(na);
        netlist[ne].b=numero(nb);
        netlist[ne].c=numero(nc);
        netlist[ne].d=numero(nd);
      }
    }
    else if (tipo=='O') {
      sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
      printf("%s %s %s %s %s\n",netlist[ne].nome,na,nb,nc,nd);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo == 'K')
    	continue;

    else if (tipo == 'M')
    {
          sscanf (p, "%10s%10s%10s%10s%10s L=%lf W=%lf %lf%lf%lf%lf%lf%lf", ntd, ntg, nts, ntb, tTipo, &netlist[ne].l, &netlist[ne].w, &netlist[ne].k, &netlist[ne].vt, &netlist[ne].lambda, &netlist[ne].gama, &netlist[ne].phi, &netlist[ne].ld);
          printf("%s %s %s %s %s\n",netlist[ne].nome,ntd,ntg,nts,ntb);
          netlist[ne].pnmos = ((!strcmp (tTipo, "PMOS"))?nmos:pmos);
          netlist[ne].td=numero (ntd);
          netlist[ne].tg=numero (ntg);
          netlist[ne].ts=numero (nts);
          netlist[ne].tb=numero (ntb);
          netlist[ne].transistorOp = corte;
          double u = (netlist[ne].pnmos == nmos?0.0025:0.0067);
          netlist[ne].cOx = 2* netlist[ne].k/u;

          linear = false;
    }
    else if (tipo=='*') { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      ne--;
    }
    else {
      printf("Elemento desconhecido: %s\n",txt);
      getch();
      exit(1);
    }
  }
  fclose(arquivo);
  /* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */
  nn=nv;
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo=='L') {
      nv++;
      if (nv>MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
    else if (tipo=='H') {
      nv=nv+2;
      if (nv>MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv-1],"jx"); strcat(lista[nv-1],netlist[i].nome);
      netlist[i].x=nv-1;
      strcpy(lista[nv],"jy"); strcat(lista[nv],netlist[i].nome);
      netlist[i].y=nv;
    }
  }
  getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d -> %s\n",i,lista[i]);
  getch();
  printf("Netlist interno final\n");
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='R' || tipo=='I' || tipo=='V' || tipo=='C') {
      printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
    }
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      printf("%s %d %d %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,netlist[i].valor);
    }
    else if (tipo=='O') {
      printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
    if (tipo=='V' || tipo=='L' || tipo=='E' || tipo=='F' || tipo=='O')
      printf("Corrente jx: %d\n",netlist[i].x);
    else if (tipo=='H')
      printf("Correntes jx e jy: %d, %d\n",netlist[i].x,netlist[i].y);
  }
  getch();


  int vezes = 0;
//  for (int cont = 0; cont <= nv+1; cont++ )
//  {
//     vAtual[cont] = 0;
//     vProximo[cont]=0;
//     transistorOpAtual [cont] = saturacao;
//  }

  while (!convergiu)
  {
    /* Monta o sistema nodal modificado */
    printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
    getch();
    /* Zera sistema */
    for (int i=0; i<=nv; i++) {
      //inicializa os vetores utilizdos na analise de convergencia
      for (int j=0; j<=nv+1; j++)
      {
        Yn[i][j]=0;
      }
    }
    /* Monta estampas */
    for (int i=1; i<=ne; i++) {
      tipo=netlist[i].nome[0];
      printf ("\ntipo :%c, i:%d\n", tipo, i);
      if (tipo=='R') {
        g=1/netlist[i].valor;
        Yn[netlist[i].a][netlist[i].a]+=g;
        Yn[netlist[i].b][netlist[i].b]+=g;
        Yn[netlist[i].a][netlist[i].b]-=g;
        Yn[netlist[i].b][netlist[i].a]-=g;
      }
      else if (tipo=='G') {
        g=netlist[i].valor;
        Yn[netlist[i].a][netlist[i].c]+=g;
        Yn[netlist[i].b][netlist[i].d]+=g;
        Yn[netlist[i].a][netlist[i].d]-=g;
        Yn[netlist[i].b][netlist[i].c]-=g;
      }
      else if (tipo=='I') {
        g=netlist[i].valor;
        Yn[netlist[i].a][nv+1]-=g;
        Yn[netlist[i].b][nv+1]+=g;
      }
      else if (tipo=='V') {
        Yn[netlist[i].a][netlist[i].x]+=1;
        Yn[netlist[i].b][netlist[i].x]-=1;
        Yn[netlist[i].x][netlist[i].a]-=1;
        Yn[netlist[i].x][netlist[i].b]+=1;
        //if ((linear) || (vezes))
          Yn[netlist[i].x][nv+1]-=netlist[i].valor;
      }
      else if (tipo=='E') {
        g=netlist[i].valor;
        Yn[netlist[i].a][netlist[i].x]+=1;
        Yn[netlist[i].b][netlist[i].x]-=1;
        Yn[netlist[i].x][netlist[i].a]-=1;
        Yn[netlist[i].x][netlist[i].b]+=1;
        Yn[netlist[i].x][netlist[i].c]+=g;
        Yn[netlist[i].x][netlist[i].d]-=g;
      }
      else if (tipo=='F') {
        g=netlist[i].valor;
        Yn[netlist[i].a][netlist[i].x]+=g;
        Yn[netlist[i].b][netlist[i].x]-=g;
        Yn[netlist[i].c][netlist[i].x]+=1;
        Yn[netlist[i].d][netlist[i].x]-=1;
        Yn[netlist[i].x][netlist[i].c]-=1;
        Yn[netlist[i].x][netlist[i].d]+=1;
      }
      else if (tipo=='H') {
        g=netlist[i].valor;
        Yn[netlist[i].a][netlist[i].y]+=1;
        Yn[netlist[i].b][netlist[i].y]-=1;
        Yn[netlist[i].c][netlist[i].x]+=1;
        Yn[netlist[i].d][netlist[i].x]-=1;
        Yn[netlist[i].y][netlist[i].a]-=1;
        Yn[netlist[i].y][netlist[i].b]+=1;
        Yn[netlist[i].x][netlist[i].c]-=1;
        Yn[netlist[i].x][netlist[i].d]+=1;
        Yn[netlist[i].y][netlist[i].x]+=g;
      }
      else if (tipo=='C' || tipo=='L') {

        //se for capacitor, a condutancia em DC tende a 0
        if (tipo=='C')
          g = 1 / FATORDC;

        //se for indutor, a condutancia em DC tende a infinito
        if (tipo=='L')
          g = 0;

        if (tipo=='L')
        {
               Yn[netlist[i].a][netlist[i].x]+=1;
               Yn[netlist[i].b][netlist[i].x]-=1;
               Yn[netlist[i].x][netlist[i].a]-=1;
               Yn[netlist[i].x][netlist[i].b]+=1;
               Yn[netlist[i].x][netlist[i].x]+=1/FATORDC;
       }

        Yn[netlist[i].a][netlist[i].a]+=g;
        Yn[netlist[i].b][netlist[i].b]+=g;
        Yn[netlist[i].a][netlist[i].b]-=g;
        Yn[netlist[i].b][netlist[i].a]-=g;
      }
     else if (tipo=='M')
     {
          g = ((double) 1.0) / FATORDC;
          gm = 0;
          gds = 0;
          io = 0;

		  #ifdef DEBUG
          	  printf ("vd %f vs %f\n", vAtual[netlist[i].td],vAtual[netlist[i].ts] );
		  #endif


          if (vAtual[netlist[i].td] < vAtual[netlist[i].ts])
          {
			  #ifdef DEBUG
        	  	  printf ("vd > vs\n");
			  #endif
        	  int aux  = netlist[i].td;
        	  netlist[i].td = netlist[i].ts;
        	  netlist[i].ts = aux;
          }


          double vgs = vAtual[netlist[i].tg]-vAtual[netlist[i].ts];
          double vds = vAtual[netlist[i].td]-vAtual[netlist[i].ts];
          double vbs = vAtual[netlist[i].tb]-vAtual[netlist[i].ts];
          double vt = netlist[i].vt + netlist[i].gama * (sqrt((netlist[i].phi - vbs)) - sqrt((netlist[i].phi)));


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
               gm = netlist[i].k * (netlist[i].w/netlist[i].l)* (2.0*(vgs - vt)) * (1 + netlist[i].lambda* vds);
               gds = netlist[i].k * (netlist[i].w/netlist[i].l)* pow ((vgs - vt),2.0) * netlist[i].lambda;
               io = netlist[i].k * (netlist[i].w/netlist[i].l) * pow((vgs - vt),2.0) * (1 + netlist[i].lambda * vds) - (gm * vgs) - (gds * vds);
          }
          else if (netlist[i].transistorOp == ohmica)
          {
               gm = netlist[i].k * (netlist[i].w/netlist[i].l)*(2.0* vds)*(1+ netlist[i].lambda* vds);
               gds = netlist[i].k * (netlist[i].w/netlist[i].l) * (2.0*(vgs - vt) - 2 * vds + 4.0* netlist[i].lambda * (vgs - vt) * vds - 3.0*netlist[i].lambda * pow (vds,2.0));
               io = netlist[i].k * (netlist[i].w/netlist[i].l) * (2.0* (vgs - vt)*vds - pow (vds,2.0)) - (gm * vgs) - (gds * vds);
          }

          double gmb = (gm*netlist[i].gama)/(2*sqrt(fabs(netlist[i].phi - (vAtual[netlist[i].tb] -vAtual[netlist[i].ts]))));

          io-= (gmb * vbs);

		  #ifdef DEBUG
          	  printf ("gm %e gmb %e  gds %e io %e \n", gm, gmb, gds, io);
		  #endif

          Yn[netlist[i].td][netlist[i].tb]+=gmb;
          Yn[netlist[i].ts][netlist[i].ts]+=gmb;
          Yn[netlist[i].td][netlist[i].ts]-=gmb;
          Yn[netlist[i].ts][netlist[i].tb]-=gmb;

          Yn[netlist[i].td][netlist[i].tg]+=gm;
          Yn[netlist[i].ts][netlist[i].ts]+=gm;
          Yn[netlist[i].td][netlist[i].ts]-=gm;
          Yn[netlist[i].ts][netlist[i].tg]-=gm;

          Yn[netlist[i].td][netlist[i].td]+=gds;
          Yn[netlist[i].ts][netlist[i].ts]+=gds;
          Yn[netlist[i].td][netlist[i].ts]-=gds;
          Yn[netlist[i].ts][netlist[i].td]-=gds;


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
      else if (tipo=='O') {
        Yn[netlist[i].a][netlist[i].x]+=1;
        Yn[netlist[i].b][netlist[i].x]-=1;
        Yn[netlist[i].x][netlist[i].c]+=1;
        Yn[netlist[i].x][netlist[i].d]-=1;
      }

      //utilizar esse if para modificar as solucoes
      //atribui o valor calculado na ultima iteracao
      //a matriz que esta sendo calculada no momento
//      if (Linear && vezes)
//      {
//        for (int i=1; i<=nv; i++)
//        {
//          Yn[i][nv+1] = vAtual[i];
//          //usar uma randomizacao do C caso queira
//        }
//      }

  #ifdef DEBUG
      /* Opcional: Mostra o sistema apos a montagem da estampa */
      printf("Sistema apos a estampa de %s\n",netlist[i].nome);
      for (k=1; k<=nv; k++) {
        for (j=1; j<=nv+1; j++)
          if (Yn[k][j]!=0) printf("%+3.3f ",Yn[k][j]);
          else printf(" ... ");
        printf("\n");
      }
      getch();
  #endif
    }
    /* Resolve o sistema */
    if (resolversistema()) {
      getch();
      //exit;
      return 1;
    }
  #ifdef DEBUG
    /* Opcional: Mostra o sistema resolvido */
    printf("Sistema resolvido:\n");
    for (i=1; i<=nv; i++) {
        for (j=1; j<=nv+1; j++)
          if (Yn[i][j]!=0) printf("%+3.1f ",Yn[i][j]);
          else printf(" ... ");
        printf("\n");
      }
    getch();
  #endif
    /* Mostra solucao */
    printf("Solucao:\n");
    strcpy(txt,"Tensao");
    for (i=1; i<=nv; i++) {
      if (i==nn+1) strcpy(txt,"Corrente");
      printf("%s %s: %g\n",txt,lista[i],Yn[i][nv+1]);
      vProximo[i] = Yn[i] [nv+1];
      printf("VProximo: %f\n", vProximo[i] );
    }
    vezes++;
    //se for nao linear, utiliza Newton-Raphson
    if (!linear)
    {
      controleConvergencia (vAtual, vProximo, iteracoes);
    }
    else
      convergiu = true;
  }

  for (int i =1; i<ne; i++)
  {
      tipo=netlist[i].nome[0];
      if (tipo == 'M')
      {
    	  CalculaCapacitancias (&netlist[i]);
    	  printf ("Cgs= %e Cgd= %e Cgb= %e \n", netlist[i].cgs, netlist[i].cgd, netlist[i].cgb);
      }

  }
  getch();
  return 0;
}
