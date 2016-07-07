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
#include <iostream>
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_TIPO 5
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-15
#define DEBUG
#define FATORDC 1e-9 //fator multiplicativo para capacitores e indutores DC
#define MAX_IT 50 //maximo de iteracoes
#define REFVAL 1 //valor de referencia utilizado nos calculos de convergencias
#define MINCONV 2 //minimo de iteracoes para considerar como estavel a solucao
#define NCONV 5 // o numero de vezes que o algoritmo pode tentar calcular uma solucao
                //mesmo que esta nao esteja convergindo

#define J dcomplex(0.0,1.0)
#define UM 0.9999999999990
#define ZERO 0.0000000000001
typedef std::complex<double> dcomplex;
//#define PI  (4*atan(1))
#define PI acos(-1.0)
#define MAX_ERRO 1e-9
#define  NMAX  50

typedef struct {
     double R,I;  //algebraic form
  } COMPLEX;

typedef COMPLEX Matc[NMAX][NMAX];
  typedef COMPLEX Vecc[NMAX];
  typedef int Veci[NMAX];
 Matc A; Vecc B, X;



using namespace std;

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

struct cmd
{
  char modo[MAX_NOME];
  int ptos, ptInicio, ptFim;
  bool temComando;
} comando;

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
  gmb = 0,
  io = 0;

complex<double> Ycomp[MAX_NOS+1][MAX_NOS+2];

ptoOperacao transistorOpAtual [MAX_NOS +1]; //TIRANDO ISSO

bool ptOperacao = true; //comeca o programa na analise do ponto de operacao
bool mantemModelo = true; //manter modelo do transistor
bool convergiu = false; //indica que a solucao convergiu
int contadorConv = 0; //conta quantas vezes os calculos deram erros menores
int iteracoes = 0; //conta o numero de iteracoes do alogoritmo
int vezNConvergiu = 0; //Conta a quantidade de vezes que o algoritmo nao convergiu

int resolversistema(void);
int resolverSistemaAC(void);
int numero(char *nome);
void controleConvergencia ( double vAtual[], double vProximo[], int iteracoes );
void CalculaCapacitancias (elemento *netlist);
void montaEstampaDC();
void montaEstampaAC(double frequencia);
inline double sind (double angulo);
inline double cosd (double angulo);


double CABS(COMPLEX Z) {
    return sqrt(Z.R*Z.R + Z.I*Z.I);
  }

  void CADD(COMPLEX Z1, COMPLEX Z2, COMPLEX *Z) {
    Z->R = Z1.R + Z2.R;
    Z->I = Z1.I + Z2.I;
  }

  void CDIF(COMPLEX Z1, COMPLEX Z2, COMPLEX *Z) {
    Z->R = Z1.R - Z2.R;
    Z->I = Z1.I - Z2.I;
  }

  void CMUL(COMPLEX Z1, COMPLEX Z2, COMPLEX *Z) {
    Z->R = Z1.R*Z2.R - Z1.I*Z2.I;
    Z->I = Z1.R*Z2.I + Z1.I*Z2.R;
  }

  void CDIV(COMPLEX Z1, COMPLEX Z2, COMPLEX *Z) {
    double d; COMPLEX C;
    d = Z2.R*Z2.R+Z2.I*Z2.I;
    if (d>2.2e-16) {
      C.R=Z2.R; C.I=-Z2.I;
      CMUL(Z1,C,Z);
      Z->R=Z->R/d; Z->I=Z->I/d;
    }
  }

  void CSwap(COMPLEX *Z1, COMPLEX *Z2) {
    COMPLEX  C;
    C.R=Z1->R; C.I=Z1->I;
    Z1->R=Z2->R; Z1->I=Z2->I;
    Z2->R=C.R; Z2->I=C.I;
  }


  /****************************************************************
  * TSCGT procedure implements the triangularization algorithm of *
  * Gauss with full pivoting at each step for a complex matrix, A *
  * and saves the made transformations in KP and LP.              *
  * ------------------------------------------------------------- *
  * INPUTS:                                                       *
  *          N:   size of complex matrix A                        *
  *          A:   complex matrix of size N x N                    *
  * OUTPUTS;                                                      *
  *          it:  =0 if A is singular, else =1.                   *
  *           C:  contains the upper triangular matrix and the    *
  *               multipliers used during the process.            *
  *          KP:  contains the column exchanges.                  *
  *          LP:  contains the line exchanges.                    *
  ****************************************************************/
  void TSCGT(double eps, int N, Matc A, int *it, Matc C, Veci KP, Veci LP) {
  int i,j,k,k0,l0;
  COMPLEX C0,C1,P0,T0;

    for (i=1; i<=N; i++)
	  for(j=1; j<=N; j++) {
        C[i][j].R=A[i][j].R;
        C[i][j].I=A[i][j].I;
      }

    *it=1; k=1;
    while (*it==1 && k<N) {
      P0.R=C[k][k].R; P0.I=C[k][k].I;
      l0=k; k0=k;
      for (i=k; i<=N; i++)
        for (j=1; j<=N; j++)
		  if (CABS(C[i][j]) > CABS(P0)) {
            P0.R=C[i][j].R; P0.I=C[i][j].I;
            l0=i; k0=j;
          }
      LP[k]=l0; KP[k]=k0;
      if (CABS(P0) < eps)
        *it=0;
      else {
        if (l0!=k)
		  for (j=k; j<=N; j++) {
            T0.R=C[k][j].R; T0.I=C[k][j].I;
            C[k][j].R=C[l0][j].R;
            C[k][j].I=C[l0][j].I;
            C[l0][j].R=T0.R;
            C[l0][j].I=T0.I;
          }
        if (k0!=k)
		  for (i=1; i<=N; i++) {
            T0.R=C[i][k].R; T0.I=C[i][k].I;
            C[i][k].R=C[i][k0].R;
            C[i][k].I=C[i][k0].I;
            C[i][k0].R=T0.R;
            C[i][k0].I=T0.I;
          }
		for (i=k+1; i<=N; i++) {
          C0.R=C[i][k].R; C0.I=C[i][k].I;
		  CDIV(C0,P0,&C[i][k]);
          for (j=k+1; j<=N; j++) {
            C0.R=C[i][j].R; C0.I=C[i][j].I;
            CMUL(C[i][k],C[k][j],&C1);
            CDIF(C0,C1,&C[i][j]);
          }
        }
        k++;
      }
    }
    if (*it==1 && CABS(C[N][N]) < eps)  *it=0;
  } //TSCGT()

  /*************************************************************
  * Function BSCGT calculates the solution of upper triangular *
  * system [S(n-1)] by the back substitution method and deduces*
  * from it the solution of system [S]. The call must be made  *
  * only after a call to TSCGT and the matrix of [S] must be   *
  * regular.                                                   *
  * ---------------------------------------------------------- *
  * INPUTS:                                                    *
  *         N : size of matrix C                               *
  *         C:  contains the upper triangular matrix and the   *
  *             multipliers used during the triangularization  *
  *             process (in output of TSCGT).     .            *
  *         W : contains the right-side coefficients           *
  *             (modified in the process).                     *
  *         KP: contains the column exchanges.                 *
  *         LP: contains the line exchanges.                   *
  * OUTPUT:                                                    *
  *         Z : system solution complex vector.                *
  *************************************************************/
  void BSCGT(int N, Matc C, Vecc W, Veci KP, Veci LP, Vecc Z) {
    int i,j,k,k0,l0;
    COMPLEX C0,C1,S,Z0;

    Z0.R=0.0; Z0.I=0.0;
    for (k=1; k<N; k++) {
      l0=LP[k];
      if (l0!=k)  CSwap(&W[k],&W[l0]);
      for (i=k+1; i<=N; i++) {
        C0.R=W[i].R; C0.I=W[i].I;
        CMUL(C[i][k],W[k],&C1); CDIF(C0,C1,&W[i]);
      }
    }
    CDIV(W[N],C[N][N],&Z[N]);
    for (i=N-1; i>0; i--) {
      S.R=Z0.R; S.I=Z0.I;
      for (j=i+1; j<=N; j++) {
        C0.R=S.R; C0.I=S.I;
        CMUL(C[i][j],Z[j],&C1); CADD(C0,C1,&S);
      }
      CDIF(W[i],S,&C0);
      CDIV(C0,C[i][i],&Z[i]);
    }
    for (k=N-1; k>0; k--) {
      k0=KP[k];
      if (k0!=k)  CSwap(&Z[k],&Z[k0]);
    }
  } // BSCGT()


  /************************************************************
  * Solve a Complex Linear System AX = B By Gauss Method with *
  * full pivoting and a correction process.                   *
  * --------------------------------------------------------- *
  * INPUTS:                                                   *
  *         eps, dta : absolute and relative precisions       *
  *         M : maximum number of iterations                  *
  *         N : size of linear system                         *
  *         A : complex matrix                                *
  *         B : right-hand side (complex vector)              *
  * OUTPUTS:                                                  *
  *         it: flag, =-1 if no convergence, =0 if matrix A   *
  *             is singular, =1 convergence ok.               *
  *         X : contains system solution (X1,X2,...,Xn)       *
  *                                                           *
  ************************************************************/
  void RSLCGTC(double eps, double dta, int M, int N, Matc A, Vecc B,
	           int *it, Vecc X)  {
    int I,Jota,L;
    double phi1,phi2;
    COMPLEX C0,C1,S,Z0;
    Matc C; Vecc B1, W, Z;
	Veci KP, LP;

    Z0.R=0.0; Z0.I=0.0;

    TSCGT(eps,N,A,it,C,KP,LP);

	if (*it==1) {

// Save B in B1 before BSCGT
      for (Jota=1; Jota<=N; Jota++) B1[Jota] = B[Jota];

      BSCGT(N,C,B,KP,LP, X);

// Restore B after BSCGT
      for (Jota=1; Jota<=N; Jota++) B[Jota] = B1[Jota];

      *it=-1; L=1;
      while (*it==-1 && L<=M) {
        phi1=0.0;
        for (I=1; I<=N; I++)
          if (CABS(X[I]) > phi1)  phi1=CABS(X[I]);
		for (I=1; I<=N; I++) {
          S.R=Z0.R; S.I=Z0.I;
          for (Jota=1; Jota<=N; Jota++) {
            C0.R=S.R; C0.I=S.I;
            CMUL(A[I][Jota],X[Jota], &C1);
            CADD(C0,C1, &S);
          }
          CDIF(B[I],S, &W[I]);
        }

        BSCGT(N,C,W,KP,LP,Z);

        for (I=1; I<=N; I++) {
          C0.R=X[I].R; C0.I=X[I].I;
          CADD(C0,Z[I], &X[I]);
        }
        phi2=0.0;
        for (I=1; I<=N; I++)
          if (CABS(Z[I]) > phi2)  phi2=CABS(Z[I]);
        if (phi2/phi1 < dta) *it=1; else L++;
	  } // while
	} // if *it==1
  } // RSLCGTC()

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
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }

    else if (tipo=='I' || tipo=='V') {
      sscanf(p,"%10s%10s%lg%lg%lg",na,nb,&netlist[ne].modulo,&netlist[ne].fase,&netlist[ne].valor);
      printf("%s %s %s %g %g %g\n",netlist[ne].nome,na,nb,netlist[ne].modulo, netlist[ne].fase, netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }

    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[ne].valor);
      printf("%s %s %s %s %s %g\n",netlist[ne].nome,na,nb,nc,nd,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
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
      printf("%s %s %s %s %s %s\n",netlist[ne].nome,ntd,ntg,nts,ntb, tTipo);
      netlist[ne].pnmos = ((tTipo[0] =='N')?nmos:pmos);
      netlist[ne].td=numero (ntd);
      netlist[ne].tg=numero (ntg);
      netlist[ne].ts=numero (nts);
      netlist[ne].tb=numero (ntb);
      netlist[ne].transistorOp = corte;
      double u = (netlist[ne].pnmos == nmos?0.05:0.02);
      netlist[ne].cOx = 2* netlist[ne].k/u;
      linear = false;
    }
    else if (tipo=='*') { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      ne--;
    }
    else if (tipo=='.')
    {
      sscanf(p,"%10s%d%d%d",comando.modo, &comando.ptos, &comando.ptInicio, &comando.ptFim);
      comando.temComando = true;
      printf("\n\n\n");
      printf("Modo: %-s\n", comando.modo);
      printf("Pontos: %d\n", comando.ptos);
      printf("Inicio: %d\n", comando.ptInicio);
      printf("Fim: %d\n\n", comando.ptFim);
      ne--;
      getch();
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
  //getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d -> %s\n",i,lista[i]);
  //getch();
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
  //getch();


  int vezes = 0;

  while (!convergiu)
  {
    montaEstampaDC();
    #ifdef DEBUG
    /* Opcional: Mostra o sistema apos a montagem da estampa */
    if (i <= nv)
      printf("Sistema apos a estampa de %s\n",netlist[i].nome);

    for (k=1; k<=nv; k++) {
      for (j=1; j<=nv+1; j++)
        if (Yn[k][j]!=0) printf("%+3.3f ",Yn[k][j]);
        else printf(" ..... ");
      printf("\n");
    }
    //getch();
    #endif
    /* Resolve o sistema */
    if (resolversistema()) {
      //getch();
      //exit;
      return 1;
    }
  #ifdef DEBUG
  /* Opcional: Mostra o sistema resolvido */
  printf("Sistema resolvido:\n");
  for (i=1; i<=nv; i++) {
      for (j=1; j<=nv+1; j++)
        if (Yn[i][j]!=0) printf("%+3.1f ",Yn[i][j]);
        else printf(" ..... ");
      printf("\n");
    }
    #endif
    /* Mostra solucao */
    printf("Solucao:\n");
    strcpy(txt,"Tensao");
    for (i=1; i<=nv; i++) {
      if (i==nn+1) strcpy(txt,"Corrente");
      printf("%s %s: %g\n",txt,lista[i],Yn[i][nv+1]);
      vProximo[i] = Yn[i] [nv+1];
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
      printf ("Gm= %e Gds= %e Gmb= %e\n", gm, gds, gmb);
      printf ("Cgs= %e Cgd= %e Cgb= %e \n", netlist[i].cgs, netlist[i].cgd, netlist[i].cgb);
    }
  }
  if (convergiu)
  {
    char modo[MAX_NOME];
    strcpy(modo, comando.modo);
    if (!strcmp(modo, "OCT"))
    {
    }
    else if (!strcmp(modo, "OCT"))
    {
    }
    else
    {
      if (comando.ptInicio == 0)
      {
        comando.ptInicio = 10;
        comando.ptFim = 100000;
        comando.ptos = 1000;
      }
      int nPtos = 0;
      double passo = ( (comando.ptFim - comando.ptInicio)/ (comando.ptos - 1) );
      for (double freq = (double) comando.ptInicio; freq <= (double)comando.ptFim; freq += passo)
      {
        printf("\nFrequencia: %f\n", freq);
        montaEstampaAC(freq);
//     for (i=1; i<=nv; i++) {
//     for (j=1; j<=nv+1; j++)
//        if (abs(Ycomp[i][j])!=0) cout << Ycomp[i][j];//printf("%+3.1f ",Yn[i][j]);
//     else printf(" ......... ");
//          printf("\n");
//     }
//        double eps=1e-350; double dta=1e-350; int M=1; int it;
//
//        for (int col =1; col <nv+1; col ++)
//        {
//             for (int row = 1; row<nv+1; row++)
//             {
//               A[col][row].R = Ycomp[col][row].real();
//               A[col][row].I = Ycomp[col][row].imag();
//             }
//             B[col].R = Ycomp[col][nv+1].real();
//             B[col].I = Ycomp[col][nv+1].imag();
//        }
//
//        for (int col =1; col <nv+1; col ++)
//        {
//             for (int row = 1; row<nv+1; row++)
//               cout << A[col][row].R << " ";
//             cout << B[col].R <<endl;
//        }

//        RSLCGTC(eps,dta,M,nv,A,B,&it,X);
//
//  printf("it: %d", it);
//     printf("\n System solutions:\n");
//  for (int I=1; I<=nv; I++) {
//    printf("  X%d = %f", I, X[I].R);
//    if (X[I].I >= 0.0) printf(" + ");
//    else printf(" - ");
//    printf("%f I\n", fabs(X[I].I));
//  }
  //getch();

        //if (resolverSistemaAC()) exit(-1);
//         Mostra solucao
//        printf("Solucao:\n");
//        strcpy(txt,"Tensao");
//        for (int col =0; col <nv+1; col++)
//        {
//          Ycomp[col][nv+1] = X[col].R + X[col].I * J;
//        }


        resolverSistemaAC();
        //Mostra solucao
        printf("Solucao:\n");
        strcpy(txt,"Tensao");

        for (i=1; i<=nv; i++) {
          if (i==nn+1) strcpy(txt,"Corrente");
          printf("MODULO %s %s: %g\n",txt,lista[i], abs(Ycomp[i][nv+1]));
          printf("FASE %s %s: %g\n",txt,lista[i],( (180.0/ PI) *  arg(Ycomp[i][nv+1] ) ) );
        }
        nPtos++;
        printf("\nFrequencia: %f\n", freq);
        getch();
      }
      printf("Foram plotados %d pontos\n", nPtos);
    }
  }
  getch();
  return 0;
}

///* Resolucao de sistema de equacoes lineares.
//   Metodo de Gauss-Jordan com condensacao pivotal */
int resolversistema(void)
{
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

int resolverSistemaAC(void)
{
  int i,j,l, a;
  complex<double> t, p;
  p = t;

  for (i=1; i<=nv; i++) {
    t = 0.0 + J * 0.0;
    a=i;
    for (l=i; l<=nv; l++) {
      if (abs(Ycomp[l][i])>abs(t)) {
          a=l;
          t=Ycomp[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=nv+1; l++) {
          p=Ycomp[i][l];
          Ycomp[i][l]=Ycomp[a][l];
          Ycomp[a][l]=p;
      }
    }
    if (abs(t)<TOLG) {
      printf("Sistema singular\n");
      //printf ("%lf\n", (abs(t) - TOLG)*1000000000);
      return 1;
    }
    for (j=nv+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
      Ycomp[i][j]/= t;
      p=Ycomp[i][j];
      if (abs(p)!=0.0)  /* Evita operacoes com zero */
        for (l=1; l<=nv; l++) {
          if (l!=i)
               Ycomp[l][j]-=Ycomp[l][i]*p;
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
void controleConvergencia ( double vAtual[], double vProximo[], int iteracoes )
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
		netlist->cgs = netlist->cOx * netlist->w * netlist->ld;
		netlist->cgd = netlist->cOx * netlist->w * netlist->ld;
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

void montaEstampaDC()
{
 /* Monta o sistema nodal modificado */
  printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
  //getch();
  /* Zera sistema */
  for (int i=0; i<=nv+1; i++) {
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
      Yn[netlist[i].x][netlist[i].b]-=1;
      Yn[netlist[i].x][netlist[i].a]+=1;
      //if ((linear) || (vezes))
        Yn[netlist[i].x][nv+1]+=netlist[i].valor;
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
    else if (tipo=='C'){
      g = FATORDC;
      Yn[netlist[i].a][netlist[i].a]+=g;
      Yn[netlist[i].b][netlist[i].b]+=g;
      Yn[netlist[i].a][netlist[i].b]-=g;
      Yn[netlist[i].b][netlist[i].a]-=g;
    }

    //se for indutor, a condutancia em DC tende a infinito
    else if (tipo=='L')
    {
       g = FATORDC;
       Yn[netlist[i].a][netlist[i].x]+=1;
       Yn[netlist[i].b][netlist[i].x]-=1;
       Yn[netlist[i].x][netlist[i].a]-=1;
       Yn[netlist[i].x][netlist[i].b]+=1;
       Yn[netlist[i].x][netlist[i].x]+= g;
    }

   else if (tipo=='M')
     {
          g =  FATORDC;
          gm = 0;
          gds = 0;
          io = 0;

          double vds = vAtual[netlist[i].td]-vAtual[netlist[i].ts];

          if ((vds < 0 && netlist[i].pnmos == nmos) || (vds > 0 && netlist[i].pnmos == pmos))
          {
           	  int aux  = netlist[i].td;
           	  netlist[i].td = netlist[i].ts;
           	  netlist[i].ts = aux;
          }

		  #ifdef DEBUG
          	  printf ("vd %f vs %f vb %f vg %f\n", vAtual[netlist[i].td],vAtual[netlist[i].ts], vAtual[netlist[i].tb], vAtual[netlist[i].tg]);
		  #endif

          double vgs = vAtual[netlist[i].tg]-vAtual[netlist[i].ts];
          vds = vAtual[netlist[i].td]-vAtual[netlist[i].ts];
          double vbs = vAtual[netlist[i].tb]-vAtual[netlist[i].ts];

          vgs *= (netlist[i].pnmos==nmos?1:-1);
          vds *= (netlist[i].pnmos==nmos?1:-1);

          vbs = (fabs(vbs)>netlist[i].phi/2.0?netlist[i].phi/2.0:vbs);

          double vt = netlist[i].vt + netlist[i].gama * (sqrt((netlist[i].phi - fabs(vbs))) - sqrt((netlist[i].phi)));

          vbs *= (netlist[i].pnmos==nmos?1:-1);

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
          io*=(netlist[i].pnmos == pmos?-1:1);

	   	 #ifdef DEBUG
		 	   printf ("gm %e gmb %e  gds %e io %e \n\n", gm, gmb, gds, io);
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
    else
    {
         cout << "Deu ruim: " << endl;
         getch();
    }
  }
}

void montaEstampaAC(double frequencia)
{
 /* Monta o sistema nodal modificado */
  //printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
  //getch();
  /* Zera sistema */
  for (int i=0; i<=nv+1; i++) {
    //inicializa os vetores utilizdos na analise de convergencia
    for (int j=0; j<=nv+1; j++)
    {
      Ycomp[i][j]= 0.0 + 0.0*J;
    }
  }
  /* Monta estampas */
  for (int i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    //printf ("\ntipo :%c, i:%d\n", tipo, i);
    if (tipo=='R') {
      g=1/netlist[i].valor;
      Ycomp[netlist[i].a][netlist[i].a]+=g;
      Ycomp[netlist[i].b][netlist[i].b]+=g;
      Ycomp[netlist[i].a][netlist[i].b]-=g;
      Ycomp[netlist[i].b][netlist[i].a]-=g;
    }
    else if (tipo=='G') {
      g=netlist[i].valor;
      Ycomp[netlist[i].a][netlist[i].c]+=g;
      Ycomp[netlist[i].b][netlist[i].d]+=g;
      Ycomp[netlist[i].a][netlist[i].d]-=g;
      Ycomp[netlist[i].b][netlist[i].c]-=g;
    }
    else if (tipo=='I') {
//      double real = netlist[i].modulo * cosd(netlist[i].fase);
//      double imag = netlist[i].modulo * sind(netlist[i].fase);
//      complex <double> gComp = real + J*imag;
      complex<double> gComp = netlist[i].modulo * cosd(netlist[i].fase) + J * netlist[i].modulo * sind(netlist[i].fase);
      Ycomp[netlist[i].a][nv+1]-=gComp;
      Ycomp[netlist[i].b][nv+1]+=gComp;
    }
    else if (tipo=='V') {
//      double real = netlist[i].modulo * cosd(netlist[i].fase);
//      double imag = netlist[i].modulo * sind(netlist[i].fase);
//      complex <double> gComp = real + J*imag;
      complex<double> gComp = netlist[i].modulo * cosd(netlist[i].fase) + J * netlist[i].modulo * sind(netlist[i].fase);
//      cout << "gComp:" <<  gComp << endl;

      Ycomp[netlist[i].a][netlist[i].x]+=1;
      Ycomp[netlist[i].b][netlist[i].x]-=1;
      Ycomp[netlist[i].x][netlist[i].b]-=1;
      Ycomp[netlist[i].x][netlist[i].a]+=1;
      Ycomp[netlist[i].x][nv+1]+=gComp;
    }
    else if (tipo=='E') {
      g=netlist[i].valor;
      Ycomp[netlist[i].a][netlist[i].x]+=1;
      Ycomp[netlist[i].b][netlist[i].x]-=1;
      Ycomp[netlist[i].x][netlist[i].a]-=1;
      Ycomp[netlist[i].x][netlist[i].b]+=1;
      Ycomp[netlist[i].x][netlist[i].c]+=g;
      Ycomp[netlist[i].x][netlist[i].d]-=g;
    }
    else if (tipo=='F') {
      g=netlist[i].valor;
      Ycomp[netlist[i].a][netlist[i].x]+=g;
      Ycomp[netlist[i].b][netlist[i].x]-=g;
      Ycomp[netlist[i].c][netlist[i].x]+=1;
      Ycomp[netlist[i].d][netlist[i].x]-=1;
      Ycomp[netlist[i].x][netlist[i].c]-=1;
      Ycomp[netlist[i].x][netlist[i].d]+=1;
    }
    else if (tipo=='H') {
      g=netlist[i].valor;
      Ycomp[netlist[i].a][netlist[i].y]+=1;
      Ycomp[netlist[i].b][netlist[i].y]-=1;
      Ycomp[netlist[i].c][netlist[i].x]+=1;
      Ycomp[netlist[i].d][netlist[i].x]-=1;
      Ycomp[netlist[i].y][netlist[i].a]-=1;
      Ycomp[netlist[i].y][netlist[i].b]+=1;
      Ycomp[netlist[i].x][netlist[i].c]-=1;
      Ycomp[netlist[i].x][netlist[i].d]+=1;
      Ycomp[netlist[i].y][netlist[i].x]+=g;
    }
    else if (tipo=='C')
    {
      complex <double> gComp = 0.0 + J * (2.0 * (double)PI *netlist[i].valor * frequencia) ;
      Ycomp[netlist[i].a][netlist[i].a]+=gComp;
      Ycomp[netlist[i].b][netlist[i].b]+=gComp;
      Ycomp[netlist[i].a][netlist[i].b]-=gComp;
      Ycomp[netlist[i].b][netlist[i].a]-=gComp;

    }

    //se for indutor, a condutancia em DC tende a infinito
    else if (tipo=='L')
    {
         //complex <double> gComp = J * ((2.0 * PI *netlist[i].valor) * frequencia) ;
       complex <double> gComp = 0.0 + J * ((2.0 * (double)PI *netlist[i].valor) * frequencia) ;
       cout << "valor indutor: " << netlist[i].valor << endl;
       Ycomp[netlist[i].a][netlist[i].x]+=1.0 + 0.0*J;
       Ycomp[netlist[i].b][netlist[i].x]-=1.0 + 0.0*J;
       Ycomp[netlist[i].x][netlist[i].a]-=1.0+ 0.0*J;
       Ycomp[netlist[i].x][netlist[i].b]+=1.0 + 0.0*J;
       Ycomp[netlist[i].x][netlist[i].x]+=gComp;
    }

   else if (tipo=='M')
   {
        complex <double> gCgs = 2 * PI * netlist[i].cgs * J;
        complex <double> gCgd = 2 * PI * netlist[i].cgd * J;
        complex <double> gCgb = 2 * PI * netlist[i].cgb * J;

         Ycomp[netlist[i].td][netlist[i].td]+=gCgd;
         Ycomp[netlist[i].tg][netlist[i].tg]+=gCgd;
         Ycomp[netlist[i].td][netlist[i].tg]-=gCgd;
         Ycomp[netlist[i].tg][netlist[i].td]-=gCgd;

         Ycomp[netlist[i].ts][netlist[i].ts]+=gCgs;
         Ycomp[netlist[i].tg][netlist[i].tg]+=gCgs;
         Ycomp[netlist[i].ts][netlist[i].tg]-=gCgs;
         Ycomp[netlist[i].tg][netlist[i].ts]-=gCgs;

         Ycomp[netlist[i].tb][netlist[i].tb]+=gCgb;
         Ycomp[netlist[i].tg][netlist[i].tg]+=gCgb;
         Ycomp[netlist[i].tb][netlist[i].tg]-=gCgb;
         Ycomp[netlist[i].tg][netlist[i].tb]-=gCgb;

   }
    else if (tipo=='O') {
       Ycomp[netlist[i].a][netlist[i].x]+=1;
       Ycomp[netlist[i].b][netlist[i].x]-=1;
       Ycomp[netlist[i].x][netlist[i].c]+=1;
       Ycomp[netlist[i].x][netlist[i].d]-=1;
    }
  }
}

inline double sind (double angulo)
{
    double temp = sin( (angulo / 180.0) * PI );
    if (fabs(temp) > UM)
        return (1.0);
    if (fabs(temp) < ZERO)
        return (0.0);

    return (temp);
}

inline double cosd (double angulo)
{
    double temp = cos( (angulo / 180.0) * PI );
    if (fabs(temp) > UM)
        return (1.0);
    if (fabs(temp) < ZERO)
        return (0.0);

    return cos( (angulo / 180.0) * PI );
}
