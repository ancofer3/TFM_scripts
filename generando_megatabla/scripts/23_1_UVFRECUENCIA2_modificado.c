/***********************************************************************/
/*  U V F R E C U E N C I A 2 . C  :  Cadenas de longitud  = [2, 14]   */
/*                             Versi�n 1.0                             */
/*  Analiza cualquier fichero con DNA y da palabras con frecuencia>k   */
/***********************************************************************/
/*  1- Codificaci�n utilizada     {A:00, C:01, G:10, T:11}             */
/*  2- Elimina l�neas que empiezan con ">..."                          */
/*  3- Se come los <CR>. Los lee y los ignora.                         */
/*  4- Reinicia automata al encontrar "N" o "n" o comentario.          */
/************************************************************************/
/*       Vicente Arnau Llombart.    9 - VI - 2004                      */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#define LETRAS              4     // Es el numero de letras del alfabeto.
#define TRES                3     // Mascara
#define ENE                 5
#define MAYOR               62

#define CERO                 0
#define UNKBYTE           1024
#define UNMEGA         1048576

void imprimir_resultados(FILE *f_out, char *header, int *tabla1, long int MEMORIA, int WORD) {
   char letra[4] = {'A','C','G','T'};
   int i, j, ll, k, aux;
   // for i in range(MEMORIA): Para cada posicion de la tabla (que representa un kmero)

   for (i=0; i < MEMORIA; i++) {
      // if tabla1[i] == 0: continue (nos saltamos los kmeros con frecuencia = 0)
      if (tabla1[i] == 0) continue;
      aux = i;
      // else: 
      // inicializamos variable kmer
      char kmer[20];
      // Hacemos una movida para reconstruir el kmero
      // for j in range(long_word): (para cada pos del kmero)
      for (j=0; j < WORD; j++) {
         // variable auxiliar ll = i
         ll = aux;
         //for k in range(j,word-1): 
               // Nos quitamos un nt de la derecha en ll (los bits que lo representan)
         for (k=j; k < WORD-1;k++) ll >>= 2;
         // sacamos el nt al final de ll
         ll = ll & 3;
         // kmer en posicion j es letra[ll]
         kmer[j] = letra[ll];
      }
      kmer[WORD] = '\0';
      //printeamos el header, el kmero y la frecuencia en i   
      fprintf(f_out,"%s\t%s\t%d\n",header,kmer,tabla1[i]);
   }
}
int main( int argc, char *argv[])
{
 FILE    *finA, *f_out;
 char    ficha1[80], ficha2[80];
 char    caracter, linea[180];

 unsigned char   ta_lut[256];
 int register   puntero;
 int register   leido;
 int     i, j, k, aux;
 int     *tabla1=NULL;  // tabla con las frecuencias de cromosoma 1.

 int     encontrado=0, contador;
 double  faux;
 int     palabras=0;
 int     WORD;         // Longitud Palabra
 long int     MEMORIA, aux_long;      // 4  exp MAX_CADENA
 long int     long1, conta;
 int     MASCARA=0;    //  Mascara  los   MAX_CAD nucleotidos

 char    letra[4]={'A', 'C', 'G', 'T'};
 int     ll, WORD_1, bases, FF;

 double  aux_double, pal=0.0;
 char    *ssaux;

 char header[200] = "seq0"; // Variable header
 

 printf (  "=============================================");
 printf ("\n=     U  V  F  R  E  C  U  E N C I A 2      =");
 printf ("\n=============================================");
 printf ("\n=     Vicente Arnau  :   11 - II - 2026     =");
 printf ("\n=============================================\n");

 if (argc<3) {
     printf("\nSINTAXIS: UVFRECU2  ChrA  k\n");
     printf("\n \t  k   : Size Word [2, 13]\n");
     exit(0);
     }

//======================    FICHEROS con ADN     ==========
 strcpy(ficha1, argv[1]);       //********* FICHERO DE ENTRADA 
 
 if ( (finA=fopen(ficha1,"r"))!=NULL)
      printf ("\n Chromosome File  :  %s\n", ficha1);
 else {
      printf ("\n What about %s ??.\n", argv[1]);
      exit(0);
      }
      	  

  WORD = atoi(argv [2]);     //**************  Tama�o de palabra  
  if ( (WORD<2) || (WORD>14) ) {
     printf ("\n k = [2, 14]\n");
     exit(0);
     }
  WORD_1 = WORD - 1;       // Acceso a estado_FINAL         


  MEMORIA = 1;             //  MEMORIA = 4 exp WORD !!
  MASCARA = 0;
  // Generamos un número binario con k*2 Stellen de todo unos
  for (i=0; i<WORD; i++){
       MEMORIA = MEMORIA * 4;
       MASCARA = MASCARA << 2; // Añadimos dos 0s a la derecha
       MASCARA = MASCARA | 3; // Convertimos en 1 esos dos ceros
       }

  aux_long = (long int) sizeof(int) * MEMORIA; // Intenta reservar memoria para el kmero más largo (?)

  if ((aux_long) == 0) {
      printf("\n Memory allocation failed.\n");
      exit(EXIT_FAILURE);
      }
  if (aux_long < UNKBYTE)    //--->> Nuevo informe sobre memoria utilizada
       printf ("\n You need %ld bytes of memory.\n", aux_long);
  else {
        if (aux_long < UNMEGA){
            aux_double = (double) aux_long;
            aux_double = aux_double/(double) UNKBYTE;
            printf ("\n You need %5.0f Kbytes of memory.\n", aux_double);
            }
        else {
               aux_double = (double) aux_long;
               aux_double = aux_double/(double) UNMEGA;
               printf ("\n You need %5.0f Mbytes of memory.\n", aux_double);
               }
        }

  printf ("\n==================================================\n");

/********** BUSQUEDA DE CADENAS DE LONGITUD N_WORD bases    ***********/

// Inicializo la tabla Y RESERVO MEMORIA:

  tabla1 = (int *) calloc (MEMORIA, sizeof(int));
  if (tabla1 == NULL)  {
      printf("\n Memory allocation failed (1).\n");
      exit(EXIT_FAILURE);
      }

  for (i=0; i<256; i++) ta_lut[i] = 4;   // inicializo tabla lut
  ta_lut['A'] = 0;  ta_lut['a'] = 0;
  ta_lut['C'] = 1;  ta_lut['c'] = 1;
  ta_lut['G'] = 2;  ta_lut['g'] = 2;
  ta_lut['T'] = 3;  ta_lut['t'] = 3;
  ta_lut['N'] = ENE;  ta_lut['n'] = ENE;
  ta_lut['>'] = MAYOR;

 //=======================================================================
 //===========  FICHERO  DE  SALIDA     ==================================

  strcpy(ficha2, argv[1]); // Le ponemos a ficha2 el nombre del input
  strcat(ficha2, "_K_"); //Le añadimos _K_
  strcat(ficha2, argv[2]); // Le añadimos el valor de K
  strcat(ficha2, ".txt");  // Fichero  de salida
  f_out=fopen(ficha2,"w");
  fprintf(f_out, "SeqID\tkmer\tcount\n");

//============================================================
//============    S E C U E N C I A      A    ================

  long1=0;             // inincializo numero de K-mers  ------
  puntero = 0;
  bases = 0;
  caracter= getc(finA);
  

  while (caracter != EOF) {
      leido = ta_lut[caracter];
      if (leido < LETRAS){
         if (bases==WORD_1){
             puntero = puntero << 2;
             puntero = puntero & MASCARA;
             puntero = puntero | leido;
             tabla1 [puntero]++;
             long1++;
             }
          else {
             bases++;
             puntero = puntero << 2;
             puntero = puntero & MASCARA;
             puntero = puntero | leido;
             }
          }
       else {
             // Hay que incluir un algo para que la primera vez se guarde el header y no imprima nada, n0(¿?)
             if (leido == MAYOR) {
               if (long1 > 0) {
                  imprimir_resultados(f_out,header,tabla1,MEMORIA,WORD);
               }
                // CAMBIO para guardarnos nombre sec
                ssaux=fgets(linea, 160, finA);
                if (ssaux != NULL) {
                  // eliminar salto de línea si existe
                  linea[strcspn(linea, "\r\n")] = '\0';
                  strcpy(header, linea);
                  
               }
                bases = 0; 
                puntero = 0;    // INICIALIZO EL AUTOMATA !!!!!
                memset(tabla1, 0, MEMORIA * sizeof(int));
                long1 = 0;
                }
             if (leido == ENE){
               bases = 0; puntero = 0;    // INICIALIZO EL AUTOMATA !!!!!
               }
            } // del else

       caracter= getc(finA);              // Leo nuevo nucleotido

      } // del While
      imprimir_resultados(f_out, header, tabla1, MEMORIA, WORD);
 fclose(finA);            // cierro fichero de entrada CROMOSOMA  A !!!

  printf("\n Output File:  %s \n",ficha2);


//===============   C I E R R O   Y   M E   V O Y     ========

 fclose(f_out);   free(tabla1);
 return (LETRAS);
 }







