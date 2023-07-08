#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gsl_rng.h"

//gsl-config --libs

gsl_rng *rnd;

int main()
{
    const double tolerancia=1e-6;
    int i,j,k,l,m,n,o,x,y,N,Npatrones,contador;
    int semilla,valor;
    double T,delta_E,p,aux;
    extern gsl_rng *rnd;
    FILE *f1;

    //pregunto por el número de patrones
    printf("¿Cuántos patrones quieres generar? ");
    scanf("%i",&Npatrones);

    //Inicializo la semilla
    srand(time(NULL));
    semilla=288411288;
    rnd=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rnd,semilla);

    //Abro arhivos y defino N y T
    f1=fopen("serecuerda.txt","w");
    N=20;
    T=0.0001;

    //defino las matrices
    int s[N][N],xi[Npatrones][N][N],theta[N][N];
    double w[N][N][N][N],a[Npatrones],solap[Npatrones];

    //Repito el programa completo 200 veces
    for (contador=0;contador<500;contador++)
    {
        //construyo la matriz de patrones aleatoriamente
        for(o=0;o<Npatrones;o++)
        {
            for (i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        x=gsl_rng_uniform_int(rnd,2); //Nº aleatorio entero entre 0 y 1
                        xi[o][i][j]=x;
                    }
                }    
        }

        //inicializo la matriz theta
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                theta[i][j]=0;
            }
        }

        //construyo la matriz de configuración inicial con 1 y 0 aleatoriamente
        for (i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                x=gsl_rng_uniform_int(rnd,2); //Nº aleatorio entero entre 0 y 1
                s[i][j]=x;
            }
        }    

        //Calculo los valores de a para los patrones dados
        for(o=0;o<Npatrones;o++)
        {
            a[o]=0;
        }

        for(o=0;o<Npatrones;o++)
        {
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    a[o]=a[o]+xi[o][i][j];
                }
            }
            a[o]=a[o]/pow(N,2);
        }

        //ahora calculo los valores de w (pesos sinápticos) para los patrones dados
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                for(k=0;k<N;k++)
                {
                    for(l=0;l<N;l++)
                    {
                        w[i][j][k][l]=0;
                    }
                }
            }
        }

        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                for(k=0;k<N;k++)
                {
                    for (l=0;l<N;l++)
                    {
                        if (i==k && j==l)
                        {
                            for(o=0;o<Npatrones;o++)
                            {
                                w[i][j][k][l]=0;
                            }
                        }
                        else
                        {
                            for(o=0;o<Npatrones;o++)
                            {
                                w[i][j][k][l]=w[i][j][k][l]+(xi[o][i][j]-a[o])*(xi[o][k][l]-a[o]);
                            }
                            w[i][j][k][l]=(1.0/pow(N,2))*w[i][j][k][l];
                        }
                    }
                }   
            }
        }

        //calculo el umbral de disparo en función de los pesos sinápticos
        for(i=0;i<N;i++)
        {    
            for(j=0;j<N;j++)
            {
                for(k=0;k<N;k++)
                {
                    for(l=0;l<N;l++)
                    {
                        theta[i][j]=theta[i][j]+0.5*w[i][j][k][l];
                    }
                }
            }
        }

        //comienzo el proceso
        for(m=0;m<20;m++)
        {
            //Paso de MonteCarlo
            for(n=0;n<=pow(N,2);n++)
            {
                //Ahora genero dos números aleatorios entre 0 y N para escoger un punto aleatorio de la matriz
                x=gsl_rng_uniform_int(rnd,N);
                y=gsl_rng_uniform_int(rnd,N);

                //calculo ahora el incremento de energía para la casilla (x,y)
                delta_E=0;
                for (i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        delta_E=delta_E+w[x][y][i][j]*s[i][j];
                    }
                }
                delta_E=(-1.*0.5*delta_E+theta[x][y])*(1.-2.*s[x][y]);

                //probabilidad de transición 
                p=exp(-1.*delta_E/T);
                if(p>1)
                {
                p=1.;
                }

                //ahora genero un número entre 0 y 1
                aux=gsl_rng_uniform(rnd);

                //si me sale menor que la probabilidad, cambio el número de la matriz sistema
                if(aux<p)
                {
                    if(s[x][y]==1)
                    {
                        s[x][y]=0;
                    }
                    else
                    {
                        s[x][y]=1;
                    }
                }
            } //fin paso de montecarlo

            //Calculo el solapamiento al final de cada paso montecarlo
            for(o=0;o<Npatrones;o++)
            {
                solap[o]=0.;
            }
            for(o=0;o<Npatrones;o++)
            {
                for(i=0;i<N;i++)
                {
                    for(j=0;j<N;j++)
                    {
                        solap[o]=solap[o]+(xi[o][i][j]-a[o])*(s[i][j]-a[o]);
                    }
                }
                solap[o]=(1/(pow(N,2)*a[o]*(1-a[o])))*solap[o];
            }   
        }
        //compruebo si ha recordado algún patrón
        for(o=0;o<Npatrones;o++)
        {
            if(fabs(solap[o])>=0.75)
            {
                valor=1;
                fprintf(f1,"%i\n",valor);
                break;
            }
            else
            {
                valor=0;
            }
        }
        if (valor==0)
            {
                fprintf(f1,"%i\n",valor);
            }
    }
    fclose(f1);
    return 0;
}
