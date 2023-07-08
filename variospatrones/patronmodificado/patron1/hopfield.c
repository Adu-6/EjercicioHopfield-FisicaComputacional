#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gsl_rng.h"

//gsl-config --libs

gsl_rng *rnd;

int main()
{
    int i,j,k,l,m,n,o,x,y,N,Npatrones;
    int semilla,pregunta,pregunta2;
    double T,delta_E,p,aux;
    extern gsl_rng *rnd;
    FILE *f1,*f2,*f3;

    //Inicializo la semilla
    srand(time(NULL));
    semilla=288411288;
    rnd=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rnd,semilla);

    //Abro arhivos y defino N
    f1=fopen("evolucionpatron.txt","w");
    f2=fopen("patron.txt","r");
    f3=fopen("solapamiento.txt","w");
    N=30;
    Npatrones=5;
    pregunta=-1;
    pregunta2=-1;

    //pregunto por la temperatura por pantalla
    printf("Introduce la Temperatura (número real): ");
    scanf("%lf",&T);

    //defino las matrices
    int s[N][N],xi[Npatrones][N][N],theta[N][N];
    double w[N][N][N][N],a[Npatrones],solap[Npatrones];

    //inicializo la matriz theta
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            theta[i][j]=0;
        }
    }
    
    //construyo la matriz patrón introducida por un txt
    for(o=0;o<Npatrones;o++)
    {
        for (i=0;i<N; i++)
        {
            for (j=0;j<N;j++)
            {
                fscanf(f2,"%i", &xi[o][i][j]);
            }
        }
    }

    //pido el tipo de configuración por pantalla
    while (pregunta!=0 && pregunta!=1)
    {
        printf("¿Quieres matriz inicial aleatoria, o patrón deformado? (introduce 0 para generación aleatoria, 1 para que sea el patrón deformado): ");
        scanf("%i",&pregunta);
    }

    if(pregunta==0)
    {
        //construyo la matriz de configuración inicial con 1 y 0 aleatoriamente
        for (i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                x=gsl_rng_uniform_int(rnd,2); //Nº aleatorio entero entre 0 y 1
                s[i][j]=x;
            }
        }    
    }
    else
    {
        //pregunto de qué patrón se quiere hacer la modificación
        while (pregunta2<0 || pregunta2>10)
        {
            printf("¿De qué patron quieres hacer la deformación (1-5)?:");
            scanf("%i",&pregunta2);
        }

        //construyo la matriz inicial como el patrón deformado aleatoriamente
        //primero establezco la matriz sistema tal cual la patrón
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                s[i][j]=xi[pregunta2-1][i][j];
            }
        }
        //ahora modifico aleatoriamente un 10 por ciento de los puntos del sistema
        for(i=0;i<90;i++)
        {
            //genero dos números aleatorios entre 0 y N para escoger un punto aleatorio de la matriz
            x=gsl_rng_uniform_int(rnd,N);
            y=gsl_rng_uniform_int(rnd,N);

            if(s[x][y]==0)
            {
                s[x][y]=1;
            }
            else
            {
                s[x][y]=0;
            }
        }
    }

    //Pinto la matriz inicial en el txt
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(j<N-1)
            {
                fprintf(f1,"%i, ",s[i][j]);
            }
            else
            {
                fprintf(f1,"%i",s[i][j]);
            }
        }
        fprintf(f1,"\n");
    }
    fprintf(f1,"\n");

    //Calculo los valores de a para el patrón dado
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

    //ahora calculo los valores de w (pesos sinápticos) para el patrón dado
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

    //Calculo el solapamiento de la matriz inicial con el patrón y lo saco por archivo
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
        fprintf(f3,"%lf\t",solap[o]);
    } 
    fprintf(f3,"\n");

    //comienzo el proceso
    for(m=0;m<500;m++)
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
            //saco el solapamiento calculado por archivo
            fprintf(f3,"%lf\t",solap[o]);
        }
        fprintf(f3,"\n");

        //saco la matriz sistema por archivo
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                if(j<N-1)
                {
                    fprintf(f1,"%i, ",s[i][j]);
                }
                else
                {
                    fprintf(f1,"%i",s[i][j]);
                }
            }
            fprintf(f1,"\n");
        }
        fprintf(f1,"\n");
    }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    return 0;
}
