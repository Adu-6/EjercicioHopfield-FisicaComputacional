#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gsl_rng.h"

//gsl-config --libs

gsl_rng *rnd;

int main()
{
    int i,j,k,l,m,n,x,y,N;
    int semilla,pregunta;
    double T,delta_E,p,a,aux,solap;
    extern gsl_rng *rnd;
    FILE *f1,*f2,*f3;

    //Inicializo la semilla
    srand(time(NULL));
    semilla=rand();
    rnd=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rnd,semilla);

    //Abro arhivos y defino N
    f1=fopen("evolucionpatron.txt","w");
    f2=fopen("patron.txt","r");
    f3=fopen("solapamiento.txt","w");
    N=30;
    pregunta=-1;

    //pregunto por la temperatura por pantalla
    printf("Introduce la Temperatura (número real): ");
    scanf("%lf",&T);

    //defino las matrices
    int s[N][N],xi[N][N],theta[N][N];
    double w[N][N][N][N];

    //inicializo la matriz theta
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            theta[i][j]=0;
        }
    }
    
    //construyo la matriz patrón introducida por un txt
    for (i=0;i<N; i++)
    {
        for (j=0;j<N;j++)
        {
            fscanf(f2,"%i", &xi[i][j]);
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
        //construyo la matriz inicial como el patrón deformado aleatoriamente
        //primero establezco la matriz sistema tal cual la patrón
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                s[i][j]=xi[i][j];
            }
        }
        //ahora modifico aleatoriamente la mitad de puntos de la matriz sistema
        for(i=0;i<450;i++)
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
    a=0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            a=a+xi[i][j];
        }
    }
    a=a/pow(N,2);

    //ahora calculo los valores de w (pesos sinápticos) para el patrón dado
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
                        w[i][j][k][l]=0;
                    }
                    else
                    {
                        w[i][j][k][l]=(1/pow(N,2))*(xi[i][j]-a)*(xi[k][l]-a);
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
    solap=0.;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            solap=solap+(xi[i][j]-a)*(s[i][j]-a);
        }
    }
    solap=(1/(pow(N,2)*a*(1-a)))*solap;
        
    fprintf(f3,"%lf\n",solap);

    //comienzo el proceso
    for(m=0;m<50;m++)
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
        solap=0.;
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                solap=solap+(xi[i][j]-a)*(s[i][j]-a);
            }
        }
        solap=(1/(pow(N,2)*a*(1-a)))*solap;

        //saco el solapamiento calculado por archivo
        fprintf(f3,"%lf\n",solap);

        //saco la matriz por archivo
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
