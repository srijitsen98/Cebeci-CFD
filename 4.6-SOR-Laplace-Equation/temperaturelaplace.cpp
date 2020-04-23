#include<stdio.h>
#include<math.h>
#define nx 10
#define ny 10
#define pi 3.14159
int main()
{
    double Aw,Ae,An,As,T[nx+1][ny+1],Tan[nx+1][ny+1],Told[nx+1][ny+1];
    double dy,dx;
    double err,w,Temp,Tempo;
    int it,i,j,loci,locj;
    dy=1.0/ny;
    dx=1.0/nx;
    w=2.0/(1+sin(pi*dx));
    //printf("\n%f %f",dx,dy);
    for(i=0;i<=nx;i++)
    {
        for(j=0;j<=ny;j++)
        {
            Tan[i][j]=exp(pi*i*dx)*sin(pi*j*dy);
        }
    }
    Aw=dy*dy/(2.0*(dx*dx+dy*dy));
    Ae=dy*dy/(2.0*(dx*dx+dy*dy));
    An=dx*dx/(2.0*(dx*dx+dy*dy));
    As=dx*dx/(2.0*(dx*dx+dy*dy));
    j=0;
    for(i=0;i<=nx;i++)
    {
        T[i][j]=0;
    }
    j=ny;
    for(i=0;i<=ny;i++)
    {
        T[i][j]=0;
    }
    i=0;
    for(j=1;j<ny;j++)
    {
        T[i][j]=sin(pi*j*dy);
    }
    i=nx;
    for(j=1;j<nx;j++)
    {
        T[i][j]=exp(pi)*sin(pi*j*dy);
    }
    for(i=1;i<nx;i++)
    {
        for(j=1;j<ny;j++)
        {
            T[i][j]=0;
        }
    }
    err=0;
    it=0;
    while(1)
    {
        err=0;
        it++;
        for(i=1;i<nx;i++)
        {
            for(j=1;j<ny;j++)
            {
                //Told[i][j]=T[i][j];
                T[i][j]=w*(Aw*T[i-1][j]+Ae*T[i+1][j]+As*T[i][j-1]+An*T[i][j+1])+(1-w)*T[i][j];
            }
        }
        err=0;
        for(i=1;i<nx;i++)
        {
            for(j=1;j<ny;j++)
            {
                if(err<fabs(T[i][j]-Ae*T[i+1][j]-Aw*T[i-1][j]-As*T[i][j-1]-An*T[i][j+1]))
                {
                    loci=i;
                    locj=j;
                    err=fabs(T[i][j]-Ae*T[i+1][j]-Aw*T[i-1][j]-As*T[i][j-1]-An*T[i][j+1]);
                }
            }
        }
        if(err<1e-6)
        {
            break;
        }
        printf("\nIteration no. %d Residual %.8f",it,err);
    }
    printf("\nx\t\t y\t\t FDS\t\t AS\t\t %% Error");
    i=2;
    for(j=1;j<ny;j+=1)
    {
        printf("\n%f\t %f\t %f\t %f\t %f",i*dx,j*dy,T[i][j],Tan[i][j],100.0*fabs(T[i][j]-Tan[i][j])/Tan[i][j]);
    }
    i=5;
    for(j=1;j<ny;j+=1)
    {
        printf("\n%f\t %f\t %f\t %f\t %f",i*dx,j*dy,T[i][j],Tan[i][j],100.0*fabs(T[i][j]-Tan[i][j])/Tan[i][j]);
    }
    i=9;
    for(j=1;j<ny;j+=1)
    {
        printf("\n%f\t %f\t %f\t %f\t %f",i*dx,j*dy,T[i][j],Tan[i][j],100.0*fabs(T[i][j]-Tan[i][j])/Tan[i][j]);
    }
    return 0;
}
