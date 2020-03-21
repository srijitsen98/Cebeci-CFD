#include<stdio.h>
#include<math.h>
#define m 100
#define n 100
#define pi 3.14159
int main()
{
    double Aw,Ae,An,As,T[m+1][n+1],Tan[m+1][n+1],Told[m+1][n+1];
    double dy,dx;
    double err,w,Temp,Tempo;
    int it,i,j,loci,locj;
    w=1.0;
    dy=1.0/n;
    dx=1.0/m;
    printf("\n%f %f",dx,dy);
    for(i=0;i<=m;i++)
    {
        for(j=0;j<=n;j++)
        {
            Tan[i][j]=exp(pi*i*dx)*sin(pi*j*dy);
        }
    }
    Aw=dy*dy/(2.0*(dx*dx+dy*dy));
    Ae=dy*dy/(2.0*(dx*dx+dy*dy));
    An=dx*dx/(2.0*(dx*dx+dy*dy));
    As=dx*dx/(2.0*(dx*dx+dy*dy));
    j=0;
    for(i=0;i<=m;i++)
    {
        T[i][j]=0;
    }
    j=n;
    for(i=0;i<=m;i++)
    {
        T[i][j]=0;
    }
    i=0;
    for(j=1;j<n;j++)
    {
        T[i][j]=sin(pi*j*dy);
    }
    i=m;
    for(j=1;j<n;j++)
    {
        T[i][j]=exp(pi)*sin(pi*j*dy);
    }
    for(i=1;i<m;i++)
    {
        for(j=1;j<n;j++)
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
        for(i=1;i<m;i++)
        {
            for(j=1;j<n;j++)
            {
                Told[i][j]=T[i][j];
                T[i][j]=w*(Aw*T[i-1][j]+Ae*T[i+1][j]+As*T[i][j-1]+An*T[i][j+1])+(1-w)*T[i][j];
                if(fabs(T[i][j]-Told[i][j])>err)
                {
                    err=fabs(T[i][j]-Told[i][j]);
                    loci=i;
                    locj=j;
                    Temp=T[i][j];
                    Tempo=Told[i][j];
                }
            }
        }
        if(err<1e-6)
        {
            break;
        }
        printf("\nIteration no. %d Error %f i %d j %d T %f Told %f",it,err,loci,locj,Temp,Tempo);
    }
    printf("\nx\t\t y\t\t FDS\t\t AS\t\t %% Error");
    i=20;
    for(j=10;j<n;j+=10)
    {
        printf("\n%f\t %f\t %f\t %f\t %f",i*dx,j*dy,T[i][j],Tan[i][j],100.0*fabs(T[i][j]-Tan[i][j])/Tan[i][j]);
    }
    i=50;
    for(j=10;j<n;j+=10)
    {
        printf("\n%f\t %f\t %f\t %f\t %f",i*dx,j*dy,T[i][j],Tan[i][j],100.0*fabs(T[i][j]-Tan[i][j])/Tan[i][j]);
    }
    i=90;
    for(j=10;j<n;j+=10)
    {
        printf("\n%f\t %f\t %f\t %f\t %f",i*dx,j*dy,T[i][j],Tan[i][j],100.0*fabs(T[i][j]-Tan[i][j])/Tan[i][j]);
    }
    return 0;
}
