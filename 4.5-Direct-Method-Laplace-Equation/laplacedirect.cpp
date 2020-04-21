#include<stdio.h>
#include<math.h>
#define dx 0.1 //spacing in x direction
#define dy 0.1 //spacing in y direction
#define L 1 //total length
#define dt 0.01 //time spacing
#define pi 3.14159
const int nx=(int)L/dx;
const int ny=(int)L/dy;
const int n=nx;
double b[n][n],c[n],x[n]; //variables that are used to return values after matrix operations
int gaussian(double A[n][n],double B[n])
{
    int i,j,k,a,b;
    double r,A1[n][n],B1[n];
    for(i=1;i<n;i++)
    {
        for(j=1;j<n;j++)
        {
            A1[i][j]=A[i][j];
        }
        B1[i]=B[i];
    }
    /*
    printf("\nA is\n");
    for(a=1;a<n;a++)
    {
        for(b=1;b<n;b++)
        {
            printf("%f ",A[a][b]);
        }
        printf("\n");
    }
    printf("\nB is \n");
    for(a=1;a<n;a++)
    {
        printf("%f\n",B[a]);
    }
    */
    for(k=1;k<n-1;k++)
    {
        for(i=k+1;i<n;i++)
        {
            r=A1[i][k]/A1[k][k];
            A1[i][k]=0;
            for(j=k+1;j<n;j++)
            {
                A1[i][j]=A1[i][j]-r*A1[k][j];
            }
            B1[i]=B1[i]-r*B1[k];
        }

    }
    for(i=n-1;i>=1;i--)
    {
        x[i]=B1[i];
        for(j=n-1;j>i;j--)
        {
            x[i]-=A1[i][j]*x[j];
        }
        x[i]=x[i]/A1[i][i];
        //printf("\nx[%d] is %f",i,x[i]);
    }

    return 0;
}
int multiplynxn(double x[n][n],double y[n][n]) //matrix multiplication function
{
    for(int i=1;i<n;i++)
    {
        for(int j=1;j<n;j++)
        {
            b[i][j]=0;
            for(int k=1;k<n;k++)
            {
                b[i][j]+=x[i][k]*y[k][j];
            }
        }
    }
    return 0;
}
int subtractnxn(double x[n][n],double y[n][n]) //matrix subtraction function
{
    for(int i=1;i<n;i++)
    {
        for(int j=1;j<n;j++)
        {
            b[i][j]=x[i][j]-y[i][j];
            //printf("\nx %f y %f b %f",x[i][j],y[i][j],b[i][j]);
        }
    }
    return 0;
}
int multiplynx1(double x[n][n],double z[n]) //matrix multiplication function
{
    for(int i=1;i<n;i++)
    {
        c[i]=0;
        for(int k=1;k<n;k++)
        {
            c[i]+=x[i][k]*z[k];
        }
    }
    return 0;
}
int subtractnx1(double x[n],double y[n]) //matrix subtraction function
{
    for(int i=1;i<n;i++)
    {
        c[i]=x[i]-y[i];
    }
    return 0;
}
int main()
{
    //double T[n];
    const int nx=(int)L/dx;
    const int ny=(int)L/dy;
    const int nt=(int)(0.1/dt);
    int i,mid,j,res,t,loc,i1,j1,k1;
    printf("\nValue of nx is %d and ny is %d",nx,ny);
    double Tan[nx+1][nx+1],T[nx+1][nx+1],A[nx][nx][nx],B[nx][nx][nx],C[nx][nx][nx],F[nx][nx],del[nx][nx][nx],phi[nx][nx],tau[nx][nx][nx],u[nx][nx],b1[nx];
    double s1, s2, s3, lamda, location,err,thetax,thetay;
    double alpha=1;
    double temp1[nx][nx],temp2[nx];

    /*Analytical Solution stated in book at x=0.3*/
    thetax=dy*dy/(2.0*(dx*dx+dy*dy));
    thetay=dx*dx/(2.0*(dx*dx+dy*dy));
    //printf("\n%f %f",thetax,thetay);
    lamda=dx/(alpha*dt);
    s1=1;
    s2=-1;
    s3=-lamda;

    for(i=0;i<=nx;i++)
    {
        for(j=0;j<=nx;j++)
        {
            T[i][j]=0;
            Tan[i][j]=exp(pi*i*dx)*sin(pi*j*dy);
        }
    }
    //Initialization of boundary values

    //West
    j=0;
    for(i=0;i<=nx;i++)
    {
        T[i][j]=0;
    }
    //East
    j=ny;
    for(i=0;i<=nx;i++)
    {
        T[i][j]=0;
    }
    //South
    i=0;
    for(j=1;j<ny;j++)
    {
        T[i][j]=sin(pi*j*dy);
    }
    //North
    i=nx;
    for(j=1;j<ny;j++)
    {
        T[i][j]=exp(pi)*sin(pi*j*dy);
    }



    //Assigning values to A,B and C
    //Left Boundary
    //A definition
    for(i1=1;i1<nx;i1++)
    {
        for(j1=1;j1<nx;j1++)
        {
            if(i1==j1)
            {
                A[1][i1][j1]=1;
            }
            else if(j1==(i1+1)||(j1==(i1-1)))
            {
                A[1][i1][j1]=-thetax;
            }
            else
            {
                A[1][i1][j1]=0;
            }
        }
    }
    //C definition
    for(i1=1;i1<nx;i1++)
    {
        for(j1=1;j1<nx;j1++)
        {
            if(i1==j1)
            {
                C[1][i1][j1]=-thetay;
            }
            else
            {
                C[1][i1][j1]=0;
            }
        }
    }
    //F definition
    F[1][1]=thetax*T[0][1]+thetay*T[1][0];
    for(i1=2;i1<nx-1;i1++)
    {
        F[1][i1]=thetay*T[i1][0];
    }
    F[1][nx-1]=thetax*T[nx][1]+thetay*T[nx-1][0];


    //internal cells
    //A definition
    for(k1=2;k1<nx-1;k1++)
    {
        for(i1=1;i1<nx;i1++)
        {
            for(j1=1;j1<nx;j1++)
            {
                if(i1==j1)
                {
                    A[k1][i1][j1]=1;
                }
                else if(j1==(i1+1)||(j1==(i1-1)))
                {
                    A[k1][i1][j1]=-thetax;
                }
                else
                {
                    A[k1][i1][j1]=0;
                }
            }
        }
    }
    //B & C definition
    for(k1=2;k1<nx-1;k1++)
    {
        for(i1=1;i1<nx;i1++)
        {
            for(j1=0;j1<nx;j1++)
            {
                if(i1==j1)
                {
                    B[k1][i1][j1]=-thetay;
                    C[k1][i1][j1]=-thetay;
                }
                else
                {
                    B[k1][i1][j1]=0;
                    C[k1][i1][j1]=0;
                }
            }
        }
    }
    //F definition
    for(k1=2;k1<nx-1;k1++)
    {
        F[k1][1]=thetax*T[0][k1];
        for(i1=2;i1<nx-1;i1++)
        {
            F[k1][i1]=0;
        }
        F[k1][nx-1]=thetax*T[nx][k1];
    }

    //Right boundary
    //A definition
    for(i1=1;i1<nx;i1++)
    {
        for(j1=1;j1<nx;j1++)
        {
            if(i1==j1)
            {
                A[nx-1][i1][j1]=1;
            }
            else if(j1==(i1+1)||(j1==(i1-1)))
            {
                A[nx-1][i1][j1]=-thetax;
            }
            else
            {
                A[nx-1][i1][j1]=0;
            }
        }
    }
    //B definition
    for(i1=1;i1<nx;i1++)
    {
        for(j1=1;j1<nx;j1++)
        {
            if(i1==j1)
            {
                B[nx-1][i1][j1]=-thetay;
            }
            else
            {
                B[nx-1][i1][j1]=0;
            }
        }
    }

    //F definition
    F[nx-1][1]=thetax*T[0][nx-1]+thetay*T[1][nx];
    for(i1=2;i1<nx-1;i1++)
    {
        F[nx-1][i1]=thetay*T[i1][nx];
    }
    F[nx-1][nx-1]=thetax*T[nx][nx-1]+thetay*T[nx-1][nx];


    //mid cell computation
    //del1=A1
    for(i1=1;i1<nx;i1++)
    {
        for(j1=1;j1<nx;j1++)
        {
                del[1][i1][j1]=A[1][i1][j1];
        }
        phi[1][i1]=F[1][i1];
    }
    printf("\n");
    for(i=2;i<=nx-1;i++)
    {
        for(i1=1;i1<nx;i1++)
        {
            for(k1=1;k1<nx;k1++)
                b1[k1]=B[i][i1][k1];
            for(k1=1;k1<nx;k1++)
            {
                for(j=1;j<nx;j++)
                {
                    temp1[k1][j]=del[i-1][j][k1];
                }
            }
            gaussian(temp1,b1);
            for(k1=1;k1<nx;k1++)
                tau[i][i1][k1]=x[k1];
        }
        multiplynxn(tau[i],C[i-1]);

        for(i1=1;i1<nx;i1++)
        {
            for(j1=1;j1<ny;j1++)
            {
                temp1[i1][j1]=b[i1][j1];
            }
        }

        subtractnxn(A[i],temp1);
        for(i1=1;i1<nx;i1++)
        {
            for(j1=1;j1<ny;j1++)
            {
                del[i][i1][j1]=b[i1][j1];
            }
        }

        multiplynx1(tau[i],phi[i-1]);
        for(i1=1;i1<nx;i1++)
        {
                temp2[i1]=c[i1];
        }
        subtractnx1(F[i],temp2);
        for(i1=1;i1<nx;i1++)
        {
                phi[i][i1]=c[i1];
        }
    }
    gaussian(del[nx-1],phi[nx-1]);
    for(k1=1;k1<nx;k1++)
    {
        u[nx-1][k1]=x[k1];
        T[k1][nx-1]=x[k1];
    }

    for(i=nx-2;i>=1;i--)
    {

        multiplynx1(C[i],u[i+1]);
        for(i1=1;i1<nx;i1++)
        {
            temp2[i1]=c[i1];
        }


        subtractnx1(phi[i],temp2);

        for(i1=1;i1<nx;i1++)
        {
                temp2[i1]=c[i1];
        }

        for(i1=1;i1<nx;i1++)
        {
            for(j1=1;j1<nx;j1++)
            {
                temp1[i1][j1]=del[i][i1][j1];
            }
            //printf("\n");
        }
        gaussian(temp1,temp2);

        for(k1=1;k1<nx;k1++)
        {
            u[i][k1]=x[k1];
            T[k1][i]=x[k1];
            //printf("\n%f",x[k1]);
        }


    }

    printf("\nFinal values of T\n");
    for(j=nx;j>=0;j--)
    {
        for(i=0;i<=nx;i++)
            printf("%f ",T[i][j]);
        printf("\n");
    }

    //final values of T
    printf("\nValues of T at x=0.2");
    i=2;
    for(j=1;j<nx;j++)
    {
        err=fabs(Tan[i][j]-T[i][j])*100.0/T[i][j];
        printf("\ny %f Numerical Temperature %f Analytical Temperature %f Error %f%%",j*dy,T[i][j],Tan[i][j],err);
    }
    printf("\nValues of T at x=0.5");
    i=5;
    for(j=1;j<nx;j++)
    {
        err=fabs(Tan[i][j]-T[i][j])*100.0/T[i][j];
        printf("\ny %f Numerical Temperature %f Analytical Temperature %f Error %f%%",j*dy,T[i][j],Tan[i][j],err);
    }
    printf("\nValues of T at x=0.9");
    i=9;
    for(j=1;j<nx;j++)
    {
        err=fabs(Tan[i][j]-T[i][j])*100.0/T[i][j];
        printf("\ny %f Numerical Temperature %f Analytical Temperature %f Error %f%%",j*dy,T[i][j],Tan[i][j],err);
    }

}
