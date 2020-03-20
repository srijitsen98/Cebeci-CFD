#include<stdio.h>
#include<math.h>
#define dt 0.01
#define dx 0.1
#define L 1.0
//#define r 1.0//alpha*dt/(dx*dx)<0.5
#define alpha 1.0
#define pi 3.14159
double b[2][2],c[2];
int inv(double a[2][2])
{
    double det;
    det=a[0][0]*a[1][1]-a[0][1]*a[1][0];
    b[0][0]=a[1][1]/(det*1.0);
    b[0][1]=-a[0][1]/(det*1.0);
    b[1][0]=-a[1][0]/(det*1.0);
    b[1][1]=a[0][0]/(det*1.0);
    return 0;
}
int multiply2x2(double x[2][2],double y[2][2])
{
    b[0][0]=x[0][0]*y[0][0]+x[0][1]*y[1][0];
    b[0][1]=x[0][0]*y[0][1]+x[0][1]*y[1][1];
    b[1][0]=x[1][0]*y[0][0]+x[1][1]*y[1][0];
    b[1][1]=x[1][0]*y[0][1]+x[1][1]*y[1][1];
    return 0;
}
int subtract2x2(double x[2][2],double y[2][2])
{
    b[0][0]=x[0][0]-y[0][0];
    b[0][1]=x[0][1]-y[0][1];
    b[1][0]=x[1][0]-y[1][0];
    b[1][1]=x[1][1]-y[1][1];
    return 0;
}
int multiply2x1(double x[2][2],double z[2])
{
    c[0]=x[0][0]*z[0]+x[0][1]*z[1];
    c[1]=x[1][0]*z[0]+x[1][1]*z[1];
    return 0;
}
int subtract2x1(double x[2],double y[2])
{
    c[0]=x[0]-y[0];
    c[1]=x[1]-y[1];
}
int main()
{
    const int n=ceil(L/dx);
    double h[n+1];
    double Tnum[6],Tan[6],time[6];
    int i,mid,t,res,loc;
    double location,d[n];
    /*
    class system_variables
    {
        public:
            double T[n+1],x[n+1],Tnew[n+1];
    }f;*/
    /*Analytical Solution stated in book at x=0.3*/
    location=0.3;
    Tan[0]=.6004;
    Tan[1]=.5799;
    Tan[2]=.5334;
    Tan[3]=.4857;
    Tan[4]=.4411;
    Tan[5]=.4000;
    for(i=0;i<=n;i++)
    {
        h[i]=dx;
    }
    double dels[n+1][2],del[n+1][2][2],delinv[n+1][2][2],A[n+1][2][2],B[n+1][2][2],C[n+1][2][2],Tau[n+1][2][2];
    double lamda[n+1],Thalfold[n+1],Rhalfold[n+1],Told[n+1],Tnew[n+1],Tw,Te,s1[n+1],s2[n+1],s3[n+1],r[n+1][2],w[n+1][2];
    double pold[n+1],pnew[n+1],temp[2][2],temp1[2],x[n+1];

    /*Initial Conditions*/
    mid=ceil((n-1)/2.0);
    printf("\nMiddle location is %d",mid);
    for(i=1;i<mid;i++)
    {
        x[i]=dx*i;
        Told[i]=2.0*x[i];
    }
    for(i=mid;i<n;i++)
    {
        x[i]=dx*i;
        Told[i]=2.0*(1-x[i]);
    }
    loc=(int)ceil(location/dx);
    printf("\n%f",loc*dx);
    printf("\nloc is %d",loc);
    res=0;
    Tnum[res]=Told[loc];
    time[res]=0;
    res++;
    printf("\nT at time %f and location %f is %f",t*dt,loc*dx,Tnew[loc]);
    Tw=0;
    Te=0;
    for(t=1;t<=10;t++)
    {
        //i=0;
        r[0][0]=Tw;
        r[0][1]=0;
        for(i=1;i<n;i++)
        {
            lamda[i]=(1/alpha)*(h[i]/dt);
            s1[i]=1;
            s2[i]=-1;
            s3[i]=-lamda[i];
            Thalfold[i]=(Told[i]+Told[i-1])/2.0;
            Rhalfold[i]=-2.0*lamda[i]*Thalfold[i]+pold[i-1]-pold[i];
            r[i][0]=Rhalfold[i];
            r[i][1]=0;
        }
        //i=n
        lamda[n]=(1/alpha)*(h[n]/dt);
        s1[n]=1;
        s2[n]=-1;
        s3[n]=-lamda[n];
        Thalfold[n]=(Told[n]+Told[n-1])/2.0;
        Rhalfold[n]=-2.0*lamda[n]*Thalfold[n]+pold[n-1]-pold[n];
        r[n][0]=Rhalfold[n];
        r[n][1]=Te;
        A[0][0][0]=1;
        A[0][0][1]=0;
        A[0][1][0]=-1;
        A[0][1][1]=-h[1]/2.0;
        C[0][0][0]=0;
        C[0][0][1]=0;
        C[0][1][0]=1;
        C[0][1][1]=-h[1]/2.0;
        for(i=1;i<n;i++)
        {
            B[i][0][0]=s3[i];
            B[i][0][1]=s2[i];
            B[i][1][0]=0;
            B[i][1][1]=0;
            A[i][0][0]=s3[i];
            A[i][0][1]=s1[i];
            A[i][1][0]=-1;
            A[i][1][1]=-h[i+1]/2.0;
            C[i][0][0]=0;
            C[i][0][1]=0;
            C[i][1][0]=1;
            C[i][1][1]=-h[i+1]/2.0;
        }
        B[n][0][0]=s3[n];
        B[n][0][1]=s2[n];
        B[n][1][0]=0;
        B[n][1][1]=0;
        A[n][0][0]=s3[n];
        A[n][0][1]=s1[n];
        A[n][1][0]=1;
        A[n][1][1]=0;
        del[0][0][0]=A[0][0][0];
        del[0][0][1]=A[0][0][1];
        del[0][1][0]=A[0][1][0];
        del[0][1][1]=A[0][1][1];
        w[0][0]=r[0][0];
        //w[0][0][1]=r[0][0][1];
        w[0][1]=r[0][1];
        //w[0][1][1]=r[0][1][1];
        for(i=1;i<=n;i++)
        {
            inv(del[i-1]);
            delinv[i-1][0][0]=b[0][0];
            delinv[i-1][0][1]=b[0][1];
            delinv[i-1][1][0]=b[1][0];
            delinv[i-1][1][1]=b[1][1];
            multiply2x2(B[i],delinv[i-1]);
            Tau[i][0][0]=b[0][0];
            Tau[i][0][1]=b[0][1];
            Tau[i][1][0]=b[1][0];
            Tau[i][1][1]=b[1][1];
            multiply2x2(Tau[i],C[i-1]);
            temp[0][0]=b[0][0];
            temp[0][1]=b[0][1];
            temp[1][0]=b[1][0];
            temp[1][1]=b[1][1];
            subtract2x2(A[i],temp);
            del[i][0][0]=b[0][0];
            del[i][0][1]=b[0][1];
            del[i][1][0]=b[1][0];
            del[i][1][1]=b[1][1];
            multiply2x1(Tau[i],w[i-1]);
            temp1[0]=c[0];
            temp1[1]=c[1];
            subtract2x1(r[i],temp1);
            w[i][0]=c[0];
            w[i][1]=c[1];
        }
        inv(del[n]);
        delinv[n][0][0]=b[0][0];
        delinv[n][0][1]=b[0][1];
        delinv[n][1][0]=b[1][0];
        delinv[n][1][1]=b[1][1];
        multiply2x1(delinv[n],w[n]);
        dels[n][0]=c[0];
        dels[n][1]=c[1];
        for(i=n-1;i>=0;i--)
        {
            multiply2x1(C[i],dels[i+1]);
            temp1[0]=c[0];
            temp1[1]=c[1];
            subtract2x1(w[i],temp1);
            temp1[0]=c[0];
            temp1[1]=c[1];
            multiply2x1(delinv[i],temp1);
            dels[i][0]=c[0];
            dels[i][1]=c[1];
        }
        for(i=0;i<=n;i++)
        {
            Tnew[i]=dels[i][0];
            pnew[i]=dels[i][1];
        }
        for(i=0;i<=n;i++)
        {
            Told[i]=Tnew[i];
            pold[i]=pnew[i];
        }
        if((t==(int)ceil(0.01/dt))||(t==(int)ceil(0.02/dt))||(t==(int)ceil(0.03/dt))||(t==(int)ceil(0.04/dt))||(t==(int)ceil(0.05/dt)))
        {
            printf("\nT at time %f and location %f is %f",t*dt,loc*dx,Tnew[loc]);
            time[res]=t*dt;
            Tnum[res]=Tnew[loc];
            res++;
        }
    }

    /*Comparison between analytical and num*/
    printf("\nt \t\t Box \t\t AS \t\t Diff \t\t Rel Error (%%)");
    for(i=0;i<6;i++)
    {
        printf("\n%f \t %f \t %f \t %f \t %f",time[i],Tnum[i],Tan[i],fabs(Tnum[i]-Tan[i]),100.0*fabs(Tnum[i]-Tan[i])/Tan[i]);
    }
}

