#include<stdio.h>
#include<math.h>
#define dt 0.01
#define dx 0.1
#define L 1.0
#define r 1.0//alpha*dt/(dx*dx)<0.5
#define alpha 1.0
#define pi 3.14159
int main()
{
    const int n=ceil(L/dx)+1;
    double a[n],b[n],c[n],w[n];
    double Tnum[6],Tan[6],time[6];
    int i,mid,t,res,loc;
    double location,lamda,d[n];
    class system_variables
    {
        public:
            double T[n],x[n],Tnew[n];
    }f;
    /*Analytical Solution stated in book at x=0.3*/
    location=0.5;
    Tan[0]=.9904;
    Tan[1]=.7743;
    Tan[2]=.6808;
    Tan[3]=.6091;
    Tan[4]=.5488;
    Tan[5]=.4959;
    /*Initial Conditions*/
    mid=ceil((n-1)/2.0);
    printf("\nMiddle location is %d",mid);
    for(i=1;i<mid;i++)
    {
        f.x[i]=dx*i;
        f.T[i]=2.0*f.x[i];
        //printf("\n%f ",f.T[i]);
    }
    for(i=mid;i<n-1;i++)
    {
        f.x[i]=dx*i;
        f.T[i]=2.0*(1-f.x[i]);
        //printf("%f ",f.T[i]);
    }
    lamda=(2.0*dx*dx)/(alpha*dt);
    loc=(int)ceil(location/dx);
    printf("\nLocation of probe is %f",loc*dx);
    printf("\nElement no. loc is %d",loc);
    res=0;
    Tnum[res]=f.T[loc];
    time[res]=0;
    res++;
    t=0;
    printf("\nT at time %f and location %f is %f",t*dt,loc*dx,f.T[loc]);
    /*Boundary*/
    f.T[0]=0;
    f.T[n-1]=0;


    for(t=1;t<=10;t++)
    {
        /*Coefficients*/
        for(i=1;i<n-1;i++)
        {
            a[i]=1;
            b[i]=-(lamda+2);
            c[i]=1;
            if(i==1)
                d[i]=-f.T[i+1]-(-2.0+lamda)*f.T[i];
            else if(i==n-2)
                d[i]=-f.T[i-1]-(-2.0+lamda)*f.T[i];
            else
                d[i]=-f.T[i+1]-(-2.0+lamda)*f.T[i]-f.T[i-1];
            //printf("\ni %d a %f b %f c %f d %f",i,a[i],b[i],c[i],d[i]);
        }

        /*TDMA*/
        for(i=2;i<n-1;i++)
        {
            w[i]=a[i]/b[i-1];
            b[i]=b[i]-w[i]*c[i-1];
            d[i]=d[i]-w[i]*d[i-1];
            a[i]=0;
            //printf("\ni %d w %f a %f b %f c %f d %f",i,w[i],a[i],b[i],c[i],d[i]);
        }
        f.Tnew[n-2]=d[n-2]/b[n-2];
        //printf("\ni %d T %f %f",n-2,f.Tnew[n-2],d[n-2]/b[n-2]);
        for(i=n-3;i>0;i--)
        {
            f.Tnew[i]=(d[i]-c[i]*f.Tnew[i+1])/b[i];
            //printf("\ni %d T %f",i,f.Tnew[i]);
        }
        //printf("\ni %d T %f",i,f.T[i]);

        for(i=0;i<n;i++)
        {
            f.T[i]=f.Tnew[i];
            //printf("\ni %d T %f",i,f.T[i]);
        }
        if((t==(int)ceil(0.01/dt))||(t==(int)ceil(0.02/dt))||(t==(int)ceil(0.03/dt))||(t==(int)ceil(0.04/dt))||(t==(int)ceil(0.05/dt)))
        {
            printf("\nT at time %f and location %f is %f",t*dt,loc*dx,f.T[loc]);
            time[res]=t*dt;
            Tnum[res]=f.T[loc];
            res++;
        }
    }
    /*Comparison between analytical and num*/
    printf("\nt \t\t FDS \t\t AS \t\t Diff \t\t Rel Error (%%)");
    for(i=0;i<6;i++)
    {
        printf("\n%f \t %f \t %f \t %f \t %f",time[i],Tnum[i],Tan[i],fabs(Tnum[i]-Tan[i]),100.0*fabs(Tnum[i]-Tan[i])/Tan[i]);
    }
}

