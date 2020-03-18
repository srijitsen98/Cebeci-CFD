#include<stdio.h>
#include<math.h>
#define dt 0.001
#define L 1.0
#define r 0.1//alpha*dt/(dx*dx)<0.5
#define alpha 1.0
#define pi 3.14159
int main()
{
    const double delx=sqrt(dt/r); //We keep alpha*dt/dx^2 much smaller than 1 to help stabilize the explicit scheme
    const int n=ceil(L/delx)+1; //Define number of points
    const double dx=L/((n-1)*1.0); //Define the length of each division. Note that if dt is 0.01, dx will come out as some weird value such that n might not be an integer. We try to ensure n is an integer using the previous two steps
    //printf("\n%f %d %f",delx,n,dx);
    double Aw,Ae,Ap;
    double Tnum[5],Tan[5],time[5];
    int i,mid,t,res,loc;
    class system_variables
    {
        public:
            double T[n],x[n],Tnew[n];
    }f;
    /*Analytical Solution stated in book at x=0.3*/
    Tan[0]=.6004;
    Tan[1]=.5966;
    Tan[2]=.5799;
    Tan[3]=.5334;
    Tan[4]=.2444;
    /*Boundary Conditions*/
    f.T[0]=0;
    f.T[n-1]=0;
    /*Initial Conditions*/
    mid=ceil((n-1)/2.0);
    printf("\nMiddle location is %d",mid);
    for(i=1;i<mid;i++)
    {
        f.x[i]=dx*i;
        f.T[i]=2.0*f.x[i];
    }
    for(i=mid;i<n-1;i++)
    {
        f.x[i]=dx*i;
        f.T[i]=2.0*(1-f.x[i]);
    }
    loc=(int)ceil(0.3/dx);
    printf("\nloc is %d",loc);
    res=0;
    Tnum[res]=f.T[loc];
    time[res]=0;
    res++;
    printf("\nT at time %f and location %f is %f",t*dt,loc*dx,f.T[loc]);
    //printf("\nT at 0.3 is %f",f.T[(int)ceil(0.3/dx)]);
    /*Coefficients*/
    Aw=alpha*dt/(dx*dx);
    Ae=alpha*dt/(dx*dx);
    Ap=1-2.0*alpha*dt/(dx*dx);
    /*Solution*/
    for(t=1;t<=1000;t++)
    {
        for(i=1;i<n-1;i++)
        {
            f.Tnew[i]=Aw*f.T[i-1]+Ap*f.T[i]+Ae*f.T[i+1];
        }
        for(i=1;i<n-1;i++)
        {
            f.T[i]=f.Tnew[i];
        }

        if((t==(int)ceil(0.005/dt))||(t==(int)ceil(0.01/dt))||(t==(int)ceil(0.02/dt))||(t==(int)ceil(0.1/dt)))
        {
            printf("\nT at time %f and location %f is %f",t*dt,loc*dx,f.T[loc]);
            time[res]=t*dt;
            Tnum[res]=f.T[loc];
            res++;
        }
    }
    /*Comparison between analytical and num*/
    printf("\nt \t\t FDS \t\t AS \t\t Diff \t\t Error");
    for(i=0;i<5;i++)
    {
        printf("\n%f \t %f \t %f \t %f \t %f",time[i],Tnum[i],Tan[i],fabs(Tnum[i]-Tan[i]),fabs(Tnum[i]-Tan[i])/Tan[i]);
    }


}
