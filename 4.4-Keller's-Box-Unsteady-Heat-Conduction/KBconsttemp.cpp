#include<stdio.h>
#include<math.h>
#define dx 0.1 //spacing in x direction
#define L 1 //total length
#define dt 0.01 //time spacing
#define pi 3.14159
double b[2][2],c[2]; //variables that are used to return values after matrix operations
int inv(double a[2][2]) //matrix inversion function
{
    double det;
    det=a[0][0]*a[1][1]-a[0][1]*a[1][0];
    b[0][0]=a[1][1]/(det*1.0);
    b[0][1]=-a[0][1]/(det*1.0);
    b[1][0]=-a[1][0]/(det*1.0);
    b[1][1]=a[0][0]/(det*1.0);
    return 0;
}
int multiply2x2(double x[2][2],double y[2][2]) //matrix multiplication function
{
    b[0][0]=x[0][0]*y[0][0]+x[0][1]*y[1][0];
    b[0][1]=x[0][0]*y[0][1]+x[0][1]*y[1][1];
    b[1][0]=x[1][0]*y[0][0]+x[1][1]*y[1][0];
    b[1][1]=x[1][0]*y[0][1]+x[1][1]*y[1][1];
    return 0;
}
int subtract2x2(double x[2][2],double y[2][2]) //matrix subtraction function
{
    b[0][0]=x[0][0]-y[0][0];
    b[0][1]=x[0][1]-y[0][1];
    b[1][0]=x[1][0]-y[1][0];
    b[1][1]=x[1][1]-y[1][1];
    return 0;
}
int multiply2x1(double x[2][2],double z[2]) //matrix multiplication function
{
    c[0]=x[0][0]*z[0]+x[0][1]*z[1];
    c[1]=x[1][0]*z[0]+x[1][1]*z[1];
    return 0;
}
int subtract2x1(double x[2],double y[2]) //matrix subtraction function
{
    c[0]=x[0]-y[0];
    c[1]=x[1]-y[1];
    return 0;
}
int main()
{
    //double T[n];
    const int nx=(int)L/dx;
    const int nt=(int)(0.1/dt);
    int i,mid,j,res,t,loc;
    printf("\nValue of nx is %d and nt is %d",nx,nt);
    double Tnew[nx+1],Told[nx+1],A[nx+1][2][2],B[nx+1][2][2],C[nx+1][2][2],r[nx+1][2],del[nx+1][2][2],w[nx+1][2],tau[nx+1][2][2],d[nx+1][2],delinv[nx+1][2][2];
    double s1, s2, s3, lamda, location,err;
    double alpha=1;
    double temp1[2][2],temp2[2];
    double pold[nx+1],pnew[nx+1];
    double Tw=0;
    double Te=0;
    double Tnum[12],time[12],Tan[12];
    /*Analytical Solution stated in book at x=0.3*/
    location=0.3;
    Tan[0]=.6004;
    Tan[1]=.5799;
    Tan[2]=.5334;
    Tan[3]=.4857;
    Tan[4]=.4411;
    Tan[5]=.4000;
    Tan[6]=.3626;
    Tan[7]=.3286;
    Tan[8]=.2977;
    Tan[9]=.2698;
    Tan[10]=.2444;
    loc=(int)ceil(location/dx);
    printf("\nLocation where analytical and numerical solution is compared = %f",loc*dx);
    printf("\nArray location is %d",loc);
    res=0; //This variable gets updated with time. Suppose you wanna compare analytical solution at x=0.3 at t=0.03s, res would be 3 as the analytical solution corresponding to t=0.03s is represented by Tan[3]
    lamda=dx/(alpha*dt);
    s1=1;
    s2=-1;
    s3=-lamda;
    //Assigning values to A,B and C
    //Left Boundary
    A[0][0][0]=1;
    A[0][0][1]=0;
    A[0][1][0]=-1;
    A[0][1][1]=-dx/2.0;
    C[0][0][0]=0;
    C[0][0][1]=0;
    C[0][1][0]=1;
    C[0][1][1]=-dx/2.0;
    //internal cells
    for(i=1;i<nx;i++)
    {
        B[i][0][0]=s3;
        B[i][0][1]=s2;
        B[i][1][0]=0;
        B[i][1][1]=0;
        A[i][0][0]=s3;
        A[i][0][1]=s1;
        A[i][1][0]=-1;
        A[i][1][1]=-dx/2.0;
        C[i][0][0]=0;
        C[i][0][1]=0;
        C[i][1][0]=1;
        C[i][1][1]=-dx/2.0;
    }
    //Right boundary
    A[nx][0][0]=s3;
    A[nx][0][1]=s1;
    A[nx][1][0]=1;
    A[nx][1][1]=0;
    B[nx][0][0]=s3;
    B[nx][0][1]=s2;
    B[nx][1][0]=0;
    B[nx][1][1]=0;
    //Initialization of values at t=0
    mid=(int)nx/2.0;    //midpoint
    for(i=0;i<=mid;i++)
    {
        Told[i]=2*i*dx;
        Tnew[i]=Told[i];
        pold[i]=(2*(i+1)*dx-Told[i])/dx;    //just a rough approximation to start of with. p is supposed to be computed in between the cells but to start off, I computed using forward difference for the first few cells and backward difference for the other cells
        pnew[i]=pold[i];
    }
    for(i=mid+1;i<=nx;i++)
    {
        Told[i]=2*(1-i*dx);
        Tnew[i]=Told[i];
        pold[i]=(Told[i]-Told[i-1])/dx;
        pnew[i]=pold[i];
    }
    Tnum[res]=Told[loc];
    time[res]=0;
    res++;
    j=0;
    err=fabs(Tnew[loc]-Tan[res])*100.0/Tnew[loc]; //relative error between analytical and numerical solution
    printf("\nT at time %f and location %f is %f, error is %f %%",j*dt,loc*dx,Tnew[loc],err);
    for(j=1;j<=nt;j++)
    {
        for(i=0;i<=nx;i++)
        {
            Told[i]=Tnew[i];
            pold[i]=pnew[i];
        }
        //Assign r
        //Left boundary
        r[0][0]=Tw;
        r[0][1]=0;
        //central nodes
        for(i=1;i<nx;i++)
        {
            r[i][0]=-lamda*(Told[i]+Told[i-1])+pold[i-1]-pold[i];
            r[i][1]=0;
        }
        //Right boundary
        r[nx][0]=-lamda*(Told[nx]+Told[nx-1])+pold[nx-1]-pold[nx];
        r[nx][1]=Te;
        //Define del
        del[0][0][0]=A[0][0][0];
        del[0][0][1]=A[0][0][1];
        del[0][1][0]=A[0][1][0];
        del[0][1][1]=A[0][1][1];
        w[0][0]=r[0][0];
        w[0][1]=r[0][1];
        for(i=1;i<=nx;i++)
        {
            //Calculate Tau
            inv(del[i-1]);
            delinv[i-1][0][0]=b[0][0];
            delinv[i-1][0][1]=b[0][1];
            delinv[i-1][1][0]=b[1][0];
            delinv[i-1][1][1]=b[1][1];
            multiply2x2(B[i],delinv[i-1]);
            tau[i][0][0]=b[0][0];
            tau[i][0][1]=b[0][1];
            tau[i][1][0]=b[1][0];
            tau[i][1][1]=b[1][1];
            multiply2x2(tau[i],C[i-1]);
            temp1[0][0]=b[0][0];
            temp1[0][1]=b[0][1];
            temp1[1][0]=b[1][0];
            temp1[1][1]=b[1][1];
            subtract2x2(A[i],temp1);
            del[i][0][0]=b[0][0];
            del[i][0][1]=b[0][1];
            del[i][1][0]=b[1][0];
            del[i][1][1]=b[1][1];
            multiply2x1(tau[i],w[i-1]);
            temp2[0]=c[0];
            temp2[1]=c[1];
            subtract2x1(r[i],temp2);
            w[i][0]=c[0];
            w[i][1]=c[1];
        }
        inv(del[nx]);
        delinv[nx][0][0]=b[0][0];
        delinv[nx][0][1]=b[0][1];
        delinv[nx][1][0]=b[1][0];
        delinv[nx][1][1]=b[1][1];

        multiply2x1(delinv[nx],w[nx]);
        d[nx][0]=c[0];
        d[nx][1]=c[1];
        for(i=nx-1;i>=0;i--)
        {
            multiply2x1(C[i],d[i+1]);
            temp2[0]=c[0];
            temp2[1]=c[1];
            subtract2x1(w[i],temp2);
            temp2[0]=c[0];
            temp2[1]=c[1];
            multiply2x1(delinv[i],temp2);
            d[i][0]=c[0];
            d[i][1]=c[1];
            Tnew[i]=d[i][0];
            pnew[i]=d[i][1];
        }
        err=fabs(Tnew[loc]-Tan[res])*100.0/Tnew[loc];
        printf("\nT at time %f and location %f is %f, error is %f %%",j*dt,loc*dx,Tnew[loc],err);
        time[res]=t*dt;
        Tnum[res]=Tnew[loc];
        res++;
    }
}

