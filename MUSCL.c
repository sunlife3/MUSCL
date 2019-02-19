//Hello World
#include "stdio.h"
#include "math.h"
#define GAMMA 1.4
#define KAPPA 1/3

void mat_multiple(double Z[3][3],double X[3][3],double Y[3][3]){
     for(int i=0;i<3;i++){
         for(int j=0;j<3;j++){
             for(int k=0;k<3;k++){
                 Z[i][j] += X[i][k]*Y[k][j];
             }
         }
     }
}

void A_mat(double A[3][3],double rho,double u,double H,double c){
    double R[3][3]={0},Eigen[3][3]={0},REigen[3][3]={0},Rinv[3][3]={0},b1=0,b2=0;
    b1 = (GAMMA-1)*u*u/(2*c*c);
    b2 = (GAMMA-1)/(c*c);
    for(int i=0;i<3;i++){//initialize A matrix
        for(int j=0;j<3;j++){
            A[i][j] = 0;
        }
    }

    R[0][0] = 1;       R[0][1] = 1;       R[0][2] = 1;
    R[1][0] = u-c;     R[1][1] = u;       R[1][2] = u+c;
    R[2][0] = H-u*c;   R[2][1] = u*u/2;   R[2][2] = H+u*c;

    Eigen[0][0] = fabs(u-c); Eigen[0][1] = 0;        Eigen[0][2] = 0;
    Eigen[1][0] = 0;         Eigen[1][1] = fabs(u);  Eigen[1][2] = 0;
    Eigen[2][0] = 0;         Eigen[2][1] = 0;        Eigen[2][2] = fabs(u+c);

    Rinv[0][0] = (b1+u/c)/2; Rinv[0][1] = -((1/c)+b2*u)/2; Rinv[0][2] = b2/2;
    Rinv[1][0] = 1-b1;       Rinv[1][1] = b2*u;            Rinv[1][2] = -b2;
    Rinv[2][0] = (b1-u/c)/2; Rinv[2][1] =  ((1/c)-b2*u)/2; Rinv[2][2] = b2/2;

    mat_multiple(REigen,R,Eigen);
    mat_multiple(A,REigen,Rinv);
}

int main(){
    double rho[104]={0},u[104]={0},p[104]={0},H[104]={0},e_tot[104]={0};
    double Q1_n[104]={0},Q2_n[104]={0},Q3_n[104]={0};
    double Q1_p[104]={0},Q2_p[104]={0},Q3_p[104]={0};
    double E1[104]={0},E2[104]={0},E3[104]={0};
    double A[3][3]={0};
    double dx=0,dt=0;
    double rho_ave=0,u_ave=0,p_ave=0,H_ave=0,c_ave=0;
    double E_HalfBack[3]={0},E_HalfFor[3]={0};

    dx = 2.0/100;
    dt = 0.001;

    //initial condition
    for(int i=0;i<104;i++){
        if(i<=52){
            rho[i]=1.0;
            p[i]=1.0;
        }else{
            rho[i]=0.125;
            p[i]=0.1;
        } 
         e_tot[i] = p[i]/(GAMMA-1) + rho[i]*u[i]*u[i]/2;
         H[i] = (e_tot[i]+p[i])/rho[i];
         u[i] = 0;
    }
    //calc initial Q vector & E vector
    for(int i=0;i<104;i++){
        Q1_n[i] = rho[i];
        Q2_n[i] = rho[i]*u[i];
        Q3_n[i] = e_tot[i];

        E1[i] = rho[i]*u[i];
        E2[i] = p[i]+rho[i]*u[i]*u[i];
        E3[i] = (GAMMA*p[i]/(GAMMA-1) + rho[i]*u[i]*u[i]/2)*u[i];
    }
    
    for(int k=0;k<500;k++){
        for(int i=2;i<=101;i++){
            double QL[3]={0},rhoL=0,uL=0,pL=0,EL[3]={0};
            double QR[3]={0},rhoR=0,uR=0,pR=0,ER[3]={0};
            //////////////////////////////////////Half Forward////////////////////////////////////////
            //calc (QL)j+1/2
            QL[0] = Q1_n[i] + ((1 - KAPPA)*(Q1_n[i] - Q1_n[i-1]) + (1+KAPPA)*(Q1_n[i+1] - Q1_n[i]))/4;
            QL[1] = Q2_n[i] + ((1 - KAPPA)*(Q2_n[i] - Q2_n[i-1]) + (1+KAPPA)*(Q2_n[i+1] - Q2_n[i]))/4;
            QL[2] = Q3_n[i] + ((1 - KAPPA)*(Q3_n[i] - Q3_n[i-1]) + (1+KAPPA)*(Q3_n[i+1] - Q3_n[i]))/4;
            //calc L state quantity
            rhoL = QL[0];
            uL = QL[1]/rhoL;
            pL = (GAMMA-1)*(QL[2]-rhoL*uL*uL/2);
            //calc (EL)j+1/2
            EL[0] = rhoL*uL;
            EL[1] = pL + rhoL*uL*uL;
            EL[2] = (GAMMA*pL/(GAMMA-1) + rhoL*uL*uL/2)*uL;

            //calc (QR)j+1/2
            QR[0] = Q1_n[i+1] - ((1 - KAPPA)*(Q1_n[i+1] - Q1_n[i]) + (1+KAPPA)*(Q1_n[i] - Q1_n[i-1]))/4;
            QR[1] = Q2_n[i+1] - ((1 - KAPPA)*(Q2_n[i+1] - Q2_n[i]) + (1+KAPPA)*(Q2_n[i] - Q2_n[i-1]))/4;
            QR[2] = Q3_n[i+1] - ((1 - KAPPA)*(Q3_n[i+1] - Q3_n[i]) + (1+KAPPA)*(Q3_n[i] - Q3_n[i-1]))/4;
            //calc R state quantity
            rhoR = QR[0];
            uR = QR[1]/rhoR;
            pR = (GAMMA-1)*(QR[2]-rhoR*uR*uR/2);
            //calc (ER)j+1/2
            ER[0] = rhoL*uR;
            ER[1] = pR + rhoR*uR*uR;
            ER[2] = (GAMMA*pR/(GAMMA-1) + rhoR*uR*uR/2)*uR;

            //calc state average(i,i+1)
            rho_ave = sqrt(rho[i]*rho[i+1]);
            u_ave = (sqrt(rho[i])*u[i] + sqrt(rho[i+1])*u[i+1])/(sqrt(rho[i]) + sqrt(rho[i+1]));
            H_ave = (sqrt(rho[i])*H[i] + sqrt(rho[i+1])*H[i+1])/(sqrt(rho[i]) + sqrt(rho[i+1]));
            c_ave = sqrt((GAMMA-1)*(H_ave - u_ave*u_ave/2));

            //calc A matrix
            A_mat(A,rho_ave,u_ave,H_ave,c_ave);

            //calc E half forward (E(n,j+1/2))
            E_HalfFor[0] = (EL[0] + ER[0] - A[0][0]*(QR[0] - QL[0]) - A[0][1]*(QR[1] - QL[1]) - A[0][2]*(QR[2] - QL[2]))/2;
            E_HalfFor[1] = (EL[1] + ER[1] - A[1][0]*(QR[0] - QL[0]) - A[1][1]*(QR[1] - QL[1]) - A[1][2]*(QR[2] - QL[2]))/2;
            E_HalfFor[2] = (EL[2] + ER[2] - A[2][0]*(QR[0] - QL[0]) - A[2][1]*(QR[1] - QL[1]) - A[2][2]*(QR[2] - QL[2]))/2;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////Half Back ///////////////////////////////////////////////////////////
            QL[0] = Q1_n[i-1] + ((1 - KAPPA)*(Q1_n[i-1] - Q1_n[i-2]) + (1+KAPPA)*(Q1_n[i] - Q1_n[i-1]))/4;
            QL[1] = Q2_n[i-1] + ((1 - KAPPA)*(Q2_n[i-1] - Q2_n[i-2]) + (1+KAPPA)*(Q2_n[i] - Q2_n[i-1]))/4;
            QL[2] = Q3_n[i-1] + ((1 - KAPPA)*(Q3_n[i-1] - Q3_n[i-2]) + (1+KAPPA)*(Q3_n[i] - Q3_n[i-1]))/4;
            //calc L state quantity
            rhoL = QL[0];
            uL = QL[1]/rhoL;
            pL = (GAMMA-1)*(QL[2]-rhoL*uL*uL/2);
            //calc (EL)j-1/2
            EL[0] = rhoL*uL;
            EL[1] = pL + rhoL*uL*uL;
            EL[2] = (GAMMA*pL/(GAMMA-1) + rhoL*uL*uL/2)*uL;

            //calc (QR)j-1/2
            QR[0] = Q1_n[i] - ((1 - KAPPA)*(Q1_n[i] - Q1_n[i-1]) + (1+KAPPA)*(Q1_n[i-1] - Q1_n[i-2]))/4;
            QR[1] = Q2_n[i] - ((1 - KAPPA)*(Q2_n[i] - Q2_n[i-1]) + (1+KAPPA)*(Q2_n[i-1] - Q2_n[i-2]))/4;
            QR[2] = Q3_n[i] - ((1 - KAPPA)*(Q3_n[i] - Q3_n[i-1]) + (1+KAPPA)*(Q3_n[i-1] - Q3_n[i-2]))/4;
            //calc R state quantity
            rhoR = QR[0];
            uR = QR[1]/rhoR;
            pR = (GAMMA-1)*(QR[2]-rhoR*uR*uR/2);
            //calc (ER)j-1/2
            ER[0] = rhoL*uR;
            ER[1] = pR + rhoR*uR*uR;
            ER[2] = (GAMMA*pR/(GAMMA-1) + rhoR*uR*uR/2)*uR;

            //calc state average(i-1,i)
            rho_ave = sqrt(rho[i]*rho[i-1]);
            u_ave = (sqrt(rho[i])*u[i] + sqrt(rho[i-1])*u[i-1])/(sqrt(rho[i]) + sqrt(rho[i-1]));
            H_ave = (sqrt(rho[i])*H[i] + sqrt(rho[i-1])*H[i-1])/(sqrt(rho[i]) + sqrt(rho[i-1]));
            c_ave = sqrt((GAMMA-1)*(H_ave - u_ave*u_ave/2));

            //calc A matrix       
            A_mat(A,rho_ave,u_ave,H_ave,c_ave);

            //calc E half Back (E(n,j-1/2))
            E_HalfBack[0] = (EL[0] + ER[0] - A[0][0]*(QR[0]-QL[0]) - A[0][1]*(QR[1] - QL[1]) - A[0][2]*(QR[2] - QL[2]))/2;
            E_HalfBack[1] = (EL[1] + ER[1] - A[1][0]*(QR[0]-QL[0]) - A[1][1]*(QR[1] - QL[1]) - A[1][2]*(QR[2] - QL[2]))/2;
            E_HalfBack[2] = (EL[2] + ER[2] - A[2][0]*(QR[0]-QL[0]) - A[2][1]*(QR[1] - QL[1]) - A[2][2]*(QR[2] - QL[2]))/2;
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //calc Q_n+1
            Q1_p[i] = Q1_n[i] - dt*(E_HalfFor[0] - E_HalfBack[0])/dx;
            Q2_p[i] = Q2_n[i] - dt*(E_HalfFor[1] - E_HalfBack[1])/dx;
            Q3_p[i] = Q3_n[i] - dt*(E_HalfFor[2] - E_HalfBack[2])/dx;
        }
        
        for(int i=2;i<=101;i++){
            //update state quantity
            rho[i] = Q1_p[i];
            u[i] = Q2_p[i]/rho[i];
            p[i] = (GAMMA-1)*(Q3_p[i] - rho[i]*u[i]*u[i]/2);
            
            e_tot[i] = Q3_p[i];
            H[i] = (e_tot[i]+p[i])/rho[i];

            Q1_n[i] = Q1_p[i];
            Q2_n[i] = Q2_p[i];
            Q3_n[i] = Q3_p[i];

            E1[i] = rho[i]*u[i];
            E2[i] = p[i] + rho[i]*u[i]*u[i];
            E3[i] = (GAMMA*p[i]/(GAMMA-1) + rho[i]*u[i]*u[i]/2)*u[i];    
        }
    }
    
    FILE *fp;
    fp = fopen("MUSCL.csv","w");
    for(int i=2;i<102;i++){
        fprintf(fp,"%f,%f,%f,%f\n",(i-2)*dx-1.0,rho[i],u[i],p[i]);
    }
    fclose(fp);
    return 0;
}   
