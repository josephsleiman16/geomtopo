#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector> 
#include <unistd.h>
#include<ctime>
using namespace std;

int main(int argc, char* argv[]){

cout << "This program generates 2 rings of Nbeads length" <<endl;
cout << "Type 1. output name; 2. Nbeads" <<endl;

ofstream write;
write.open(argv[1]);
int N=2;
int Nbeads=atoi(argv[2]);

double pos_k1[N][Nbeads][3];

write <<"LAMMPS data file from restart file: timestep = 0, procs = 1"<<endl;
write <<endl;

write<< 1*N*Nbeads << " atoms"<<endl;
write<< 1*N*Nbeads << " bonds"<<endl;
write<< 1*N*Nbeads << " angles"<<endl;
write <<endl;

write << 4 << " atom types" <<endl;
write << 1 << " bond types" <<endl;
write << 1 << " angle types" <<endl;
write <<endl;

int Lx=100;
int Ly=100;
int Lz=100;
write << -Lx/2.0 << " " << Lx/2.0  << " xlo xhi "<<endl;
write << -Ly/2.0 << " " << Ly/2.0 << " ylo yhi "<<endl;
write << -Lz/2.0 << " " << Lz/2.0 << " zlo zhi "<<endl;
write <<endl;

write << "Masses "<<endl;
write <<endl;
write << " 1 1 " <<endl;
write << " 2 1 " <<endl;
write << " 3 1 " <<endl;
write << " 4 1 " <<endl;
write <<endl;

write << "Atoms "<<endl;
write <<endl;

int n;
for(n=0;n<N;n++){
    double shiftz= 0; //rand()*1.0/RAND_MAX*Lz;
    for(int m=0; m<Nbeads; m++){

    double t = 2.0*M_PI*m/(double)(Nbeads);
    double r_out=Nbeads*1.0/10.0;
        
    if(n==0){
    pos_k1[n][m][0]=r_out*cos(t);
    pos_k1[n][m][1]=r_out*sin(t);
    pos_k1[n][m][2]=shiftz;
    }
        
    if(n==1){
       pos_k1[n][m][0]=0.0;
       pos_k1[n][m][1]=-r_out+r_out*cos(t);
       pos_k1[n][m][2]=r_out*sin(t);
       }
      
        //NEED TO WRITE:
    //INDEX MOLECULE TYPE X Y Z IX IY IZ
    write << n*Nbeads+m+1 << " " << n+1 <<" " << n+1 << " " << pos_k1[n][m][0] << " " <<  pos_k1[n][m][1] << " "  << pos_k1[n][m][2] << " 0 0 0 " <<endl;
        //cout << m <<endl;
}
}

    write <<endl;
    write <<endl;
//CREATE BONDS BETWEEN BEADS
//THIS IS A LINEAR POLYMER SO N-1 BONDS
int nbonds=1;
write << "Bonds"<<endl;
write <<endl;
for(int n=0;n<N;n++){
    for(int m=0;m<Nbeads;m++){
        if(m<Nbeads-1)write << nbonds << " " << 1  << " " << n*Nbeads+m+1 << " " << n*Nbeads+m+2<<endl;
        if(m==Nbeads-1)write << nbonds << " " << 1  << " " << n*Nbeads+m+1 << " " << n*Nbeads+1<<endl;
nbonds++;
}
}
    
write <<endl;

//CREATE ANGLES BETWEEN BEADS
//THIS IS A LINEAR POLYMER SO N-2 ANGLES
int nangles=1;
write << "Angles"<<endl;
write <<endl;
for(int n=0;n<N;n++){
    for(int m=0;m<Nbeads;m++){
if(m<Nbeads-2)write << nangles << " " << 1  << " " << n*Nbeads+m+1<< " " << n*Nbeads+m+2 << " " << n*Nbeads+m+3<<endl;
if(m==Nbeads-2)write << nangles << " " << 1  << " " << n*Nbeads+m+1 << " " << n*Nbeads+m+2 << " " << n*Nbeads+1<<endl;
if(m==Nbeads-1)write << nangles << " " << 1  << " " << n*Nbeads+m+1 << " " << n*Nbeads+1 << " " << n*Nbeads+2<<endl;
nangles++;
}
}

return 0;
}


