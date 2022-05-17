using namespace std;

#include<iostream>
#include<stdlib.h>
#include<cmath>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
#include<iomanip>
#include <string.h>
#include "tools.h"

int id,type,Ntot;
double Lminx,Lmaxx,Lminy,Lmaxy,Lminz,Lmaxz;
int ix,iy,iz;
double x,y,z;
double Lx,Ly,Lz;

long int k=0;
long int Kmax;
int Kstart;
const int Nmaxx=10;
const int Npolymax=2000;

int Nmax,Nbeads;
vector< vector< vector<double>>> atoms(Nmaxx, vector<vector<double> >(Npolymax, vector<double>(3)));
atom *atoms_self;
int DeltaT=10000;
double Wr=0.0;

string dummy; 

double Polymer[Nmaxx][Npolymax][3];
double Polymer0[Nmaxx][Npolymax][3];
double COM[Nmaxx][3];
double COM0[Nmaxx][3];
int lp;

double dotProduct(vector<double> A, vector<double> B);

vector<double> crossProduct(vector<double> A, vector<double> B);

double lengthVector(vector<double> A);

vector<double> makeUnit(vector<double> A);

vector<double> subtraction(vector<double> A, vector<double> B);

double ComputeLocalUnSignedWrithe_self(int n1, int n2, int n3);
double ComputeLocalUnSignedWrithe_total(int n1, int n2, int n3);

int main(int argc, char* argv[]){

cout << "Write argv[1]: datain -- \nargv[2] = dataout no extension \targv[3]=max # of frames -- argv[4]=nbeads; argv[5]= nchains argv[6]=window size"<< endl;
    
Kmax=atoi(argv[3]);
//
//
Nbeads=atoi(argv[4]);
Nmax=atoi(argv[5]);
//
//
lp=atof(argv[6]);
//
//OUTFILE
//
stringstream writeFileWr;
writeFileWr <<"Writhe_" << argv[2]<< ".lp"<<lp<<".dat";
ofstream writeWr(writeFileWr.str().c_str());
writeWr << "#Timestep BeadNumber 3DLocalWrithe_self_chain1 3DLocalWrithe_total_chain1 3DLocalWrithe_diff_chain1 3DLocalWrithe_self_chain2 3DLocalWrithe_total_chain2 3DLocalWrithe_diff_chain2 "<< endl;

stringstream writeFileWrmax;
writeFileWrmax <<"WritheMAX_" << argv[2]<< ".lp"<<lp<<".dat";
ofstream writeWrmax(writeFileWrmax.str().c_str());
writeWrmax << "#Timestep Chain Max3DLocalWrithe ArgMax3DLocalWrithe Max3DLocalWrithe_diff ArgMax3DLocalWrithe_diff" << endl;

//
/////////////////////////////////////////////////////
//MAIN LOOP OVER TIME!!!
////////////////////////////////////////////////////
for(k=0;k<Kmax;k++){
ifstream indata;
stringstream readFile;
readFile.str("");
readFile.clear();
long int timet;
long long int sum = int(k);
if(k==0)readFile << argv[1] <<sum;
//had to fix this as previously forgot "sum" and "" around 0000, so was only reading first config repeatedly - joseph
else readFile << argv[1] <<sum <<"0000";
indata.open(readFile.str().c_str());
cout << readFile.str().c_str()<<endl;
if(!indata.good()){cout << "INCORRECT INPUT FILE NAME"<<endl; cin.get();cin.get();}

//read 10 lines
for(int i=0;i<10;i++){
if(i==1) {
indata >> timet;
timet = k*DeltaT/10000.0;
cout << "time " << timet <<endl;
}
if(i==3) {
indata >> Ntot;
cout << "Ntot " << Ntot<<endl;
}
if(i==5) {
indata >> Lminx >> Lmaxx;
cout << "L " << Lminx<< " " << Lmaxx <<endl;
Lx = Lmaxx-Lminx;
cout << "Lx " << Lx <<endl;
}
if(i==6) {
indata >> Lminy >> Lmaxy;
cout << "L " << Lminy << " " << Lmaxy <<endl;
Ly = Lmaxy-Lminy;
cout << "Ly " << Ly <<endl; 
}
if(i==7) {
indata >> Lminz >> Lmaxz;
cout << "L " << Lminz << " " << Lmaxz <<endl;
Lz = Lmaxz-Lminz;
cout << "Lz " << Lz <<endl;
}

else getline(indata,dummy);
cout << dummy<<endl;
}

//////////////////////////
//READ FILES
////////////////////////////
int id,mol,type;
cout << "reading .. "<<endl;
for(int n=0; n<Ntot; n++){

	indata >> id>> type >> x>>y>>z>>ix>>iy>>iz;
        int npoly = floor((double)(id-1)/(double)(Nbeads));
	int nb = (id-1)%Nbeads;

	Polymer[npoly][nb][0]=x+Lx*ix;
	Polymer[npoly][nb][1]=y+Ly*iy;
	Polymer[npoly][nb][2]=z+Lz*iz;
	//COM
	COM[npoly][0]+=Polymer[npoly][nb][0]*1.0/Nbeads;
	COM[npoly][1]+=Polymer[npoly][nb][1]*1.0/Nbeads;
	COM[npoly][2]+=Polymer[npoly][nb][2]*1.0/Nbeads;
	
	if(k==0){
	Polymer0[npoly][nb][0]=Polymer[npoly][nb][0];
	Polymer0[npoly][nb][1]=Polymer[npoly][nb][1];
	Polymer0[npoly][nb][2]=Polymer[npoly][nb][2];
	}
}
//set initial points
if(k==0){
for(int n=0; n<Nmax; n++){
COM0[n][0]=COM[n][0];
COM0[n][1]=COM[n][1];
COM0[n][2]=COM[n][2];
}
}
cout << "done "<<endl;

//////////////
// Compute 3D Writhe
/////////////////
cout << "Computing 3D Writhe " <<endl;
double localWr_self[Nmax][Nbeads];
double localWr_total[Nmax][Nbeads];
double localWr_diff[Nmax][Nbeads];

for(int nn=0;nn<Nmax;nn++){ // iterating over polymer chains

    double wrmax=0;
    int wrargmax=0;
    double wrmax_diff=0;
    int wrargmax_diff=0;

    cout << "doing polymer .. " << nn << endl;
    for(int kk=0;kk<Nbeads;kk++){ // iterating over beads
        double wr_self, wr_total, wr_diff;

        wr_self = ComputeLocalUnSignedWrithe_self(nn,kk,lp);
        wr_total = ComputeLocalUnSignedWrithe_total(nn,kk,lp);
        wr_diff = wr_total - wr_self;

        localWr_self[nn][kk] = wr_self;
        localWr_total[nn][kk] = wr_total;
        localWr_diff[nn][kk] = wr_diff;

        if (wr_total > wrmax){
            wrmax = wr_total;
            wrargmax = kk;
        }

        if (wr_diff > wrmax_diff){
            wrmax_diff = wr_diff;
            wrargmax_diff = kk;
        }
    }
    writeWrmax << k  << " " << nn << " " << wrmax << " " << wrargmax << " " << wrmax_diff << " " << wrargmax_diff << endl;
}

//write
for(int n=0;n<Nbeads;n++){
    writeWr << k << " " << n << " " << localWr_self[0][n] << " " << localWr_total[0][n] << " " << localWr_diff[0][n] << " " << localWr_self[1][n] << " " << localWr_total[1][n] << " " << localWr_diff[1][n] << endl;
}
writeWr << "\n\n" ;
} //close loop over read/time
return 0;
} // close main function

double ComputeLocalUnSignedWrithe_self(int r, int k, int lp){
    atoms_self = new atom[Nbeads];
    for (int i=0;i<Nbeads;i++) {
        atoms_self[i].x=Polymer[r][i][0];
        atoms_self[i].y=Polymer[r][i][1];
        atoms_self[i].z=Polymer[r][i][2];
    }
    // calculate writhe for segment k-lp/2 < i < k+lp/2 over entire contour - local method of non-local entanglement
    double Wr=0.0;
    xyz one,two, three, four,
    r12,r34,r23,r13,r14,r24,
    n1,n2,n3,n4,cvec;
    double omega, n1n2, n2n3, n3n4, n4n1;
    int i1=k-(lp/2);int i2=k+(lp/2);
    
    for (int i=i1;i<i2;i++) {  //loop over segments within k-lp/2 < i < k+lp/2
        int ib=i;
        if(i<0)ib=i+Nbeads;
        if(i>=Nbeads)ib=i%Nbeads;
        
        int j1=0; int j2=Nbeads;
        for (int j=j1;j<j2;j++) {  //loop over entire contour 0 < j < Nbeads
            int jb=j;
            if(j<0)jb=j+Nbeads;
            if(j>=Nbeads)jb=j%Nbeads;
            
            if (ib==jb||jb==ib-1||jb==ib+1||(ib==0&&jb==Nbeads-1)||(ib==Nbeads-1&&jb==0)) {
                omega=0.0;
            } else {
                //
                if (jb==0) {
                    three=atoms_self[Nbeads-1];//D[L-1];
                } else{
                    three=atoms_self[jb-1];
                }
                if(ib==0){
                    one=atoms_self[Nbeads-1];
                }else{
                    one=atoms_self[ib-1];
                }
                
                two=atoms_self[ib]; four=atoms_self[jb];
                r12=two-one;
                r34=four-three;
                r23=three-two;
                r13=three-one;
                r14=four-one;
                r24=four-two;
                
                n1=r13.cross(r14); //if (n1.length()<SMALL) {cout<<"error in writhe 1 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n1.make_unit();
                
                n2=r14.cross(r24); //if (n2.length()<SMALL) {cout<<"error in writhe 2 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n2.make_unit();
                
                n3=r24.cross(r23); //if (n3.length()<SMALL) {cout<<"error in writhe 3 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n3.make_unit();
                
                n4=r23.cross(r13); //if (n4.length()<SMALL) {cout<<"error in writhe 4 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n4.make_unit();
                
                n1n2=n1.dot(n2); fix_roundoff11(n1n2);
                n2n3=n2.dot(n3); fix_roundoff11(n2n3);
                n3n4=n3.dot(n4); fix_roundoff11(n3n4);
                n4n1=n4.dot(n1); fix_roundoff11(n4n1);
                
                cvec=r34.cross(r12);
                
                omega = ( asin( n1n2 ) + asin( n2n3 ) + asin( n3n4 ) + asin( n4n1 ) )*1;
                
            }
            Wr+= omega/(4.0*PI);
        } //loop over j
    } //loop over beads i
    Wr*=2.0;
    delete [] atoms_self;
    return Wr;
}

double ComputeLocalUnSignedWrithe_total(int nn, int k, int lp){
    // atoms = new atom[Nmaxx][Npolymax];

    for (int j=0;j<Nmax;j++) { // loop over polymers
    	for (int i=0;i<Nbeads;i++) { // loop over beads
        atoms[j][i][0]=Polymer[j][i][0];
        atoms[j][i][1]=Polymer[j][i][1];
        atoms[j][i][2]=Polymer[j][i][2];
    	}
    }

    double Wr=0.0;
    vector<double> one, two, three, four,
    r12,r34,r23,r13,r14,r24,
    n1,n2,n3,n4,cvec;
    double omega, n1n2, n2n3, n3n4, n4n1;
    
    
    for (int jj=0;jj<Nmax;jj++){ // loop over polymer chains

    int i1=k-(lp/2);int i2=k+(lp/2);
    
    for (int i=i1;i<i2;i++) {  //loop over segments within k-lp/2 < i < k+lp/2
        int ib=i;
        if(i<0)ib=i+Nbeads;
        if(i>=Nbeads)ib=i%Nbeads;
        
        int j1=0; int j2=Nbeads;
        for (int j=j1;j<j2;j++) {  //loop over entire contour
            int jb=j;
            if(j<0)jb=j+Nbeads;
            if(j>=Nbeads)jb=j%Nbeads;
            
            if (ib==jb||jb==ib-1||jb==ib+1||(ib==0&&jb==Nbeads-1)||(ib==Nbeads-1&&jb==0)) {
                omega=0.0;
            } else {
                //
                ////MAIN BIT IM UNSURE ABOUT - joseph ///////////
                if (jb==0) {
                    three=atoms[jj][Nbeads-1];//D[L-1];
                } else{
                    three=atoms[jj][jb-1];
                }
                if(ib==0){
                    one=atoms[nn][Nbeads-1];
                }else{
                    one=atoms[nn][ib-1];
                }
                
                two=atoms[nn][ib]; four=atoms[jj][jb];
                ////////////////////////////////////////////
                r12=subtraction(two,one);
                r34=subtraction(four,three);
                r23=subtraction(three,two);
                r13=subtraction(three,one);
                r14=subtraction(four,one);
                r24=subtraction(four,two);

                n1=crossProduct(r13,r14); //if (n1.length()<SMALL) {cout<<"error in writhe 1 "<<n1.length()<<" ("<<i<<")"<<endl;}

                n1=makeUnit(n1);
                
                n2=crossProduct(r14,r24); //if (n2.length()<SMALL) {cout<<"error in writhe 2 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n2=makeUnit(n2);
                
                n3=crossProduct(r24,r23); //if (n3.length()<SMALL) {cout<<"error in writhe 3 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n3=makeUnit(n3);
                
                n4=crossProduct(r23,r13); //if (n4.length()<SMALL) {cout<<"error in writhe 4 "<<n1.length()<<" ("<<i<<")"<<endl;}
                n4=makeUnit(n4);

                n1n2=dotProduct(n1, n2); fix_roundoff11(n1n2);
                // exit(EXIT_FAILURE);
                n2n3=dotProduct(n2, n3); fix_roundoff11(n2n3);
                n3n4=dotProduct(n3, n4); fix_roundoff11(n3n4);
                n4n1=dotProduct(n4, n1); fix_roundoff11(n4n1);
                cvec=crossProduct(r34,r12);
                omega = ( asin( n1n2 ) + asin( n2n3 ) + asin( n3n4 ) + asin( n4n1 ) )*1;
                
            }
            Wr+= omega/(4.0*PI);
        } //loop over j
    } //loop over beads i

    }

    Wr*=2.0;
    // cout << k << " " << nn << endl;
    // delete [] atoms;
    return Wr;
}

double dotProduct(vector<double> A, vector<double> B){
    double p=0.0;
    for(int i=0; i<3;i++){
        p+=A[i]*B[i];
    }
    return p;
}

vector<double> crossProduct(vector<double> A, vector<double> B){
    vector<double> C(3);
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
    return C;
}

double lengthVector(vector<double> A){ // return length of vector
    return sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
}

vector<double> makeUnit(vector<double>A){
    double l;
    l = lengthVector(A);
    A[0]/=l;
    A[1]/=l;
    A[2]/=l;
    return A;
}

vector<double> subtraction(vector<double> A, vector<double> B){
    vector<double> C(3);
    C[0] = A[0] - B[0];
    C[1] = A[1] - B[1];
    C[2] = A[2] - B[2];
    return C;
}