#include<iostream>
#include<stdlib.h>
#include<cmath> 
#include<sstream>
#include<fstream>
#include<string>
#include<iomanip>
#include<numeric>
#include<functional>
using namespace std;

const double PI  = 3.141592653589793238462;

//int j=0;
int k=0;
int Kmax;
int Kstart;
const int Tmax=10000;
const int Nmax=2; //HOPF LINK - dmi
const int Nbeadsmax=1000;

int N,Nbeads;
int DeltaT=10000;
double dt=0.01;
double Lx,Ly,Lz;

double Rg[Nmax];
double Rg2[Nmax];
double Rg4[Nmax];
double Rgave;
double Rg2ave;
double Rg4ave;
double RAD;
//Observables
double position[Nmax][Nbeadsmax][3];

int main(int argc, char* argv[]){
srand(time(NULL));

cout << "Write argv[1]: datain; argv[2] = dataout; argv[3]=Tmax; argv[4]:Nbeads argv[5]:LD Radius"<< endl;
cout << "Num of files to convert:";
Kmax=int(atoi(argv[3])); 
cout << Kmax << endl;
cout << "Number of Polymers:";
N=1;
cout << N <<endl;
Nbeads=int(atoi(argv[4]));
cout<<"Number of Beads in each polymer:" << Nbeads << endl;
RAD=double(atoi(argv[5]));

// stringstream writeFileKN;
// writeFileKN <<"KN_"<<argv[2];
// ofstream writeKN(writeFileKN.str().c_str());

// stringstream writeFileXYZ;
// writeFileXYZ <<"XYZ_"<<argv[2];
// ofstream writeXYZ(writeFileXYZ.str().c_str());

// stringstream writeFileRG;
// writeFileRG <<"RG_"<<argv[2];
// ofstream writeRG(writeFileRG.str().c_str());

stringstream writeFileLD;
writeFileLD <<"LD_"<<argv[2];
ofstream writeLD(writeFileLD.str().c_str());
writeLD << "time bead_index localdens_chain1_Single localdens_chain2_total localdens_chain1_diff ...chain2..." <<endl; 

stringstream writeFileLDMAX;
writeFileLDMAX <<"LDMAX_"<<argv[2];
ofstream writeLDMAX(writeFileLDMAX.str().c_str());
writeLDMAX << "time bead_index max argmax max_diff argmax_diff" <<endl; 

// stringstream writeFilecurv;
// writeFilecurv <<"CURVATURE_"<<argv[2];
// ofstream writecurv(writeFilecurv.str().c_str());

// stringstream writeFileLCMAX;
// writeFileLCMAX <<"LCMAX_"<<argv[2];
// ofstream writeLCMAX(writeFileLCMAX.str().c_str());

// stringstream writeFilecurvW;
// writeFilecurvW <<"WCURVATURE_"<<argv[2];
// ofstream writecurvW(writeFilecurvW.str().c_str());

// stringstream writeFileLCWMAX;
// writeFileLCWMAX <<"LCWMAX_"<<argv[2];
// ofstream writeLCWMAX(writeFileLCWMAX.str().c_str());

// stringstream writeFileCOR;
// writeFileCOR <<"COR_"<<argv[2];
// ofstream writeCOR(writeFileCOR.str().c_str());

/////////////////////////////////////////////////////
//MAIN LOOP OVER TIME!!!
////////////////////////////////////////////////////
for(k=0;k<Kmax;k++){
ifstream indata;
stringstream readFile;
readFile.str("");
readFile.clear();
long long int sum = int(k);
if(sum==0) readFile << argv[1] <<sum;
if(sum>0)readFile << argv[1] <<sum << "0000";
indata.open(readFile.str().c_str());
cout << readFile.str().c_str()<<endl;
if(!indata){cout <<"file "<< readFile.str().c_str() << " is not there"<<endl; return 0;}

long int id,type,mol;
double num1,num2,num3;
double x,y,z;
string dummy;
long long int time;
long int Ntot;
double l1,l2;

//read 10 lines
for(int i=0;i<10;i++){
if(i==1) {
indata >> time;
time = k*DeltaT*dt;
cout << "time " << time <<endl;
}
if(i==3) {
indata >> Ntot;
cout << "Ntot " << Ntot<<endl;
}
if(i==5) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Lx = l2-l1;
cout << "Lx " << Lx <<endl;
}
if(i==6) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Ly = l2-l1;
cout << "Ly " << Ly <<endl;
}
if(i==7) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Lz = l2-l1;
cout << "Lz " << Lz<<endl;
}

else getline(indata,dummy);
cout << dummy<<endl;
}
 
//////////////////////////
//READ FILES
////////////////////////////
for(int n=0; n<Ntot; n++){

    indata >> id>>type>>x>>y>>z>>num1>>num2>>num3;

    int Nring = floor((double)(id-1)/(1.0*Nbeads));
    
    position[Nring][(id-1)%Nbeads][0] = (x + Lx*num1);
    position[Nring][(id-1)%Nbeads][1] = (y + Ly*num2);
    position[Nring][(id-1)%Nbeads][2] = (z + Lz*num3);

}
//////////////////////////////
//Calculating Decorrelation time
//OR Compiling xyz coordinates
////////////////////////////
// int nr=0;
// //ring
// //writeXYZ << Nbeads+1<< endl;
// for(int n=0;n<Nbeads;n++)writeXYZ << position[nr][n][0] << " " << position[nr][n][1]<<" " << position[nr][n][2] <<endl;
// writeXYZ << endl;

//////////////////////////////
//WRITE KNOT file - for KYMOKNOT
////////////////////////////
// int nr=0;
// //ring
// for(int nr=0;nr<N;nr++){
// writeKN << Nbeads+1<< endl;
// for(int n=0;n<Nbeads;n++)writeKN << position[nr][n][0] << " " << position[nr][n][1]<<" " << position[nr][n][2] <<endl;
// writeKN<<position[nr][0][0]<<" "<<position[nr][0][1]<<" " <<position[nr][0][2]<<endl;
// }

// /////////////////////////////
// // COMPUTE GYRATION RADIUS
// ////////////////////////////
// double rg=0;
// double com[3]={0.0,0.0,0.0};

// for(int n=0;n<Nbeads;n++){
// com[0]+=position[0][n][0];
// com[1]+=position[0][n][1];
// com[2]+=position[0][n][2];
// }
// com[0]/=Nbeads;
// com[1]/=Nbeads;
// com[2]/=Nbeads;
// for(int n=0;n<Nbeads;n++){
// for(int d=0;d<3;d++)rg+=pow(position[0][n][d]-com[d],2.0);
// }
// rg/=Nbeads;


// cout << "Rg " << rg << endl; #this is rg^2

// writeRG << time <<" "<< rg << endl;

/////////////////////////////
// COMPUTE LOCAL DENSITY FUNCTION
//CREATE A TEXT BLOCK WITH 
// time, bead_id, LD_chain1_single, LD_chain1_multi, LD_chain1_difference, LD_chain2_single,  LD_chain2_multi,  LD_chain2_difference
//single = looks at self density
//multi = looks at global density
////////////////////////////
double LDchainsingle[Nmax][Nbeads];
double LDchainmulti[Nmax][Nbeads];
double LDdiff[Nmax][Nbeads];

//parameter for spherical radius - work out what value is optimal
double threshold=RAD;
for(int nn=0;nn<Nmax;nn++){

//calculating the local density around each monomer n to the other monomers j
int max = 0;
int argmax = 0;
int max_diff = 0;
int argmax_diff = 0;
for(int n=0;n<Nbeads;n++){
    int ldsingle = 0;
    int ldmulti = 0;
	for(int jj=0;jj<Nmax;jj++){    
   		for(int j=0;j<Nbeads;j++){
                double eucmulti = sqrt(((position[nn][n][0]-position[jj][j][0])*(position[nn][n][0]-position[jj][j][0])) + ((position[nn][n][1]-position[jj][j][1])*(position[nn][n][1]-position[jj][j][1])) + ((position[nn][n][2]-position[jj][j][2])*(position[nn][n][2]-position[jj][j][2])));
                if (eucmulti<=threshold) ldmulti++;
                           
 		if(jj==nn){
 		double eucsingle = sqrt(((position[nn][n][0]-position[jj][j][0])*(position[nn][n][0]-position[jj][j][0])) + ((position[nn][n][1]-position[jj][j][1])*(position[nn][n][1]-position[jj][j][1])) + ((position[nn][n][2]-position[jj][j][2])*(position[nn][n][2]-position[jj][j][2])));
                if (eucsingle<=threshold) ldsingle++;
                }
    		}
    	}

        if (ldmulti > max){
            max = ldmulti;
            argmax = n;
        }

        int lddiff = ldmulti - ldsingle;
        if (lddiff > max_diff){
            max_diff = lddiff;
            argmax_diff = n;
        }


    
    LDchainsingle[nn][n]=ldsingle;
    LDchainmulti[nn][n]=ldmulti;
    LDdiff[nn][n]=ldmulti-ldsingle;
}

writeLDMAX << k <<" "<< nn << " " << max <<" "<< argmax << " " << max_diff <<" "<< argmax_diff << endl;
}

//write
for(int n=0;n<Nbeads;n++){
writeLD << k <<" "<< n <<" "<< LDchainsingle[0][n] << " " <<LDchainmulti[0][n] << " " << LDdiff[0][n] << " " << LDchainsingle[1][n] << " " << LDchainmulti[1][n] << " " << LDdiff[1][n] << endl;
}
writeLD << endl;
writeLD << endl;

// /////////////////////////////
// // COMPUTE LOCAL CURVATURE FUNCTION
// ////////////////////////////
// double maxlc = 0.0;
// int argmaxlc = 0;
// for(int nn=0;nn<Nmax;nn++){

// for(int n=0;n<Nbeads;n++){
//     //tangents of bead i to bead i-1 and bead i to bead i+1
//     double tangent1[3]={(position[nn][n%Nbeads][nn]-position[nn][(n-1)%Nbeads][0]),(position[nn][n%Nbeads][1]-position[nn][(n-1)%Nbeads][1]),(position[nn][n%Nbeads][2]-position[nn][(n-1)%Nbeads][2])};
//     double tangent2[3]={(position[nn][(n+1)%Nbeads][0]-position[nn][n%Nbeads][0]),(position[nn][(n+1)%Nbeads][1]-position[nn][n%Nbeads][1]),(position[nn][(n+1)%Nbeads][2]-position[nn][n%Nbeads][2])};
//     if(n==0){
//         //account for negative modulus error
//         tangent1[0]=(position[nn][n%Nbeads][0]-position[nn][Nbeads-1][0]);
//         tangent1[1]=(position[nn][n%Nbeads][1]-position[nn][Nbeads-1][1]);
//         tangent1[2]=(position[nn][n%Nbeads][2]-position[nn][Nbeads-1][2]);
//     } 
//     //norm of tangent vectors
//     double tangent1norm = sqrt((tangent1[0]*tangent1[0])+(tangent1[1]*tangent1[1])+(tangent1[2]*tangent1[2]));
//     double tangent2norm = sqrt((tangent2[0]*tangent2[0])+(tangent2[1]*tangent2[1])+(tangent2[2]*tangent2[2]));
//     //dotproduct between tangent vectors
//     double dotprodval = 0;
//     for(int i=0; i<3; i++)dotprodval+=(tangent1[i]*tangent2[i]);
//     //local curvature of bead n
//     double lc = acos(dotprodval/(tangent1norm*tangent2norm));
//     writecurv << k <<" "<< nn << " " <<n <<" "<< lc << endl;
//     if(lc>maxlc) {
//         maxlc = lc;
//         argmaxlc = n;
//     }
// }
// //file with max local curvature per timestep
 // writeCOR << argmaxlc << endl;
// writeLCMAX << k <<" "<< nn << " " << maxlc <<" "<< argmaxlc << endl;
// writecurv << endl;
// writecurv << endl;
// }
// /////////////////////////////////////////////////////////
// // COMPUTE LOCAL CURVATURE FUNCTION WITH WINDOW AVERAGING
// /////////////////////////////////////////////////////////
// double maxlcw = 0.0;
// int argmaxlcw = 0;
// int window = 10; //half of window averaging region
// for(int nn=0;nn<Nmax;nn++){

// for(int n=0;n<Nbeads;n++){
//     double lcwindow = 0;
//     for(int z=n-window;z<n+window;z++){
//         //tangents of bead i to bead i-1 and bead i to bead i+1
//         double tangent1[3]={(position[nn][z%Nbeads][0]-position[nn][(z-1)%Nbeads][0]),(position[nn][z%Nbeads][1]-position[nn][(z-1)%Nbeads][1]),(position[nn][z%Nbeads][2]-position[nn][(z-1)%Nbeads][2])};
//         double tangent2[3]={(position[nn][(z+1)%Nbeads][0]-position[nn][z%Nbeads][0]),(position[nn][(z+1)%Nbeads][1]-position[nn][z%Nbeads][1]),(position[nn][(z+1)%Nbeads][2]-position[nn][z%Nbeads][2])};
//         //account for negative modulus error
//         if(z==0){
//             tangent1[0]=(position[nn][z%Nbeads][0]-position[nn][Nbeads-1][0]);
//             tangent1[1]=(position[nn][z%Nbeads][1]-position[nn][Nbeads-1][1]);
//             tangent1[2]=(position[nn][z%Nbeads][2]-position[nn][Nbeads-1][2]);
//         } 
// 	else if(z==-1){
//             tangent1[0]=(position[nn][Nbeads+z][0]-position[nn][Nbeads+(z-1)][0]);
//             tangent1[1]=(position[nn][Nbeads+z][1]-position[nn][Nbeads+(z-1)][1]);
//             tangent1[2]=(position[nn][Nbeads+z][2]-position[nn][Nbeads+(z-1)][2]);

//             tangent2[0]=(position[nn][(z+1)%Nbeads][0]-position[nn][Nbeads+z][0]);
//             tangent2[1]=(position[nn][(z+1)%Nbeads][1]-position[nn][Nbeads+z][1]);
//             tangent2[2]=(position[nn][(z+1)%Nbeads][2]-position[nn][Nbeads+z][2]);
//         }
//         else if(z<-1){
//             tangent1[0]=(position[nn][Nbeads+z][0]-position[nn][Nbeads+(z-1)][0]);
//             tangent1[1]=(position[nn][Nbeads+z][1]-position[nn][Nbeads+(z-1)][1]);
//             tangent1[2]=(position[nn][Nbeads+z][2]-position[nn][Nbeads+(z-1)][2]);

//             tangent2[0]=(position[nn][Nbeads+(z+1)][0]-position[nn][Nbeads+z][0]);
//             tangent2[1]=(position[nn][Nbeads+(z+1)][1]-position[nn][Nbeads+z][1]);
//             tangent2[2]=(position[nn][Nbeads+(z+1)][2]-position[nn][Nbeads+z][2]);
//         }
//         //norm of tangent vectors
//         double tangent1norm = sqrt((tangent1[0]*tangent1[0])+(tangent1[1]*tangent1[1])+(tangent1[2]*tangent1[2]));
//         double tangent2norm = sqrt((tangent2[0]*tangent2[0])+(tangent2[1]*tangent2[1])+(tangent2[2]*tangent2[2]));
//         //dotproduct between tangent vectors
//         double dotprodval = 0;
//         for(int i=0; i<3; i++)dotprodval+=(tangent1[i]*tangent2[i]);
//         lcwindow += acos(dotprodval/(tangent1norm*tangent2norm));
//     }
//     //local curvature of bead n
//     double lc = lcwindow/(2*window);
//     writecurvW << k <<" "<< n <<" "<< lc << endl;
//     if(lc>maxlcw){
//         maxlcw = lc;
//         argmaxlcw = n;
//     }
// }
//file with max local curvature per timestep
// // writeCOR << argmaxlc << endl;
// writeLCWMAX << k <<" "<< maxlcw <<" "<< argmaxlcw << endl;
// writecurvW << endl;
// writecurvW << endl;
// }

}//closes time (loop over k)

return 0;
}
