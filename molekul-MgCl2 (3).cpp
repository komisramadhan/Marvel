%%writefile molekul-MgCl2.cpp
#include "iostream"
#include "cmath"
#include "iomanip"
#include "fstream"

using namespace std;

int main(){
    // deklarasi variabel
    float pi = 3.14;
    
    // deklarasi molekul air
    float Mgx, Mgy, Mgz;
    float Cl1x, Cl1y, Cl1z;
    float Cl2x, Cl2y, Cl2z;

    // masukkan nilai untuk setiap variabel MgCl2
    Mgx = Mgy = Mgz = 0.0;
    Cl1x = sin(90*pi/180) * 2.350;
    Cl1y = 0.0;
    Cl1z = cos(90*pi/180) * 2.350;

    Cl2x = - Cl1x;
    Cl2y = Cl1y;
    Cl2z = Cl1z;

    // memasukkan nilai tersebut ke dalam file
    ofstream file;
    file.open("MgCl2.xyz")
    file << "3\n" << endl;
    file << setw(3) << "Mg" << setw(3) << " " 
    << fixed << setprecision(3) << Mgx << setw(3) << " " \
    << fixed << setprecision(3) << Mgy << setw(3) << " " \
    << fixed << setprecision(3) << Mgz << "\n";
    

    file << setw(3) << "Cl1" << setw(3) << " " \
    << fixed << setprecision(3) << Cl1x << setw(3) << " " \
    << fixed << setprecision(3) << Cl1y << setw(3) << " " \
    << fixed << setprecision(3) << Cl1z << "\n";

    file << setw(3) << "Cl2" << setw(3) << " " \
    << fixed << setprecision(3) << Cl1x << setw(3) << " " \
    << fixed << setprecision(3) << Cl2y << setw(3) << " " \
    << fixed << setprecision(3) << Cl2z << "\n"; 

    file.close();
    
    cout << "program selesai";
    
    return 0;

}

%%writefile ep-MgCl2.cpp

#include "iostream"
#include "cmath"
#include "iomanip"
#include "fstream"

using namespace std;

float densitas;          /* densitas ini memiliki satuan g cm^-3
                            dan ini harus diinput dari user */
// converter dari nilai cm^-3 ke A^3
float cm3_to_A3 = 1.0E-24;  
float numb_mol;
float Na = 6.022E23;            // bilangan avogadro    
float Mr_MgCl2 = 95.0;            // massa relatif dari MgCl2  
int numb_lat;                   // panjang sel simulasi    
int N;                          // jumlah molekul 
float volum;                    // volume dari molekul MgCl2 
float lx, ly, lz;               // panjang sel simulasi terhitung
float lat;                      // jarak antar molekul MgCl2 

// variabel perhitungan energi potensial 
float sigma = 2.725;            // satuan (A)
float epsilon = 3.725;          // satuan (kj/mol)
float rij, rcut, rcut2;
float Ep, Ep_LJ, A12, B6;
float rij6, rij12;
float dx, dy, dz;       


int main(){
    
    rcut = 3 * sigma;
    rcut2 = pow(rcut,2);

    cout << "program menghitung energi potensial dari MgCl2" << endl;
    cout << "Masukkan nilai densitas dari MgCl2 (eg. 1.0): ";
    cin >> densitas;

    numb_mol = densitas * (Na / Mr_MgCl2) * cm3_to_A3;
    cout << "Masukkan panjang sel simulasi (eg.5):";
    cin >> numb_lat;
    N = (pow(numb_lat,3)) * 3;
    volum = (float)N / numb_mol;

    // panjang sel simulasi secara perhitungan
    lx = pow(volum,(1.0/3.0));
    ly = lx;
    lz = lx;
    if(lx < (2*rcut)){
        cout << "sorry pak, sel simulasi kamu terlalu besar";
        exit(0);
    } else{
        cout << "panjang sel simulasi baru: " << lx << endl;
    }
    lat = lx / (float)numb_lat;

     // deklarasi variabel
    float pi = 3.14;
    
    // deklarasi molekul MgCl2 
    float Mgx[N], Mgy[N], Mgz[N];
    float Cl1x[N], Cl1y[N], Cl1z[N];
    float Cl2x[N], Cl2y[N], Cl2z[N];

    /* Cl    Cl
        \  /
         Mg
      molekul MgCl2 itu bentuk rigid
    */

    float rB_Mgx = 0.0;
    float rB_Mgy = 0.0;
    float rB_Mgz = 0.0;

    float rB_Cl1x = sin((104.52/2.0)*pi/180.0) * 1.0;
    float rB_Cl1y = 0.0;
    float rB_Cl1z = cos((104.52/2.0)*pi/180.0) * 1.0;

    float rB_Cl2x = - rB_Cl1x;
    float rB_Cl2y = rB_Cl1y;
    float rB_Cl2z = rB_Cl1z;

    int counter = 0;
    for(int i = 0; i < numb_lat; i++){ //x
        for(int j = 0; j < numb_lat; j++){ //y
            for(int k = 0; k < numb_lat; k++){ //z
                
                Mgx[counter] = rB_Mgx + (i*lat);
                Mgy[counter] = rB_Mgy + (j*lat);
                Mgz[counter] = rB_Mgz + (k*lat);

                Cl1x[counter] = Mgx[counter] + rB_Cl1x;
                Cl1y[counter] = Mgy[counter] + rB_Cl1y;
                Cl1z[counter] = Mgz[counter] + rB_Cl1z;

                Cl2x[counter] = Mgx[counter] + rB_Cl2x;
                Cl2y[counter] = Mgy[counter] + rB_Cl2y;Ep = 0.0;
    for(int a = 0; a < (N-1); a++){
        for(int b = (0+1); b < N; b++){
            dx = Mgx[a] - Mgx[b];
            dy = Mgy[a] - Mgy[b];
            dz = Mgz[a] - Mgz[b];

            dx = dx - round(dx/lx) * lx;
            dy = dy - round(dy/ly) * ly;
            dz = dz - round(dz/lz) * lz;

            rij = pow(dx,2) + pow(dy,2) + pow(dz,2);
            if(rij < rcut2){
                rij6 = pow(rij,3);
                rij12 = pow(rij6,2);
                Ep = (A12/rij12) - (B6/rij6);
                Ep_LJ = (A12/rij12) - (B6/rij6);
                Ep += Ep_LJ;
            }

        }
    }

    cout << "Energi potensial dari molekul MgCl2 ruah: " \
         << Ep/(float)N << "kj/mol";

    // hasil iterasi dimasukkan ke dalam file xyz
    ofstream file;
    file.open("MgCl2-hitung_ep.xyz");
    file << N << "\n" << endl;

    int Nw = N/3;

    // iterasi dalam output array 
    for(int m = 0; m < Nw; m++){
        file << setw(3) << "Mg" << setw(3) << " " \
             << fixed << setprecision(3) << Mgx[m] << setw(3) << " " \
             << fixed << setprecision(3) << Mgy[m] << setw(3) << " " \
             << fixed << setprecision(3) << Mgz[m] << "\n";

        file << setw(3) << "Cl" << setw(3) << " " \
             << fixed << setprecision(3) << Cl1x[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl1y[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl1z[m] << "\n";

        file << setw(3) << "Cl" << setw(3) << " " \
             << fixed << setprecision(3) << Cl2x[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl2y[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl2z[m] << "\n";
    }

    file.close();

    return 0;
}
                Cl2z[counter] = Mgz[counter] + rB_Cl2z;

                counter += 1;
            }
        }
    }
    // main hitung energi potensial
    A12 = 4.0 * epsilon * pow(sigma,12);
    B6 = 4.0 * epsilon * pow(sigma, 6);

    Ep = 0.0;
    for(int a = 0; a < (N-1); a++){
        for(int b = (0+1); b < N; b++){
            dx = Mgx[a] - Mgx[b];
            dy = Mgy[a] - Mgy[b];
            dz = Mgz[a] - Mgz[b];

            dx = dx - round(dx/lx) * lx;
            dy = dy - round(dy/ly) * ly;
            dz = dz - round(dz/lz) * lz;

            rij = pow(dx,2) + pow(dy,2) + pow(dz,2);
            if(rij < rcut2){
                rij6 = pow(rij,3);
                rij12 = pow(rij6,2);
                Ep = (A12/rij12) - (B6/rij6);
                Ep_LJ = (A12/rij12) - (B6/rij6);
                Ep += Ep_LJ;
            }

        }
    }

    cout << "Energi potensial dari molekul MgCl2 ruah: " \
         << Ep/(float)N << "kj/mol";

    // hasil iterasi dimasukkan ke dalam file xyz
    ofstream file;
    file.open("MgCl2-hitung_ep.xyz");
    file << N << "\n" << endl;

    int Nw = N/3;

    // iterasi dalam output array 
    for(int m = 0; m < Nw; m++){
        file << setw(3) << "Mg" << setw(3) << " " \
             << fixed << setprecision(3) << Mgx[m] << setw(3) << " " \
             << fixed << setprecision(3) << Mgy[m] << setw(3) << " " \
             << fixed << setprecision(3) << Mgz[m] << "\n";

        file << setw(3) << "Cl" << setw(3) << " " \
             << fixed << setprecision(3) << Cl1x[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl1y[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl1z[m] << "\n";

        file << setw(3) << "Cl" << setw(3) << " " \
             << fixed << setprecision(3) << Cl2x[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl2y[m] << setw(3) << " " \
             << fixed << setprecision(3) << Cl2z[m] << "\n";
    }

    file.close();

    return 0;
}
    
    
