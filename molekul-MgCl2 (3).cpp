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
