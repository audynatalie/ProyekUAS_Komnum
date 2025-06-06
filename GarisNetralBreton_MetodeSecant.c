#include <stdio.h>
#include <math.h>

// mendefinisikan konstanta untuk masalah ini
#define FC_PRIME 300.0    // kekuatan tekan beton dalam kg/cm^2
#define FY 4000.0         // kekuatan tarik baja dalam kg/cm^2
#define ES 2000000.0      // modulus elastisitas baja dalam kg/cm^2
#define BETA1 0.85        // faktor untuk blok tegangan persegi panjang ekuivalen
#define EPSILON_CU 0.003  // regangan beton maksimum

// dimensi kolom
#define B 30.0            // lebar kolom dalam cm
#define H 50.0            // tinggi kolom dalam cm

// detail tulangan
#define BAR_DIAMETER 2.5  // diameter satu batang tulangan dalam cm (untuk Phi 25)
#define NUM_BARS_TOP 3    // jumlah batang di zona kompresi
#define NUM_BARS_BOTTOM 3 // jumlah batang di zona tarik

// menghitung luas area dari satu batang tulangan
#define AREA_BAR (M_PI * (BAR_DIAMETER / 2.0) * (BAR_DIAMETER / 2.0))

// kedalaman efektif d dan d' (mengasumsikan 4 cm penutup untuk pusat batang)
#define D (H - 4.0 - (BAR_DIAMETER / 2.0))  // kedalaman efektif untuk baja tarik
#define D_PRIME (4.0 + (BAR_DIAMETER / 2.0)) // kedalaman untuk baja kompresi

// fungsi untuk menghitung tegangan dalam baja berdasarkan regangan
double get_steel_stress(double epsilon_s) {
    double epsilon_y = FY / ES;
    
    if (fabs(epsilon_s) >= epsilon_y) {
        return (epsilon_s > 0) ? FY : -FY;  // tegangan maksimum pada regangan lebih besar dari regangan alir
    } else {
        return ES * epsilon_s;  // tegangan sesuai dengan hukum Hooke untuk regangan kecil
    }
}

// mendefinisikan fungsi f(c) untuk mencari akar
double f(double c) {
    if (c <= 0 || c > H) {
        return NAN; // nilai c tidak valid
    }
    
    double a = BETA1 * c;
    
    // gaya kompresi beton (Cc)
    double Cc = 0.85 * FC_PRIME * B * a;
    
    // regangan pada baja kompresi (epsilon_s_prime)
    double epsilon_s_prime = EPSILON_CU * (c - D_PRIME) / c;
    
    // tegangan pada baja kompresi (f_s_prime)
    double fs_prime = get_steel_stress(epsilon_s_prime);
    
    // gaya baja kompresi (Cs)
    double Cs = NUM_BARS_TOP * AREA_BAR * fs_prime;
    
    // regangan pada baja tarik (epsilon_s)
    double epsilon_s = EPSILON_CU * (D - c) / c;
    
    // tegangan pada baja tarik (f_s)
    double fs = get_steel_stress(epsilon_s);
    
    // gaya baja tarik (Ts)
    double Ts = NUM_BARS_BOTTOM * AREA_BAR * fs;
    
    // jumlah gaya (Cc + Cs - Ts)
    return Cc + Cs - Ts;
}

// implementasi metode secant
double secantMethod(double x0, double x1, double tolerance, int maxIterations) {
    int iteration = 0;
    double x_next;
    double fx0, fx1;
    
    printf("\n--- metode secant untuk perhitungan garis netral ---\n");
    printf("Iterasi\tc_i (cm)\tf(c_i) (kg)\tError (cm)\n");
    printf("---------------------------------------------------\n");
    
    while (iteration < maxIterations) {
        fx0 = f(x0);
        fx1 = f(x1);
        
        // memeriksa apakah nilai fungsi tidak valid
        if (isnan(fx0) || isnan(fx1)) {
            printf("error: nilai fungsi tidak valid ditemukan.\n");
            return NAN;
        }
        
        // menghindari pembagian dengan nol
        if (fabs(fx1 - fx0) < 1e-12) {
            printf("error: pembagian dengan nol. cobalah tebakan awal yang berbeda.\n");
            return NAN;
        }
        
        x_next = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        double current_error = fabs(x_next - x1);
        
        // mencetak hasil iterasi dengan format yang rapi
        printf("%-10d%-15.4f%-20.4f%-15.6f\n", 
               iteration + 1, x_next, f(x_next), current_error);
        
        if (current_error < tolerance) {
            printf("\nkonvergen setelah %d iterasi.\n", iteration + 1);
            return x_next;
        }
        
        x0 = x1;
        x1 = x_next;
        iteration++;
    }
    
    printf("\njumlah iterasi maksimum tercapai. tidak ada konvergensi.\n");
    return x_next;
}

int main() {
    // tebakan awal untuk 'c' (kedalaman garis netral)
    double x0 = 1.0;         // tebakan awal 1 (cm)
    double x1 = 10.0;        // tebakan awal 2 (cm)
    double tolerance = 0.0001; // toleransi untuk konvergensi (cm)
    int maxIterations = 100;   // jumlah iterasi maksimum
    
    printf("kolom beton bertulang - analisis garis netral\n");
    printf("dimensi: %.0f x %.0f cm\n", B, H);
    printf("tulangan: %d phi %.0f (atas) + %d phi %.0f (bawah)\n", 
           NUM_BARS_TOP, BAR_DIAMETER*10, NUM_BARS_BOTTOM, BAR_DIAMETER*10);
    printf("f'c = %.0f kg/cm², fy = %.0f kg/cm²\n", FC_PRIME, FY);
    
    double root = secantMethod(x0, x1, tolerance, maxIterations);
    
    if (!isnan(root)) {
        printf("\nhasil akhir:\n");
        printf("kedalaman garis netral (c) = %.4f cm\n", root);
        printf("nilai fungsi f(c) = %.6f kg\n", f(root));
        
        // perhitungan tambahan
        double a = BETA1 * root;
        printf("kedalaman blok tegangan ekuivalen (a) = %.4f cm\n", a);
    }
    
    return 0;
}

