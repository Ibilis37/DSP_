#include <vector>
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include<algorithm>
#include<bitset>
#define PI 3.14159265
using namespace std;
using namespace std::complex_literals;
typedef complex<double> Complex;
const std::complex<double> i(0.0, 1.0);


vector<Complex> read_file(string filename) {
    vector<Complex> signal;
    std::ifstream in(filename);
    if (in.is_open()) {
        double a, b;
        while (in >> a >> b)
            signal.push_back((a + b * i));
    }
    else
    {
	    cout << "File is not open" << endl;
    }
    in.close();
    return signal;
}

void write(vector<Complex> x, string filename) {
    ofstream out;
    out.open(filename);
    if (out.is_open()) {
        for (int i = 0; i < x.size(); i++)
            out << setprecision(20) << x[i].real() << " " << setprecision(20) << x[i].imag() << endl;
    }
    else { cout << "Error while writing" << endl; }
    out.close();
}

vector<Complex> DFT(vector<Complex> x, bool reverse) {
    int N = x.size();
    int flag = pow(-1, reverse+1);
    vector<Complex> y(N);
    for (int k = 0; k < N; k++) {
        y[k] = 0;
        for (int j = 0; j < N; j++) {
            y[k] += x[j] * exp(flag * 2 * PI / N * j * k * i);
        }
        y[k] = y[k] / sqrt(N);
    }
    return y;
}

void Permutation(vector<Complex>& x)
{
    int N = x.size();
    vector<bool> arr(N, false);
    for (int k = 0; k < N; k++)
    {
        string binary = std::bitset<32>(k).to_string(); // получчаем строку - двоичное представление числа
        binary.erase(0, -log2(N) + binary.size());
        reverse(binary.begin(), binary.end()); // инвертируем строку
        int k2 = std::bitset<32>(binary).to_ulong(); // получаем число, соответствующее инвертированной строке
        // переставляем элементы массива
        if (!(arr[k] || arr[k2]))
        {
            Complex a = x[k];
            x[k] = x[k2];
            x[k2] = a;
            arr[k] = true;

        }
    }
}

vector<Complex> FFT_(vector<Complex> x, bool reverse) {
    int N = x.size();
    Permutation(x);
    if (reverse) {
        for (int p = 0; p < N; p++)
            x[p] = conj(x[p]);
    }

    vector<Complex> y(N);
    int n = log2(N);
    for (int k = 1; k < n + 1; k++) {
        for (int j = 0; j < pow(2, n - k); j++) {
            for (int l = 0; l < (pow(2, k - 1)); l++) {
                Complex omega = exp(-2 * PI * l / pow(2, k) * i);
                y[j * pow(2, k) + l] = x[j * pow(2, k) + l] + omega * x[j * pow(2, k) + l + pow(2, k - 1)];
                y[j * pow(2, k) + l + pow(2, k - 1)] = x[j * pow(2, k) + l] - omega * x[j * pow(2, k) + l + pow(2, k - 1)];
            }
        }
    }

    for (int p = 0; p < N ; p++)
        y[p] = y[p] / sqrt(N);
    if (reverse) {
        for (int p = 0; p < N; p++)
            y[p] = conj(y[p]);
    }
    return y;
}
std::vector<Complex> FFT(std::vector<Complex> x, bool reverse) {
    int N = x.size();
    int n = log2(N);
    Complex h;
    Complex y1;
    int ind; //индексы для перестановок
    int new_ind;

    //перестановка
    for (int i = 1; i < N; ++i) {
        ind = i;
        new_ind = 0;
        for (int j = 0; j < n; ++j)
        {
            new_ind = (new_ind << 1) | (ind & 1);
            ind >>= 1;
        }
        if (new_ind > i) swap(x[i], x[new_ind]);
    }

    //сам бпф
    for (int k = 1; k < n + 1; ++k) {
        for (int j = 0; j < pow(2, n - k); ++j) {
            for (int l = 0; l < pow(2, k - 1); ++l) {
                h = Complex(cos(2 * PI * l / pow(2, k)), pow(-1, reverse+1) * sin(2 * PI * l / pow(2, k))) * x[j * pow(2, k) + l + pow(2, k - 1)]; //w_n^k*y1(k)
                y1 = x[j * pow(2, k) + l] + h;
                x[j * pow(2, k) + l + pow(2, k - 1)] = (x[j * pow(2, k) + l] - h) / sqrt(2);
                x[j * pow(2, k) + l] = y1 / sqrt(2);
            }
        }
    }
    return x;
}

vector<Complex> Convolution(vector<Complex> x, vector<Complex> y) {
    if (x.size() > y.size()) x.swap(y);
    int N = x.size() + y.size();
    vector<Complex> c(N);

    for (int n = 0; n < N; n++) {
        c[n] = 0;
        for (int k = 0; k < x.size(); k++) {
            if (n - k >= y.size() || k > n) {
                continue;
            }
            else { c[n] += x[k] * y[n - k]; }
        }
    }
    return c;
}

vector<Complex> FConvolution(vector<Complex> x, vector<Complex> y) {
    Complex zero = (0, 0);
    int M = x.size();
    int L = y.size();
    int N2 = M > L ? 2 * M : 2 * L;
    int N = pow(2, ceil(log2(N2)));
    vector<Complex> u(N);
    x.resize(N, zero);
    y.resize(N, zero);

    FFT(x, false);
    FFT(y, false);
    for (int i = 0; i < N; i++) {
        u[i] = sqrt(N) * x[i] * y[i];
    }
    FFT(u, true);
    return u;
}

double RMSE(vector<Complex> x, vector<Complex> y) {
    double err = 0;
    for (int i = 0; i < x.size(); i++)
        err += abs(pow(x[i] - y[i], 2));
    err = sqrt(err/x.size());
    return err;
}

void main() {

    // 3.1
    vector<Complex> signal = read_file("C:\\Users\\Диана\\Desktop\\data.txt");
    cout << "3.1 - RMSE of difference between odft(dft(signal)) and signal - " << RMSE(signal, DFT(DFT(signal, false), true)) << endl;
    cout <<"3.1 - RMSE of differente between offt(fft(signal)) and signal - " << RMSE(signal, FFT(FFT(signal, false), true)) << endl;
    // 3.2
    cout << "3.2 - RMSE of difference between fft and dft - " << RMSE(DFT(signal, false), FFT(signal, false)) << endl;
    // 3.3
    vector<Complex> python_fft = read_file("C:\\Users\\Диана\\Desktop\\fft_python.txt");
    cout <<"3.3 - RMSE of difference between fft and python_fft " << RMSE(python_fft, FFT(signal, false)) << endl;
    // 4
    int N = 10;
    vector<Complex> time_dft_fft(N);
    double time_dft;
	double time_fft;
    for (int j = 0; j < N; j++) {
        vector<Complex> dft(pow(2,j));
        vector<Complex> fft(pow(2,j));
        auto start = chrono::steady_clock::now();
        DFT(dft, false);
        auto end = chrono::steady_clock::now();
        time_dft =  chrono::duration<double, std::nano>(end - start).count();

        start = std::chrono::steady_clock::now();
        FFT(fft, false);
        end = std::chrono::steady_clock::now();
        time_fft =  std::chrono::duration<double, std::nano>(end - start).count();
        time_dft_fft[j] = Complex(time_dft, time_fft);
    }
    write(time_dft_fft, "C:\\Users\\Диана\\Desktop\\time_dft_fft.txt");

    //5-7

    vector<Complex> x = read_file("C:\\Users\\Диана\\Desktop\\x.txt");
    vector<Complex> y = read_file("C:\\Users\\Диана\\Desktop\\y.txt");
    vector<Complex> python_conv = read_file("C:\\Users\\Диана\\Desktop\\python_conv.txt");

    cout << "5 - RMSE of convolution - " << RMSE(python_conv, Convolution(x, y)) << endl;
    cout << "6 - RMSE of fconvolution - " << RMSE(python_conv, FConvolution(x, y)) << endl;
    
    // 8
	N = 10;
    vector<Complex> time_conv1(N);
    vector<Complex> time_fconv1(N);
    vector<Complex> time_conv2(N);
    vector<Complex> time_fconv2(N);
    vector<Complex> a;
    vector<Complex> b;
    vector<Complex> c;
    a.resize(256, 1);

    for (int i = 0; i < N; i++) {
        b.resize(pow(2, i), 1);

        auto start = chrono::steady_clock::now();
        c = Convolution(a, b);
        auto end = chrono::steady_clock::now();
        time_conv1[i] = Complex(i, chrono::duration<double, nano>(end - start).count());

        start = chrono::steady_clock::now();
        c = FConvolution(a, b);
        end = chrono::steady_clock::now();
        time_fconv1[i] = Complex(i, chrono::duration<double, nano>(end - start).count());
    }

    for (int i = 0; i < N; i++) {
        a.resize(pow(2, i), 1);
        b.resize(pow(2, i), 1);

        auto start = chrono::steady_clock::now();
        c = Convolution(a, b);
        auto end = chrono::steady_clock::now();
        time_conv2[i] = Complex(i, chrono::duration<double, nano>(end - start).count());

        start = chrono::steady_clock::now();
        c = FConvolution(a, b);
        end = chrono::steady_clock::now();
        time_fconv2[i] = Complex(i, chrono::duration<double, nano>(end - start).count());
    }

    write(time_conv1, "C:\\Users\\Диана\\Desktop\\time_conv1.txt");
    write(time_conv2, "C:\\Users\\Диана\\Desktop\\time_conv2.txt");
    write(time_fconv1, "C:\\Users\\Диана\\Desktop\\time_fconv1.txt");
    write(time_fconv2, "C:\\Users\\Диана\\Desktop\\time_fconv2.txt");
}