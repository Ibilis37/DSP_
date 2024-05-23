# DSP_Lab2
Задание
### Задание 1. 
  > Реализовать на С или С++ алгоритмы непосредственного вычисления ДПФ и ОДПФ по формулам (1) и (2) для комплексного входного сигнала с двойной точностью (double). Входные данные загружать из текстового файла (разделитель – пробел), сгенерированного, например, в MATLAB.
```c++
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
```
### Задание 2. 
  >Реализовать на С или С++ алгоритмы прямого и обратного БПФ для
  комплексного входного сигнала длиной n 2 , n – любое натуральное число:
  а) с прореживанием по времени и двоично-инверсными перестановками
  (вариант 1);
```c++
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

std::vector<Complex> FFT(std::vector<Complex> x, bool reverse) {
    int N = x.size();
    int n = log2(N);
    Complex h;
    Complex y1;
    Permutation(x);
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

```
### Задание 3.
  > Убедиться в корректности работы алгоритмов:
  а) проверить выполнение равенства X  ОДПФ ДПФX, а также равенства
  X  ОБПФ БПФX;
  б) сравнить результаты ДПФ(Х) и БПФ(Х);
  в) сравнить результаты работы реализованного алгоритма, например, с
  результатами процедуры fft, встроенной в MATLAB.
  (рекомендуется для сравнения использовать значение ошибки)
```c++
double RMSE(vector<Complex> x, vector<Complex> y) {
    double err = 0;
    for (int i = 0; i < x.size(); i++)
        err += abs(pow(x[i] - y[i], 2));
    err = sqrt(err/x.size());
    return err;
}
int main{
    // 3.1
    vector<Complex> signal = read_file("C:\\Users\\Диана\\Desktop\\data.txt");
    cout << "3.1 - RMSE of difference between odft(dft(signal)) and signal - " << RMSE(signal, DFT(DFT(signal, false), true)) << endl;
    cout <<"3.1 - RMSE of differente between offt(fft(signal)) and signal - " << RMSE(signal, FFT(FFT(signal, false), true)) << endl;
    // 3.2
    cout << "3.2 - RMSE of difference between fft and dft - " << RMSE(DFT(signal, false), FFT(signal, false)) << endl;
    // 3.3
    vector<Complex> python_fft = read_file("C:\\Users\\Диана\\Desktop\\fft_python.txt");
    cout <<"3.3 - RMSE of difference between fft and python_fft " << RMSE(python_fft, FFT(signal, false)) << endl;
}
```
>  **3.1 - RMSE of difference between odft(dft(signal)) and signal - 2.46384e-08 \
	3.1 - RMSE of differente between offt(fft(signal)) and signal - 5.89683e-16 \
	3.2 - RMSE of difference between fft and dft - 3.68293e-08 \
	3.3 - RMSE of difference between fft and python_fft 1.21502e-09**
> 
### Задание 4.
  > Проанализировать зависимость времени выполнения БПФ и непосредственного
  вычисления ДПФ от длины N преобразования.
```c++
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
```
### Задание 5.
  > Реализовать на С или С++ процедуру прямого вычисления свертки двух
  последовательностей по формуле (3). Входные данные загружать из текстового
  файла (разделитель – пробел), сгенерированного, например, в MATLAB.
```c++
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
```
> **5 - RMSE of convolution - 1.5773e-15**
### Задание 6. 
  > Реализовать процедуру нахождения дискретной свертки, основанную на БПФ.
  При вычислении БПФ использовать результаты п. 2 задания.
```c++
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
```
> **6 - RMSE of fconvolution - 7.72303**
### Задание 7. 
  > Убедится в корректности работы процедуры из п. 5 и п. 6 задания, сравнив
  полученные результаты с результатами работы встроенной функций MATLAB
  conv. (рекомендуется для сравнения использовать значение ошибки)
```c++
    vector<Complex> x = read_file("C:\\Users\\Диана\\Desktop\\x.txt");
    vector<Complex> y = read_file("C:\\Users\\Диана\\Desktop\\y.txt");
    vector<Complex> python_conv = read_file("C:\\Users\\Диана\\Desktop\\python_conv.txt");

    cout << "5 - RMSE of convolution - " << RMSE(python_conv, Convolution(x, y)) << endl;
    cout << "6 - RMSE of fconvolution - " << RMSE(python_conv, FConvolution(x, y)) << endl;
```
### Задание 8.
  > Сравнить производительность алгоритмов вычисления свертки по
  определению (3) и с помощью БПФ в двух случаях: когда размер одной из
  последовательностей фиксирован, и когда меняются длины обеих
  последовательностей.
```c++
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
```
