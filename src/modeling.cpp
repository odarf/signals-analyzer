#include "modeling.h"
#include "qmath.h"

#include <iostream>
#include <QMainWindow>

#include <cstdlib>

#include <math.h>
#include <random>
#include <time.h>

using namespace std;

modeling::modeling() {}

QVector<double> modeling::linearIncrement(int length, double k){
    QVector<double> linearInc(length);
    for(int x = 0; x<length; x++){
        linearInc[x] = -k * x;
    }
    return linearInc;
}

QVector<double> modeling::linearDecrement(int length, double k){
    QVector<double> linearDec(length);
    for(int x = 0; x<length; x++){
        linearDec[x] = k * x + length;
    }
    return linearDec;
}

QVector<double> modeling::exponentialIncrement(int length, double alpha, double beta){
    QVector<double> exponentialInc(length);
    for(int x = 0; x<length; x++){
        exponentialInc[x] = beta * exp(-alpha * x);
    }
    return exponentialInc;
}

QVector<double> modeling::exponentialDecrement(int length, double alpha, double beta){
    QVector<double> exponentialDec(length);
    for(int x = 0; x<length; x++){
        exponentialDec[x] = beta * exp(alpha * x);
    }
    return exponentialDec;
}

float modeling::embedRandom(){
    int const max = 1;
    int const min = -1;
    return static_cast<float>(rand() * (1.0 / (static_cast<float>(RAND_MAX) + 1.0)) * (max - min) + min);
}

float modeling::randomGenerator(float coefficient){
    float const left = -1.0 * coefficient;
    float const right = 1.0 * coefficient;
    default_random_engine engine{random_device()()};
    uniform_real_distribution<float> distribution{left, right};
    return distribution(engine);
}

QVector<double> modeling::harmonic(int length, int amplitude, int frequency){
    QVector<double> harmonic(length);
    double delta_t = 0.001;
    for(int x = 0; x<length; x++){
        harmonic[x] = amplitude * sin(2 * 3.14 * frequency * x * delta_t);
    }
    return harmonic;
}

QVector<double> modeling::polyharmonic(int length, int amplitude[], int frequency[]){
    QVector<double> polyharmonic(length);
    for(int count=0; count<3; count++){
        for(int i=0; i<5; i++){
            for(int x=0; x<length; x++){
                polyharmonic[x] += amplitude[count] * sin(2 * 3.14 * frequency[count] * x);// * delta_t);
            }
        }
    }
    return polyharmonic;
}

QVector<double> modeling::cardiogram(){
    int const N = 1000;
    int const M = 200;
    QVector<double> x(N), h(M), y(N+M);

    float alpha = 10;
    int frequency = 4;
    int const tempN = N + M - 1;

    for(int i = 0; i<h.length(); i++){
        h[i] = sin(2 * 3.14 * frequency * (i*0.005)) * exp(-alpha * (i*0.005));
    }

    x[200] = 120;
    x[400] = 130;
    x[600] = 110;

    for(auto i(0); i<tempN; ++i){
        int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
        int const jmx = (i < N - 1)  ? i           : N-1;
        for(auto j(jmn); j<=jmx; ++j){
            y[i] += x[j] * h[i-j];
        }
    }
    return y;
}
