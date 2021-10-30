#include "analysis.h"
#include "mainwindow.h"

#include <algorithm>
#include <math.h>

analysis::analysis(){}

double analysis::minValue(QVector<double> x){
    float min = *std::min_element(x.constBegin(), x.constEnd());
    return min;
}

double analysis::maxValue(QVector<double> x){
    float max = *std::max_element(x.constBegin(), x.constEnd());
    return max;
}

double analysis::mean(QVector<double> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++) { sum+=x[i]; }
    return sum/x.length();
}

double analysis::dispersion(QVector<double> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++)
    { sum+=pow(x[i] - mean(x), 2); }
    return sum/x.length();
}

double analysis::standartDeviation(QVector<double> x){
    return pow(dispersion(x), 2);
}

double analysis::sqrtDeviation(QVector<double> x){
    double sqrtDeviation = 0;
    for (int i = 0; i<x.length(); i++) { sqrtDeviation += pow(x[i], 2); }
    return sqrtDeviation/x.length();
}

double analysis::sqrtError(QVector<double> x){
    return pow(sqrtDeviation(x), 2);
}

double analysis::assymetry(QVector<double> x){
    double assymetry = 0;
    for (int i = 0; i <x.length(); i++) { assymetry += pow(x[i] - mean(x), 3); }
    return assymetry/x.length();
}

double analysis::assymetryCoeff(QVector<double> x){
    return assymetry(x) / pow(standartDeviation(x) ,3);
}

double analysis::excess(QVector<double> x){
    double excess = 0;
    for (int i = 0; i<x.length()-1; i++) { excess += pow(x[i] - mean(x), 4); }
    return excess/x.length();
}

double analysis::curtosis(QVector<double> x){
    return excess(x) / (pow(standartDeviation(x), 4) - 3);
}

bool analysis::stationarity(QVector<double> x){
    return true; //дописать!
}

QVector<double> analysis::autocovariance(QVector<double> x){
    double avg = mean(x);
    double n = x.length();
    QVector<double> autocovariance;

    for(int i = 0; i<n; i++){
        double sum1 = 0;
        double sum2 = 0;
        for(int j = 0; j<n; j++){
            sum1 += (x[j] - avg) * (x[j+1] - avg);
        }
        for(int k = 0; k<n; k++){
            sum2 += pow((x[k] - avg), 2);
        }
        autocovariance.append(sum1/sum2);
    }
    return autocovariance;
}

QVector<double> analysis::covariance(QVector<double> x, QVector<double> y){
    double avgX = mean(x);
    double avgY = mean(y);
    int n = x.length();
    QVector<double> covariance;

    for(int i = 0; i<n; i++){
        double sum = 0;
        for(int j = 0; j<n-i; j++){
            sum += (x[j] - avgX) * (y[j + i] - avgY);
        }
        covariance.append(sum/n);
    }
    return covariance;
}

/*QVector<double> analysis::density(QVector<double> x){
    //дописать!

}*/

