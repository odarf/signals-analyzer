#include "analysis.h"
#include "mainwindow.h"

#include <algorithm>
#include <math.h>
#include <QDebug>

analysis::analysis(){}

double analysis::minValue(QVector<double> x){
    double min = *std::min_element(x.constBegin(), x.constEnd());
    return min;
}

double analysis::maxValue(QVector<double> x){
    double max = *std::max_element(x.constBegin(), x.constEnd());
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

bool analysis::isStationary(QVector<double> x){
    int m = 1000;
    int e = 10;
    int n = x.length();
    QVector<double> mArr;
    QVector<double> dArr;
    double mean = 0;
    double dispersion = 0;
    double avg = 0;

    for(int i = 0; i<m; i++){
        for(int k = i * (n/m); k < i*n/m+n/m; k++){
            mean += x[k];
        }
        avg = mean / n;
        mArr.append(avg);

        for(int k = i*m; k<i*n/m+n/m; k++){
            dispersion += pow((x[k] - avg), 2);
        }
        dArr.append(dispersion / n);

        for(int j = 0; j<i-1; j++){
            if(i != 0){
                if(std::max(mArr[i], mArr[j]) - std::min(mArr[i], mArr[j]) / std::max(mArr[i], mArr[j]) * 100 > e){
                    return false;
                }
                if(std::max(dArr[i], dArr[j]) - std::min(dArr[i], dArr[j]) / std::max(dArr[i], dArr[j]) * 100 > e){
                    return false;
                }
                //return (std::max(mArr[i], mArr[j]) - std::min(mArr[i], mArr[j]) / std::max(mArr[i], mArr[j]) * 100 > e) ? false : NULL;
                //return (std::max(dArr[i], dArr[j]) - std::min(dArr[i], dArr[j]) / std::max(dArr[i], dArr[j]) * 100 > e) ? false : NULL;
            }
        }
    }
    return true;
}

QVector<double> analysis::autocorrelation(QVector<double> x){
    double avg = mean(x);
    int n = x.length();
    QVector<double> autocovariance;

    for(int i = 0; i<n; i++){
        double sum1 = 0;
        double sum2 = 0;
        for(int j = 0; j<n-i; j++){
            sum1 += (x[j] - avg) * (x[j+i] - avg);
        }
        for(int k = 0; k<n; k++){
            sum2 += pow((x[k] - avg), 2);
        }
        autocovariance.append(sum1/sum2);
    }
    return autocovariance;
}

QVector<double> analysis::covariance(QVector<double> firstProcess, QVector<double> secondProcess){
    double avgX = mean(firstProcess);
    double avgY = mean(secondProcess);
    int n = firstProcess.length();
    QVector<double> covariance;

    for(int i = 0; i<n; i++){
        double sum = 0;
        for(int j = 0; j<n-i; j++){
            sum += (firstProcess[j] - avgX) * (secondProcess[j + i] - avgY);
        }
        covariance.append(sum/n);
    }
    return covariance;
}

QVector<double> analysis::fourierAmplitude(QVector<double> inputData){
    int length = inputData.length();
    QVector<double> outputData;

    for(int i = 0; i<length/2; i++){
        double real = 0;
        double imagine = 0;

        for(int j = 0; j<length; j++){
            real += inputData[j] * std::cos((2 * 3.14 * i * j) / length);
            imagine += inputData[j] * std::sin((2 * 3.14 * i * j) / length);
        }
        real /= length;
        imagine /= length;
        outputData.append(std::sqrt(std::pow(real, 2) + std::pow(imagine, 2)));
    }
    return outputData;
}

QVector<double> analysis::fourierSpectrum(QVector<double> inputData, double window){
    int length = inputData.length();
    double w = (length - length * window) / 2;
    QVector<double> outputData;

    for(int i = 0; i<w; i++){
        inputData[i] = 0;
        inputData[length - i - 1] = 0;
    }
    for(int i = 0; i<length/2; i++){
        double real = 0;
        double imagine = 0;
        for(int j = 0; j<length; j++){
            real += inputData[j] * std::cos((2 * 3.14 * i * j) / length);
            imagine += inputData[j] * std::sin((2 * 3.14 * i * j) / length);
        }
        real /= length;
        imagine /= length;
        outputData.append(std::sqrt(std::pow(real, 2) + std::pow(imagine, 2)));
    }
    return outputData;
}
