#include "analysis.h"
#include "mainwindow.h"

#include <algorithm>
#include <math.h>

analysis::analysis(){}

QVector<float> analysis::minMax(QVector<float> x){
    QVector<float> min_max(2);
    min_max[0] = *std::min_element(x.constBegin(), x.constEnd());
    min_max[1] = *std::max_element(x.constBegin(), x.constEnd());
    return min_max;
}

float analysis::averageValue(QVector<double> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++) { sum+=x[i]; }
    return sum/x.length();
}

float analysis::dispersion(QVector<double> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++)
    { sum+=pow(x[i] - averageValue(x), 2); }
    return sum/x.length();
}
