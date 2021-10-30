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

double analysis::averageValue(QVector<double> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++) { sum+=x[i]; }
    return sum/x.length();
}

double analysis::dispersion(QVector<double> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++)
    { sum+=pow(x[i] - averageValue(x), 2); }
    return sum/x.length();
}
