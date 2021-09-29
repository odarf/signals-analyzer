#include "analysis.h"
#include "mainwindow.h"

#include <algorithm>

analysis::analysis()
{

}

QVector<float> minMax(QVector<float> x){
    QVector<float> min_max(2);
    min_max[0] = *std::min_element(x.constBegin(), x.constEnd());
    min_max[1] = *std::max_element(x.constBegin(), x.constEnd());
    return min_max;
}

float averageValue(QVector<float> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++) { sum+=x[i]; }
    return sum/x.length();
}

float dispersion(QVector<float> x){
    float sum = 0.0;
    for (int i=0; i<x.length()-1; i++)
    { sum+=(x[i] - averageValue(x)) * (x[i] - averageValue(x)); }
    return sum/x.length();
}
