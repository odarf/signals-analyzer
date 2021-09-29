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

QVector<float> averageValue(QVector<float> a){

    return 1.0;
}

float get_statistics(){
    QVector<float> x(1000);
    float mean, dispersion, deviation, sk, epsilon, asympt, excess = 0.0;

    return 1.0;
}
