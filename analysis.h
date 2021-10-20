#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
private:
    QVector<float> min_max;
public:
    analysis();
    float get_statistics();
    QVector<float> minMax(QVector<float>);
    float averageValue(QVector<double>);
    float dispersion(QVector<double>);
};

#endif // ANALYSIS_H
