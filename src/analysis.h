#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
private:
    QVector<float> min_max;
public:
    analysis();
    double get_statistics();
    static double minValue(QVector<double>);
    static double maxValue(QVector<double>);
    double averageValue(QVector<double>);
    double dispersion(QVector<double>);
};

#endif // ANALYSIS_H
