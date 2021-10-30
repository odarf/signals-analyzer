#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
private:
public:
    analysis();
    double minValue(QVector<double>);
    double maxValue(QVector<double>);
    double mean(QVector<double>);
    double dispersion(QVector<double>);
    double standartDeviation(QVector<double>);
    double sqrtDeviation(QVector<double>);
    double sqrtError(QVector<double>);
    double assymetry(QVector<double>);
    double assymetryCoeff(QVector<double>);
    double excess(QVector<double>);
    double curtosis(QVector<double>);
    bool stationarity(QVector<double>);
    QVector<double> autocovariance(QVector<double>);
    QVector<double> covariance(QVector<double>, QVector<double>);
    QVector<double> density(QVector<double>);
};

#endif // ANALYSIS_H
