#ifndef PROCESSING_H
#define PROCESSING_H

#include "analysis.h"
#include "mainwindow.h"
#include <QMainWindow>
#include <iostream>
#include <math.h>
#include <time.h>
#include <random>

using namespace std;

class processing
{
public:
    processing();
    QVector<double> offset(QVector<double>, double);
    QVector<double> spikes(QVector<double>);
    float shift(float);
    QVector<double> aim(QVector<double>);
    QVector<double> antiShift(QVector<double>);
    QVector<double> antiSpike(QVector<double>);
    QVector<double> pickoutTrend(QVector<double>);
    QVector<double> trendAddRandom(QVector<double>);
    QVector<double> antiTrend(QVector<double>);
    QVector<double> pdfTaskEight();
    QVector<double> pdfTaskEight2(QVector<double>);
    QVector<double> pdfTaskNine(QVector<double>, int);
};

#endif // PROCESSING_H
