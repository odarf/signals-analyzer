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
};

#endif // PROCESSING_H
