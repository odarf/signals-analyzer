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

    ///Смещение на заданный коэффициент
    QVector<double> offset(QVector<double> data, double coeff);

    ///Неправдоподобные выбросы прменяются к данным
    QVector<double> spikes(QVector<double> data);

    ///Неправдоподобные выбросы
    float shift(float n);

    ///"Прицел"
    QVector<double> aim(QVector<double> data);

    ///Устранение смещения
    QVector<double> antiShift(QVector<double> data);

    ///Устранение неправдободобных выбросов
    QVector<double> antiSpike(QVector<double> spikedData);

    ///Выделение трендовой составляющей
    QVector<double> pickoutTrend(QVector<double> inputData);

    ///Добавление случайного шума к данным
    QVector<double> trendAddRandom(QVector<double> inputData);

    ///Выделение случайного шума из данных
    QVector<double> antiTrend(QVector<double> inputData);

    ///Отображение стандартного отклонения в замисимости от количества накоплений
    QVector<double> pdfTaskEight();

    ///Отношение стандартного отклонения(1) к стандартонму отклонению(i)
    QVector<double> pdfTaskEight2(QVector<double> inputStandartDev);

    ///Различное количество накоплений
    QVector<double> pdfTaskNine(QVector<double> inputData, int count);
};

#endif // PROCESSING_H
