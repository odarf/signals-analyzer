#ifndef MODELING_H
#define MODELING_H

#include <iostream>
#include <QMainWindow>


class modeling
{
public:
    modeling();

    ///Линейный возрастающий тренд (размер, коэффициент)
    QVector<double> linearIncrement(int ,double);

    ///Линейный убывающий тренд (размер, коэффицинт)
    QVector<double> linearDecrement(int, double);

    ///Экспоненциальный возрастающий процесс (размер, альфа, бета)
    QVector<double> exponentialIncrement(int, double, double);

    ///Экспоненциальный убываюищй процесс (размер, альфа, бета)
    QVector<double> exponentialDecrement(int, double, double);

    ///Встроенный рандом(?)
    float embedRandom();

    ///Мой рандом(?)
    float randomGenerator(float);

    ///Модель кардиограммы
    QVector<double> cardiogram();
};

#endif // MODELING_H
