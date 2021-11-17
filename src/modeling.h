#ifndef MODELING_H
#define MODELING_H

#include <iostream>
#include <QMainWindow>


class modeling
{
public:
    modeling();

    ///Линейный возрастающий тренд (размер, коэффициент)
    QVector<double> linearIncrement(int length, double k);

    ///Линейный убывающий тренд (размер, коэффицинт)
    QVector<double> linearDecrement(int length, double k);

    ///Экспоненциальный возрастающий процесс (размер, альфа, бета)
    QVector<double> exponentialIncrement(int length, double alpha=-0.01, double beta=11);

    ///Экспоненциальный убываюищй процесс (размер, альфа, бета)
    QVector<double> exponentialDecrement(int length, double alpha=-0.01, double beta=11);

    ///Встроенный рандом(?)
    float embedRandom();

    ///Мой рандом(?)
    float randomGenerator(float coefficient);

    ///Гармонический процесс (размер, амплитуда, частота)
    QVector<double> harmonic(int length, int amplitude, int frequency);

    ///Полигармонический процесс (размер, амплитуда, частота)
    QVector<double> polyharmonic(int length, int amplitude[], int frequency[]);

    ///Модель кардиограммы
    QVector<double> cardiogram();
};

#endif // MODELING_H
