#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
private:
public:
    analysis();
    ///Поиск минимального значения во входном массиве
    double minValue(QVector<double> x);

    ///Поиск максимального значения во входном массиве
    double maxValue(QVector<double> x);

    ///Математическое ожидания(среднее)
    double mean(QVector<double> x);

    ///Дисперсия
    double dispersion(QVector<double> x);

    ///Стандартное отклонение
    double standartDeviation(QVector<double> x);

    ///Среднеквадратическое отклонение
    double sqrtDeviation(QVector<double> x);

    ///Среднеквадратическая ошибка
    double sqrtError(QVector<double> x);

    ///Ассиметрия
    double assymetry(QVector<double> x);

    ///Коэффициент ассиметрии
    double assymetryCoeff(QVector<double> x);

    ///Эксцесс
    double excess(QVector<double> x);

    ///Куртозис
    double curtosis(QVector<double> x);

    ///Проверка на стационарность
    bool isStationary(QVector<double> x);

    ///Автокорреляция
    QVector<double> autocorrelation(QVector<double> x);

    ///Ковариация
    QVector<double> covariance(QVector<double> firstProcess, QVector<double> secondProcess);

    ///Амплитуда Фурье
    QVector<double> fourierAmplitude(QVector<double> inputData);

    ///Спектр Фурье
    QVector<double> fourierSpectrum(QVector<double> inputData, double window);

    ///Расчёт частоты для оси х графика спектра Фурье
    QVector<double> calculateFrequency(double delta_t, int N);

    ///Расчёт импульсной реакции фильтра низких частот Поттера
    QVector<double> lowpassFilterPotter(double fc, int m);

};

#endif // ANALYSIS_H
