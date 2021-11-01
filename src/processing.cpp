#include "processing.h"
#include "mainwindow.h"
#include <math.h>
#include <qmath.h>

processing::processing(){}


QVector<double> processing::offset(QVector<double> data, double coeff){
    for(int i = 0; i<data.length(); i++){ data[i] += coeff; }
    return data;
}

float processing::shift(float n){
    return MainWindow().randomGenerator(n) * 1000.0;
}

QVector<double> processing::spikes(QVector<double> data){
    int randomValue = 0;
    for(int i=0; i<5; i++){
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> uid(0,1000);
        randomValue = uid(rng);
        data[randomValue] = shift(2); //спайки
    }
    return data;
}

QVector<double> processing::aim(QVector<double> data){
    QVector<double> output;
    int coeff = 10;
    double temp = 0;
    for(int i = 0; i<data.length(); i++){
        temp = data[i] + coeff;
        output.append(temp);
    }
    return output;
}

QVector<double> processing::antiShift(QVector<double> data){
    QVector<double> aimData = aim(data);
    QVector<double> output;
    analysis analysis;
    float mean = analysis.mean(aimData);
    double temp = 0;
    for(int i = 0; i<data.length(); i++){
        temp = aimData[i] - mean;
        output.append(temp);
    }
    return output;
}

QVector<double> processing::antiSpike(QVector<double> spikedData){
    QThread::sleep(1);
    auto output = spikedData;
    for(int i = 1; i<output.length()-1; i++){
        if (qAbs(output[i] / output[i-1]) > 2){
            output[i] = (output[i-1] + output[i+1]) / 2.0;
        }
        else{ output[i] = output[i]; }
    }
    return output;
}

QVector<double> processing::pickoutTrend(QVector<double> inputData){
    double L = 10;
    double sum = 0;
    QVector<double> outputData;
    for(int i = 0; i<inputData.length() - L; i++){
        sum = 0;
        for(int j = i; j<i+L; j++){
            sum += inputData[j];
        }
        outputData.append(sum/L);
    }
    return outputData;
}

QVector<double> processing::trendAddRandom(QVector<double> inputData){
    MainWindow main;
    QVector<double> outputData;
    for(int i = 0; i<inputData.length(); i++){
        inputData[i] += main.randomGenerator(5);
        outputData.append(inputData[i]);
    }
    return outputData;
}

QVector<double> processing::antiTrend(QVector<double> inputData){
    QVector<double> pickTrend = pickoutTrend(inputData);
    QVector<double> outputData;
    for(int i = 0; i<pickTrend.length(); i++){
        outputData.append(inputData[i] - pickTrend[i]);
    }
    return outputData;
}

QVector<double> processing::pdfTaskEight(QVector<double> inputData){
    QVector<double> standartDev(inputData.length());
    QVector<double> temp;
    analysis analysis;
    for(int i = 0; i<10; i++){
        temp.append(trendAddRandom(inputData));
        for(int j = 0; j<standartDev.length(); j++){
            standartDev[j] += temp[j] * 0.1;
        }
        temp.clear();
    }
    return standartDev;
}
