#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "analysis.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QtCharts>

#include <math.h>
#include <time.h>

#include <random>

using namespace QtCharts;
using namespace std;

const int LENGTH = 1000;
const double DELTA = 50.0;

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

float MainWindow::embedRandom(){
    float value = 0;
    int max = 1;
    int min = -1;
    int c = 1;
    value = rand() * (c * max - c * min) + c * min;
    return value;
}

// ---------------- Мой генератор случайных чисел ----------------
float MainWindow::randomGenerator(float coeff){
    float left = -1.0 * coeff;
    float right = 1.0 * coeff;
    default_random_engine engine{random_device()()};
    uniform_real_distribution<float> distribution{left, right};
    return distribution(engine);
}

// ---------------- Генератор спайков ----------------
float MainWindow::shift(float n){
    return randomGenerator(n) * 1000.0;
}

void MainWindow::on_pushButton_clicked(){
    //int selectedTask =  ui->cbSelectTask->currentText().toInt();
    processing processing;
    analysis analysis;

    cout << embedRandom() << endl;

    randomCoeff = (float)ui->sbCoeff->value();

    QVector<double> yLinInc(LENGTH), yLinDec(LENGTH),
                    yExpInc(LENGTH), yExpDec(LENGTH),
                    ySinMod(LENGTH), yMyRandom(LENGTH),
                    yTemp(LENGTH),   yEmbedRandom(LENGTH),
                    ySin(LENGTH),    yCombinated,
                    xLarge(3*LENGTH),
                    x(LENGTH),       ySmoothed(LENGTH),
                    xSin(LENGTH),    xTemp(LENGTH);

    double a = -1;
    double beta = 1;
    double alpha = -0.01;
    int amplitude[] = {10, 100, 15};
    int frequency[] = {4, 37, 173};
    double delta_t = 0.001;

    for(int X=0; X<LENGTH-1; X++){
        x[X] = X;

        yLinInc[X] = -a * X;
        yLinDec[X] = a * X + LENGTH;
        yExpInc[X] = (double)(beta * exp(-alpha * X));
        yExpDec[X] = (double)(beta * exp(alpha * X));
        yMyRandom[X] = randomGenerator(randomCoeff);
        yEmbedRandom[X] = embedRandom();
        ySin[X] = amplitude[0] * sin(2 * 3.14 * frequency[0] * X * delta_t);
    }

    for(int i = 0; i<3000; i++){
        xLarge[i] = i;
    }
    yCombinated.append(yLinInc);
    yCombinated.append(yLinDec);
    yCombinated.append(yExpInc);

    for(int count=0; count<3; count++){
        for(int i=0; i<5; i++){
            for(int X=0; X<LENGTH-1; X++){
                xSin[X] = X;
                ySinMod[X] += amplitude[0] * sin(2 * 3.14 * frequency[0] * X * delta_t) + DELTA; //график со смещением с одним значением амплитуды и частоты
            }
        }
    }

    for(int i=0; i<5; i++){
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> uni(0,1000);
        int randomValue = uni(rng);
        ySinMod[randomValue] = shift(2); //спайки
    }

    int L = 10;
    double sum = 0.0;
    int k = 0;
    for(int i = 0; i<LENGTH-1; i++){
        xTemp[i] = i;
        yTemp[i] = yLinDec[i];
        yTemp[i] += randomGenerator(randomCoeff);
    }
    for(int m = 0; m<LENGTH-L; m++){
        for(k = m; k<m+L-1; k++){
            sum += yTemp[k];
        }
        ySmoothed[m] = sum / L;
    }

// ---------------- Рисование ----------------
    ui->graphLinInc->addGraph();
    ui->graphLinInc->graph(0)->setData(x, yLinInc);
    ui->graphLinInc->xAxis->setRange(0, LENGTH+1);
    ui->graphLinInc->yAxis->setRange(analysis.minValue(yLinInc), analysis.maxValue(yLinInc)+1);
    ui->graphLinInc->replot();

    ui->graphLinDec->addGraph();
    ui->graphLinDec->graph(0)->setData(x, yLinDec);
    ui->graphLinDec->xAxis->setRange(0, LENGTH+1);
    ui->graphLinDec->yAxis->setRange(analysis.minValue(yLinDec), analysis.maxValue(yLinDec)+1);
    ui->graphLinDec->replot();

    ui->graphExpInc->addGraph();
    ui->graphExpInc->graph(0)->setData(x, yExpInc);
    ui->graphExpInc->xAxis->setRange(0, LENGTH+1);
    ui->graphExpInc->yAxis->setRange(analysis.minValue(yExpInc), analysis.maxValue(yExpInc)+1);
    ui->graphExpInc->replot();

    ui->graphExpDec->addGraph();
    ui->graphExpDec->graph(0)->setData(x, yExpDec);
    ui->graphExpDec->xAxis->setRange(0, LENGTH+1);
    ui->graphExpDec->yAxis->setRange(analysis.minValue(yExpDec), analysis.maxValue(yExpDec));
    ui->graphExpDec->replot();

    ui->graphCombinated->addGraph();
    ui->graphCombinated->graph(0)->setData(xLarge, yCombinated);
    ui->graphCombinated->xAxis->setRange(0, 3*LENGTH);
    ui->graphCombinated->yAxis->setRange(0, 3000);
    ui->graphCombinated->replot();

    ui->graphMyRandom->addGraph();
    ui->graphMyRandom->graph(0)->setData(x, yMyRandom);
    ui->graphMyRandom->xAxis->setRange(0, LENGTH+1);
    ui->graphMyRandom->yAxis->setRange(-randomCoeff-0.5, randomCoeff+0.5);
    ui->graphMyRandom->replot();

    ui->graphEmbedRandom->addGraph();
    ui->graphEmbedRandom->graph(0)->setData(x, yEmbedRandom);
    ui->graphEmbedRandom->xAxis->setRange(0, LENGTH+1);
    ui->graphEmbedRandom->yAxis->setRange(analysis.minValue(yEmbedRandom), analysis.maxValue(yEmbedRandom));
    ui->graphEmbedRandom->replot();

    ui->widget_6->addGraph();
    ui->widget_6->graph(0)->setData(xSin, ySinMod);
    ui->widget_6->xAxis->setRange(0, LENGTH);
    ui->widget_6->yAxis->setRange(-2500, 2500);
    ui->widget_6->replot();
    QThread::sleep(1);

    ui->graphHarmonicSin->addGraph();
    ui->graphHarmonicSin->graph(0)->setData(xSin, ySin);
    ui->graphHarmonicSin->xAxis->setRange(0, LENGTH);
    ui->graphHarmonicSin->yAxis->setRange(analysis.minValue(ySin) - 1, analysis.maxValue(ySin) + 1);
    ui->graphHarmonicSin->replot();
//-------------------- Статистики для гармонического графика ---------------------------
    ui->labelMinimumValue->setText(QString::number(analysis.minValue(ySin)));
    ui->labelMaximumValue->setText(QString::number(analysis.maxValue(ySin)));
    ui->labelMeanValue->setText(QString::number(analysis.mean(ySin)));
    ui->labelDispersionValue->setText(QString::number(analysis.dispersion(ySin)));
    ui->labelStdDevValue->setText(QString::number(analysis.standartDeviation(ySin)));
    ui->labelSqrtDevValue->setText(QString::number(analysis.sqrtDeviation(ySin)));
    ui->labelAssymetryValue->setText(QString::number(analysis.assymetry(ySin)));
    ui->labelAssymetryCoeffValue->setText(QString::number(analysis.assymetryCoeff(ySin)));
    ui->labelExcessValue->setText(QString::number(analysis.excess(ySin)));
    ui->labelCurtosisValue->setText(QString::number(analysis.curtosis(ySin)));
    ui->labelIsStationarity->setText(analysis.isStationary(yEmbedRandom) ? "Стационарен" : "Не стационарен");

    ui->widget_7->addGraph();
    ui->widget_7->graph(0)->setData(xTemp, ySmoothed);
    ui->widget_7->xAxis->setRange(0, LENGTH);
    ui->widget_7->yAxis->setRange(0, 1000);
    ui->widget_7->replot();

    ui->graphLinInc->clearGraphs();
    ui->graphLinDec->clearGraphs();
    ui->graphExpInc->clearGraphs();
    ui->graphExpDec->clearGraphs();
    ui->graphMyRandom->clearGraphs();
    ui->widget_6->clearGraphs();
    ui->widget_7->clearGraphs();
}
