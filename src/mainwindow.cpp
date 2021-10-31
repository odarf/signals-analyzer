#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "analysis.h"
#include "qmath.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QtCharts>
#include <QtWidgets>
#include <QtWidgets/QHBoxLayout>
#include <QHBoxLayout>

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

void MainWindow::on_pushButton_clicked(){
    //int selectedTask =  ui->cbSelectTask->currentText().toInt();
    processing processing;
    analysis analysis;

    randomCoeff = (float)ui->sbCoeff->value();

    QVector<double> yLinInc(LENGTH), yLinDec(LENGTH),
                    yExpInc(LENGTH), yExpDec(LENGTH),
                    ySinMod(LENGTH), yMyRandom(LENGTH),
                    yTemp(LENGTH),   yEmbedRandom1(LENGTH),
                    ySin(LENGTH),    yCombinated,
                    xLarge(3*LENGTH),yEmbedRandom2(LENGTH),
                    x(LENGTH),       ySmoothed(LENGTH),
                    xSin(LENGTH),    xTemp(LENGTH);

    double a = -1;
    double beta = 11;
    double alpha = -0.01;
    int amplitude[] = {10, 100, 15};
    int frequency[] = {4, 37, 173};
    double delta_t = 0.001;

    for(int X=0; X<LENGTH; X++){
        x[X] = X;

        yLinInc[X] = -a * X;
        yLinDec[X] = a * X + LENGTH;
        yExpInc[X] = (double)(beta * exp(-alpha * X));
        yExpDec[X] = (double)(beta * exp(alpha * X));
        yMyRandom[X] = randomGenerator(randomCoeff);
        yEmbedRandom1[X] = embedRandom();
        yEmbedRandom2[X] = embedRandom();
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
            for(int X=0; X<LENGTH; X++){
                ySinMod[X] += amplitude[count] * sin(2 * 3.14 * frequency[count] * X * delta_t) + DELTA; //график со смещением с одним значением амплитуды и частоты
            }
        }
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

// ---------------- Task 1 ----------------
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
// ----------------------------------------

// ---------------- Task 2 ----------------
    ui->graphMyRandom->addGraph();
    ui->graphMyRandom->graph(0)->setData(x, yMyRandom);
    ui->graphMyRandom->xAxis->setRange(0, LENGTH+1);
    ui->graphMyRandom->yAxis->setRange(-randomCoeff-0.5, randomCoeff+0.5);
    ui->graphMyRandom->replot();

    ui->graphEmbedRandom->addGraph();
    ui->graphEmbedRandom->graph(0)->setData(x, yEmbedRandom1);
    ui->graphEmbedRandom->xAxis->setRange(0, LENGTH+1);
    ui->graphEmbedRandom->yAxis->setRange(analysis.minValue(yEmbedRandom1), analysis.maxValue(yEmbedRandom1));
    ui->graphEmbedRandom->replot();
// ----------------------------------------

// ---------------- Task 3 ----------------
    ui->graphTaskThree->addGraph();
    ui->graphTaskThree->graph(0)->setData(x, yMyRandom);
    ui->graphTaskThree->xAxis->setRange(0, LENGTH);
    ui->graphTaskThree->yAxis->setRange(analysis.minValue(yMyRandom) - 1, analysis.maxValue(yMyRandom) + 1);
    ui->graphTaskThree->replot();
//---------- Статистики для Task3-----------
    ui->labelMinimumValue->setText(QString::number(analysis.minValue(yMyRandom)));
    ui->labelMaximumValue->setText(QString::number(analysis.maxValue(yMyRandom)));
    ui->labelMeanValue->setText(QString::number(analysis.mean(yMyRandom)));
    ui->labelDispersionValue->setText(QString::number(analysis.dispersion(yMyRandom)));
    ui->labelStdDevValue->setText(QString::number(analysis.standartDeviation(yMyRandom)));
    ui->labelSqrtDevValue->setText(QString::number(analysis.sqrtDeviation(yMyRandom)));
    ui->labelAssymetryValue->setText(QString::number(analysis.assymetry(yMyRandom)));
    ui->labelAssymetryCoeffValue->setText(QString::number(analysis.assymetryCoeff(yMyRandom)));
    ui->labelExcessValue->setText(QString::number(analysis.excess(yMyRandom)));
    ui->labelCurtosisValue->setText(QString::number(analysis.curtosis(yMyRandom)));
    ui->labelIsStationarity->setText(analysis.isStationary(yMyRandom) ? "Стационарен" : "Не стационарен"); //использовать только для случайных процессов
// -----------------------------------------

// ---------------- Task 4 ----------------
    ui->graphTaskFourRandom1->addGraph();
    ui->graphTaskFourRandom1->graph(0)->setData(x, yEmbedRandom1);
    ui->graphTaskFourRandom1->xAxis->setRange(0, LENGTH);
    ui->graphTaskFourRandom1->yAxis->setRange(analysis.minValue(yEmbedRandom1) - 1, analysis.maxValue(yEmbedRandom1) + 1);
    ui->graphTaskFourRandom1->replot();

    ui->graphTaskFourRandom2->addGraph();
    ui->graphTaskFourRandom2->graph(0)->setData(x, yEmbedRandom2);
    ui->graphTaskFourRandom2->xAxis->setRange(0, LENGTH);
    ui->graphTaskFourRandom2->yAxis->setRange(analysis.minValue(yEmbedRandom2) - 1, analysis.maxValue(yEmbedRandom2) + 1);
    ui->graphTaskFourRandom2->replot();

    ui->graphCovar->addGraph();
    ui->graphCovar->graph(0)->setData(x, analysis.covariance(yEmbedRandom1, yEmbedRandom2));
    ui->graphCovar->xAxis->setRange(0, LENGTH);
    ui->graphCovar->yAxis->setRange(analysis.minValue(analysis.covariance(yEmbedRandom1, yEmbedRandom2)) - 1, analysis.maxValue(analysis.covariance(yEmbedRandom1, yEmbedRandom2)) + 1);
    ui->graphCovar->replot();

    ui->graphAutocovar->addGraph();
    ui->graphAutocovar->graph(0)->setData(x, analysis.autocovariance(yEmbedRandom1));
    ui->graphAutocovar->xAxis->setRange(0, LENGTH);
    ui->graphAutocovar->yAxis->setRange(analysis.minValue(analysis.autocovariance(yEmbedRandom1)), 0.2);
    ui->graphAutocovar->replot();

    QBarSet *barSet = new QBarSet("");
    QBarSeries *series = new QBarSeries();
    int M = 20;
    QVector<double> hist(20);
    double kHist = analysis.minValue(yEmbedRandom1);
    double shag = (qAbs(analysis.minValue(yEmbedRandom1)) + analysis.maxValue(yEmbedRandom1))/M;

    //0..20, 20..40, 40..60, 60..80, 80..100, ... гистограмма
    for(int i = 0; i<M; i++){
        for(int j = 0; j<yEmbedRandom1.length(); j++){
            if(yEmbedRandom1[j] >= kHist && yEmbedRandom1[j] <= (kHist+shag)){
                hist[i] = hist[i] + 1;
            }
        }
        kHist += shag;
        *barSet << hist[i];
    }

    series->append(barSet);
    QChart *chart = new QChart();
    QValueAxis *axisY = new QValueAxis();
    QValueAxis *axisX = new QValueAxis();
    axisY->setRange(0, round(analysis.maxValue(hist)));
    axisX->setRange(0, analysis.maxValue(yEmbedRandom1));
    axisX->setTickCount(M);
    axisX->setTickInterval(shag);
    chart->addAxis(axisY, Qt::AlignLeft);
    chart->addAxis(axisX, Qt::AlignBottom);
    chart->legend()->setVisible(false);
    series->attachAxis(axisY);
    series->attachAxis(axisX);
    chart->addSeries(series);
    chart->setAnimationOptions(QChart::SeriesAnimations);
    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    ui->widgetTest->setContentsMargins(0,0,0,0);
    QHBoxLayout *lay = new QHBoxLayout(ui->widgetTest);
    lay->addWidget(chartView);
// ----------------------------------------

// ---------------- Task 5 ----------------
    ui->graphTaskFiveSinMod->addGraph();
    ui->graphTaskFiveSinMod->graph(0)->setData(x, ySinMod);
    ui->graphTaskFiveSinMod->xAxis->setRange(0, LENGTH);
    ui->graphTaskFiveSinMod->yAxis->setRange(analysis.minValue(ySinMod)-1, analysis.maxValue(ySinMod)+1);
    ui->graphTaskFiveSinMod->replot();
    //QThread::sleep(1);
// ----------------------------------------

// ---------------- Task 6 ----------------
    int offsetCoeff = 200;
    QVector<double> offsetted = processing.offset(yLinInc, offsetCoeff);

    ui->graphTaskSixOrig1->addGraph();
    ui->graphTaskSixOrig1->graph(0)->setData(x, yLinInc);
    ui->graphTaskSixOrig1->xAxis->setRange(0, LENGTH);
    ui->graphTaskSixOrig1->yAxis->setRange(analysis.minValue(yLinInc)-1, analysis.maxValue(yLinInc)+1);
    ui->graphTaskSixOrig1->replot();

    ui->graphTaskSixShifted1->addGraph();
    ui->graphTaskSixShifted1->graph(0)->setData(x, offsetted);
    ui->graphTaskSixShifted1->xAxis->setRange(0, LENGTH);
    ui->graphTaskSixShifted1->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
    ui->graphTaskSixShifted1->replot();

    ui->graphTaskSixOrig2->addGraph();
    ui->graphTaskSixOrig2->graph(0)->setData(x, ySin);
    ui->graphTaskSixOrig2->xAxis->setRange(0, LENGTH);
    ui->graphTaskSixOrig2->yAxis->setRange(analysis.minValue(ySin)-1, analysis.maxValue(ySin)+1);
    ui->graphTaskSixOrig2->replot();

    offsetCoeff = 5;
    offsetted = processing.offset(ySin, 5);

    ui->graphTaskSixShifted2->addGraph();
    ui->graphTaskSixShifted2->graph(0)->setData(xTemp, offsetted);
    ui->graphTaskSixShifted2->xAxis->setRange(0, LENGTH);
    ui->graphTaskSixShifted2->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
    ui->graphTaskSixShifted2->replot();

    QVector<double> spiked = processing.spikes(ySin);

    ui->graphTaskSixSpiked->addGraph();
    ui->graphTaskSixSpiked->graph(0)->setData(x, spiked);
    ui->graphTaskSixSpiked->xAxis->setRange(0, LENGTH);
    ui->graphTaskSixSpiked->yAxis->setRange(analysis.minValue(spiked)-1, analysis.maxValue(spiked)+1);
    ui->graphTaskSixSpiked->replot();
// ----------------------------------------

// ---------------- Task 7 ----------------
    offsetCoeff = 5;
    offsetted = processing.offset(ySin, offsetCoeff);

    ui->graphTaskSevenShifted->addGraph();
    ui->graphTaskSevenShifted->graph(0)->setData(x, offsetted);
    ui->graphTaskSevenShifted->xAxis->setRange(0, LENGTH);
    ui->graphTaskSevenShifted->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
    ui->graphTaskSevenShifted->replot();

    QVector<double> antiShifted = processing.antiShift(offsetted);

    ui->graphTaskSevenAntiShifted->addGraph();
    ui->graphTaskSevenAntiShifted->graph(0)->setData(x, antiShifted);
    ui->graphTaskSevenAntiShifted->xAxis->setRange(0, LENGTH);
    ui->graphTaskSevenAntiShifted->yAxis->setRange(analysis.minValue(antiShifted)-offsetCoeff, analysis.maxValue(antiShifted)+offsetCoeff);
    ui->graphTaskSevenAntiShifted->replot();

    ui->graphTaskSevenSpiked->addGraph();
    ui->graphTaskSevenSpiked->graph(0)->setData(x, spiked);
    ui->graphTaskSevenSpiked->xAxis->setRange(0, LENGTH);
    //ui->graphTaskSevenSpiked->yAxis->setRange(analysis.minValue(spiked)-1, analysis.maxValue(spiked)+1);
    ui->graphTaskSevenSpiked->yAxis->setRange(-11, 11);
    ui->graphTaskSevenSpiked->replot();

    QVector<double> antispiked = processing.antiSpike(spiked);

    ui->graphTaskSevenAntiSpiked->addGraph();
    ui->graphTaskSevenAntiSpiked->graph(0)->setData(x, antispiked);
    ui->graphTaskSevenAntiSpiked->xAxis->setRange(0, LENGTH);
    ui->graphTaskSevenAntiSpiked->yAxis->setRange(analysis.minValue(antispiked)-1, analysis.maxValue(antispiked)+1);
    ui->graphTaskSevenAntiSpiked->replot();


// ----------------------------------------

    ui->widget_7->addGraph();
    ui->widget_7->graph(0)->setData(xTemp, ySmoothed);
    ui->widget_7->xAxis->setRange(0, LENGTH);
    ui->widget_7->yAxis->setRange(0, 1000);
    ui->widget_7->replot();


// ---------- Очистка графиков ----------
    ui->graphLinInc->clearGraphs();
    ui->graphLinDec->clearGraphs();
    ui->graphExpInc->clearGraphs();
    ui->graphExpDec->clearGraphs();
    ui->graphCombinated->clearGraphs();
    ui->graphMyRandom->clearGraphs();
    ui->graphEmbedRandom->clearGraphs();
    ui->graphTaskThree->clearGraphs();
    ui->graphTaskFourRandom1->clearGraphs();
    ui->graphTaskFourRandom2->clearGraphs();
    ui->graphCovar->clearGraphs();
    ui->graphAutocovar->clearGraphs();
    ui->graphTaskFiveSinMod->clearGraphs();
    ui->graphTaskSixOrig1->clearGraphs();
    ui->graphTaskSixShifted1->clearGraphs();
    ui->graphTaskSixOrig2->clearGraphs();
    ui->graphTaskSixShifted2->clearGraphs();
    ui->graphTaskSixSpiked->clearGraphs();
    ui->graphTaskSevenShifted->clearGraphs();
    ui->graphTaskSevenAntiShifted->clearGraphs();
    ui->graphTaskSevenSpiked->clearGraphs();
    ui->graphTaskSevenAntiSpiked->clearGraphs();
    ui->graphTaskSevenShifted->clearGraphs();
// --------------------------------------
}
