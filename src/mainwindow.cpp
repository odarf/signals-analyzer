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
#include <cstdlib>

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
    //float value = 0;
    int max = 1;
    int min = -1;
    //int c = 1;
    //QRandomGenerator *rnd = QRandomGenerator::global();
    //value = rnd->bounded(min, max) * (c * max - c * min) + c * min;
    return static_cast<float>(rand() * (1.0 / (static_cast<float>(RAND_MAX) + 1.0)) * (max - min) + min);
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
                ySinMod[X] += amplitude[count] * sin(2 * 3.14 * frequency[count] * X);// * delta_t) + DELTA; //график со смещением с одним значением амплитуды и частоты
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
    ui->graphCovar->yAxis->setRange(analysis.minValue(analysis.covariance(yEmbedRandom1, yEmbedRandom2)), analysis.maxValue(analysis.covariance(yEmbedRandom1, yEmbedRandom2)));
    ui->graphCovar->replot();

    ui->graphAutocorr->addGraph();
    ui->graphAutocorr->graph(0)->setData(x, analysis.autocovariance(yEmbedRandom1));
    ui->graphAutocorr->xAxis->setRange(0, LENGTH);
    ui->graphAutocorr->yAxis->setRange(analysis.minValue(analysis.autocovariance(yEmbedRandom1)), 0.2);
    ui->graphAutocorr->replot();

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
    ui->graphTaskFiveHarmo->addGraph();
    ui->graphTaskFiveHarmo->graph(0)->setData(x, ySin);
    ui->graphTaskFiveHarmo->xAxis->setRange(0, LENGTH);
    ui->graphTaskFiveHarmo->yAxis->setRange(analysis.minValue(ySin)-1, analysis.maxValue(ySin)+1);
    ui->graphTaskFiveHarmo->replot();

    ui->graphTaskFivePolyharmo->addGraph();
    ui->graphTaskFivePolyharmo->graph(0)->setData(x, ySinMod);
    ui->graphTaskFivePolyharmo->xAxis->setRange(0, LENGTH);
    ui->graphTaskFivePolyharmo->yAxis->setRange(analysis.minValue(ySinMod)-1, analysis.maxValue(ySinMod)+1);
    ui->graphTaskFivePolyharmo->replot();

    ui->graphTaskFiveHarmoAutocorr->addGraph();
    ui->graphTaskFiveHarmoAutocorr->graph(0)->setData(x, analysis.autocovariance(ySin));
    ui->graphTaskFiveHarmoAutocorr->xAxis->setRange(0, LENGTH);
    ui->graphTaskFiveHarmoAutocorr->yAxis->setRange(analysis.minValue(analysis.autocovariance(ySin))-1, analysis.maxValue(analysis.autocovariance(ySin))+1);
    ui->graphTaskFiveHarmoAutocorr->replot();

    ui->graphTaskFivePolyAutocorr->addGraph();
    ui->graphTaskFivePolyAutocorr->graph(0)->setData(x, analysis.autocovariance(ySinMod));
    ui->graphTaskFivePolyAutocorr->xAxis->setRange(0, LENGTH);
    ui->graphTaskFivePolyAutocorr->yAxis->setRange(analysis.minValue(analysis.autocovariance(ySinMod))-1, analysis.maxValue(analysis.autocovariance(ySinMod))+1);
    ui->graphTaskFivePolyAutocorr->replot();

    ui->graphTaskFiveCrosscorr->addGraph();
    ui->graphTaskFiveCrosscorr->graph(0)->setData(x, analysis.covariance(ySin, ySinMod));
    ui->graphTaskFiveCrosscorr->xAxis->setRange(0, LENGTH);
    ui->graphTaskFiveCrosscorr->yAxis->setRange(analysis.minValue(analysis.covariance(ySin, ySinMod))-1, analysis.maxValue(analysis.covariance(ySin, ySinMod))+1);
    ui->graphTaskFiveCrosscorr->replot();

    QBarSet *barSet2 = new QBarSet("");
    QBarSeries *series2 = new QBarSeries();
    M = 20;
    QVector<double> hist2(20);
    kHist = analysis.minValue(ySin);
    shag = (qAbs(analysis.minValue(ySin)) + analysis.maxValue(ySin))/M;

    //0..20, 20..40, 40..60, 60..80, 80..100, ... гистограмма
    for(int i = 0; i<M; i++){
        for(int j = 0; j<ySin.length(); j++){
            if(ySin[j] >= kHist && ySin[j] <= (kHist+shag)){
                hist2[i] = hist2[i] + 1;
            }
        }
        kHist += shag;
        *barSet2 << hist2[i];
    }

    series2->append(barSet2);
    QChart *chart2 = new QChart();
    QValueAxis *axisY2 = new QValueAxis();
    QValueAxis *axisX2 = new QValueAxis();
    axisY2->setRange(0, round(analysis.maxValue(hist2)));
    axisX2->setRange(0, analysis.maxValue(ySin));
    axisX2->setTickCount(M);
    axisX2->setTickInterval(shag);
    chart2->addAxis(axisY2, Qt::AlignLeft);
    chart2->addAxis(axisX2, Qt::AlignBottom);
    chart2->legend()->setVisible(false);
    series2->attachAxis(axisY2);
    series2->attachAxis(axisX2);
    chart2->addSeries(series2);
    chart2->setAnimationOptions(QChart::SeriesAnimations);
    QChartView *chartView2 = new QChartView(chart2);
    chartView2->setRenderHint(QPainter::Antialiasing);
    ui->graphTaskFiveDensity->setContentsMargins(0,0,0,0);
    QHBoxLayout *lay2 = new QHBoxLayout(ui->graphTaskFiveDensity);
    lay2->addWidget(chartView2);
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
    ui->graphTaskSevenSpiked->yAxis->setRange(-100, 100);
    ui->graphTaskSevenSpiked->replot();

    QVector<double> antispiked = processing.antiSpike(spiked);

    ui->graphTaskSevenAntiSpiked->addGraph();
    ui->graphTaskSevenAntiSpiked->graph(0)->setData(x, antispiked);
    ui->graphTaskSevenAntiSpiked->xAxis->setRange(0, LENGTH);
    ui->graphTaskSevenAntiSpiked->yAxis->setRange(-100, 100);
    ui->graphTaskSevenAntiSpiked->replot();


// ----------------------------------------

// ---------------- Task 8 ----------------
    QVector<double> trAdRa = processing.trendAddRandom(yLinInc);
    QVector<double> pot = processing.pickoutTrend(trAdRa);
    QVector<double> antr = processing.antiTrend(trAdRa);

    ui->graphTaskEightTrendAdd->addGraph();
    ui->graphTaskEightTrendAdd->graph(0)->setData(x, trAdRa);
    ui->graphTaskEightTrendAdd->xAxis->setRange(0, LENGTH);
    ui->graphTaskEightTrendAdd->yAxis->setRange(analysis.minValue(trAdRa)-1, analysis.maxValue(trAdRa)+1);
    ui->graphTaskEightTrendAdd->replot();

    ui->graphTaskEightPickTrend->addGraph();
    ui->graphTaskEightPickTrend->graph(0)->setData(x, pot);
    ui->graphTaskEightPickTrend->xAxis->setRange(0, LENGTH);
    ui->graphTaskEightPickTrend->yAxis->setRange(analysis.minValue(pot)-1, analysis.maxValue(pot)+1);
    ui->graphTaskEightPickTrend->replot();

    ui->graphTaskEightAntiTrend->addGraph();
    ui->graphTaskEightAntiTrend->graph(0)->setData(x, antr);
    ui->graphTaskEightAntiTrend->xAxis->setRange(0, LENGTH);
    ui->graphTaskEightAntiTrend->yAxis->setRange(analysis.minValue(antr)-1, analysis.maxValue(antr)+1);
    ui->graphTaskEightAntiTrend->replot();

    QVector<double> stdDev = processing.pdfTaskEight();

    ui->graphTaskEightStdDev->addGraph();
    ui->graphTaskEightStdDev->graph(0)->setData(x, stdDev);
    ui->graphTaskEightStdDev->xAxis->setRange(0, stdDev.length()-1);
    ui->graphTaskEightStdDev->yAxis->setRange(0, analysis.maxValue(stdDev)+0.05);
    ui->graphTaskEightStdDev->replot();

    QVector<double> stdDevDiv = processing.pdfTaskEight2(stdDev);

    ui->graphTaskEightStdDev_2->addGraph();
    ui->graphTaskEightStdDev_2->graph(0)->setData(x, stdDevDiv);
    ui->graphTaskEightStdDev_2->xAxis->setRange(0, stdDevDiv.length()-1);
    ui->graphTaskEightStdDev_2->yAxis->setRange(analysis.minValue(stdDevDiv)-0.5, analysis.maxValue(stdDevDiv)+0.5);
    ui->graphTaskEightStdDev_2->replot();

// ----------------------------------------

// ---------------- Task 8 ----------------
    QVector<double> taskNine1 = processing.pdfTaskNine(ySin, 1);
    QVector<double> taskNine10 = processing.pdfTaskNine(ySin, 10);
    QVector<double> taskNine100 = processing.pdfTaskNine(ySin, 100);
    QVector<double> taskNine1000 = processing.pdfTaskNine(ySin, 1000);


    ui->graphTaskNine1->addGraph();
    ui->graphTaskNine1->graph(0)->setData(x, taskNine1);
    ui->graphTaskNine1->xAxis->setRange(0, taskNine1.length()+1);
    ui->graphTaskNine1->yAxis->setRange(analysis.minValue(taskNine1), analysis.maxValue(taskNine1));
    ui->graphTaskNine1->replot();

    ui->graphTaskNine10->addGraph();
    ui->graphTaskNine10->graph(0)->setData(x, taskNine10);
    ui->graphTaskNine10->xAxis->setRange(0, taskNine10.length()+1);
    ui->graphTaskNine10->yAxis->setRange(analysis.minValue(taskNine10), analysis.maxValue(taskNine10));
    ui->graphTaskNine10->replot();

    ui->graphTaskNine100->addGraph();
    ui->graphTaskNine100->graph(0)->setData(x, taskNine100);
    ui->graphTaskNine100->xAxis->setRange(0, taskNine100.length()+1);
    ui->graphTaskNine100->yAxis->setRange(analysis.minValue(taskNine100), analysis.maxValue(taskNine100));
    ui->graphTaskNine100->replot();

    ui->graphTaskNine1000->addGraph();
    ui->graphTaskNine1000->graph(0)->setData(x, taskNine1000);
    ui->graphTaskNine1000->xAxis->setRange(0, taskNine1000.length()+1);
    ui->graphTaskNine1000->yAxis->setRange(analysis.minValue(taskNine1000), analysis.maxValue(taskNine1000));
    ui->graphTaskNine1000->replot();

// ----------------------------------------

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
    ui->graphAutocorr->clearGraphs();
    ui->graphTaskFiveHarmo->clearGraphs();
    ui->graphTaskFivePolyharmo->clearGraphs();
    ui->graphTaskFiveHarmoAutocorr->clearGraphs();
    ui->graphTaskFivePolyAutocorr->clearGraphs();
    ui->graphTaskFiveCrosscorr->clearGraphs();
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
