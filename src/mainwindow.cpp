#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "modeling.h"
#include "analysis.h"
#include "inou.h"
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

#include <fstream>

using namespace QtCharts;
using namespace std;

const int LENGTH = 1000;
const double DELTA = 50.0;

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_pushButton_clicked(){
    modeling  modeling;
    processing processing;
    analysis analysis;

    double k = -1;

    //double beta = 11;
    //double alpha = -0.01;

    int amplitude[] = {10, 100, 15};
    int frequency[] = {4, 37, 173};
    float randomCoeff = (float)ui->sbCoeff->value();

    QVector<double> x(LENGTH),             xLarge(3*LENGTH),
                    yLinInc,               yLinDec,                 yExpInc,
                    yExpDec,               yHarmonic,               yPolyharmonic,
                    yCombinated,           yCardio,
                    yMyRandom(LENGTH),     yEmbedRandom1(LENGTH),   yEmbedRandom2(LENGTH);

    yLinInc.append(modeling.linearIncrement(LENGTH, k));
    yLinDec.append(modeling.linearDecrement(LENGTH, k));
    yExpInc.append(modeling.exponentialIncrement(LENGTH));
    yExpDec.append(modeling.exponentialDecrement(LENGTH));
    yHarmonic.append(modeling.harmonic(LENGTH, 10, 4));
    yPolyharmonic.append(modeling.polyharmonic(LENGTH, amplitude, frequency));
    yCardio.append(modeling.cardiogram());

    for(int X=0; X<LENGTH; X++){
        x[X] = X;
        yMyRandom[X] = modeling.randomGenerator(randomCoeff);
        yEmbedRandom1[X] = modeling.embedRandom();
        yEmbedRandom2[X] = modeling.embedRandom();
    }

    for(int i = 0; i<3000; i++){
        xLarge[i] = i;
    }

//------Для первого задания ------
    yCombinated.append(yLinInc);
    yCombinated.append(yLinDec);
    yCombinated.append(yExpInc);
//--------------------------------

    double fc = 0.3;
    double fc2 = 0.6;
    int m = 64;
    int Loper = 2 * m + 1;
    QVector<double> lpfWeights = analysis.lowpassFilterPotter(fc, m);
    QVector<double> lpfWeights2 = analysis.lowpassFilterPotter(fc2, m);
    QVector<double> hpfWeights;
    QVector<double> pfWeights;
    QVector<double> rfWeights;

// ---------------- Рисование ----------------
    int currentTask = ui->tabWidget->currentIndex() + 1;
    switch(currentTask){
        case 1:
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
            break;
        case 2:
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
            break;
        case 3:
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
            break;
        case 4:{
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
            ui->graphAutocorr->graph(0)->setData(x, analysis.autocorrelation(yEmbedRandom1));
            ui->graphAutocorr->xAxis->setRange(0, LENGTH);
            ui->graphAutocorr->yAxis->setRange(analysis.minValue(analysis.autocorrelation(yEmbedRandom1)), 0.2);
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
            break;
        }
        case 5:{
        // ---------------- Task 5 ----------------
            ui->graphTaskFiveHarmo->addGraph();
            ui->graphTaskFiveHarmo->graph(0)->setData(x, yHarmonic);
            ui->graphTaskFiveHarmo->xAxis->setRange(0, LENGTH);
            ui->graphTaskFiveHarmo->yAxis->setRange(analysis.minValue(yHarmonic)-1, analysis.maxValue(yHarmonic)+1);
            ui->graphTaskFiveHarmo->replot();

            ui->graphTaskFivePolyharmo->addGraph();
            ui->graphTaskFivePolyharmo->graph(0)->setData(x, yPolyharmonic);
            ui->graphTaskFivePolyharmo->xAxis->setRange(0, LENGTH);
            ui->graphTaskFivePolyharmo->yAxis->setRange(analysis.minValue(yPolyharmonic)-1, analysis.maxValue(yPolyharmonic)+1);
            ui->graphTaskFivePolyharmo->replot();

            ui->graphTaskFiveHarmoAutocorr->addGraph();
            ui->graphTaskFiveHarmoAutocorr->graph(0)->setData(x, analysis.autocorrelation(yHarmonic));
            ui->graphTaskFiveHarmoAutocorr->xAxis->setRange(0, LENGTH);
            ui->graphTaskFiveHarmoAutocorr->yAxis->setRange(analysis.minValue(analysis.autocorrelation(yHarmonic))-1, analysis.maxValue(analysis.autocorrelation(yHarmonic))+1);
            ui->graphTaskFiveHarmoAutocorr->replot();

            ui->graphTaskFivePolyAutocorr->addGraph();
            ui->graphTaskFivePolyAutocorr->graph(0)->setData(x, analysis.autocorrelation(yPolyharmonic));
            ui->graphTaskFivePolyAutocorr->xAxis->setRange(0, LENGTH);
            ui->graphTaskFivePolyAutocorr->yAxis->setRange(analysis.minValue(analysis.autocorrelation(yPolyharmonic))-1, analysis.maxValue(analysis.autocorrelation(yPolyharmonic))+1);
            ui->graphTaskFivePolyAutocorr->replot();

            ui->graphTaskFiveCrosscorr->addGraph();
            ui->graphTaskFiveCrosscorr->graph(0)->setData(x, analysis.covariance(yHarmonic, yPolyharmonic));
            ui->graphTaskFiveCrosscorr->xAxis->setRange(0, LENGTH);
            ui->graphTaskFiveCrosscorr->yAxis->setRange(analysis.minValue(analysis.covariance(yHarmonic, yPolyharmonic))-1, analysis.maxValue(analysis.covariance(yHarmonic, yPolyharmonic))+1);
            ui->graphTaskFiveCrosscorr->replot();

            QBarSet *barSet2 = new QBarSet("");
            QBarSeries *series2 = new QBarSeries();
            int M = 20;
            QVector<double> hist2(20);
            double kHist = analysis.minValue(yHarmonic);
            double shag = (qAbs(analysis.minValue(yHarmonic)) + analysis.maxValue(yHarmonic))/M;

            //0..20, 20..40, 40..60, 60..80, 80..100, ... гистограмма
            for(int i = 0; i<M; i++){
                for(int j = 0; j<yHarmonic.length(); j++){
                    if(yHarmonic[j] >= kHist && yHarmonic[j] <= (kHist+shag)){
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
            axisX2->setRange(0, analysis.maxValue(yHarmonic));
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
            break;
        }
        case 6:{
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
            ui->graphTaskSixOrig2->graph(0)->setData(x, yHarmonic);
            ui->graphTaskSixOrig2->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixOrig2->yAxis->setRange(analysis.minValue(yHarmonic)-1, analysis.maxValue(yHarmonic)+1);
            ui->graphTaskSixOrig2->replot();

            offsetCoeff = 5;
            offsetted = processing.offset(yHarmonic, 5);

            ui->graphTaskSixShifted2->addGraph();
            ui->graphTaskSixShifted2->graph(0)->setData(x, offsetted);
            ui->graphTaskSixShifted2->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixShifted2->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
            ui->graphTaskSixShifted2->replot();

            QVector<double> spiked = processing.spikes(yHarmonic);

            ui->graphTaskSixSpiked->addGraph();
            ui->graphTaskSixSpiked->graph(0)->setData(x, spiked);
            ui->graphTaskSixSpiked->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixSpiked->yAxis->setRange(analysis.minValue(spiked)-1, analysis.maxValue(spiked)+1);
            ui->graphTaskSixSpiked->replot();
        // ----------------------------------------
            break;
        }
        case 7:{
        // ---------------- Task 7 ----------------
            int offsetCoeff = 5;
            QVector<double> offsetted = processing.offset(yHarmonic, offsetCoeff);
            QVector<double> spiked = processing.spikes(yHarmonic);

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
            break;
        }
        case 8:{
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
            break;
        }
        case 9:{
        // ---------------- Task 9 ----------------
            QVector<double> taskNine1 = processing.pdfTaskNine(yHarmonic, 1);
            QVector<double> taskNine10 = processing.pdfTaskNine(yHarmonic, 10);
            QVector<double> taskNine100 = processing.pdfTaskNine(yHarmonic, 100);
            QVector<double> taskNine1000 = processing.pdfTaskNine(yHarmonic, 1000);


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
            break;
        }
        case 10:{
        // --------------- Task 10 ----------------
            QVector<double> dataFile = inou().load("C:/data.dat");
            QVector<double> xfile(dataFile.length());
            for (int i=0; i<xfile.length(); i++) { xfile[i] = i; }


            ui->graphFromFile->addGraph();
            ui->graphFromFile->graph(0)->setData(xfile, dataFile);
            ui->graphFromFile->xAxis->setRange(0, xfile.length());
            ui->graphFromFile->yAxis->setRange(analysis.minValue(dataFile)-10, analysis.maxValue(dataFile)+10);
            ui->graphFromFile->replot();

            QVector<double> fouramp = analysis.fourierAmplitude(dataFile);
            QVector<double> fourierFreq = analysis.calculateFrequency(0.002, dataFile.length());


            ui->graphFourierAmplitude->addGraph();
            ui->graphFourierAmplitude->graph(0)->setData(fourierFreq, fouramp);
            ui->graphFourierAmplitude->xAxis->setRange(0, 250);
            ui->graphFourierAmplitude->yAxis->setRange(analysis.minValue(fouramp), analysis.maxValue(fouramp)+5);
            ui->graphFourierAmplitude->replot();

            QVector<double> fourspec = analysis.fourierSpectrum(dataFile, 0.91);

            ui->graphFourierSpectrum->addGraph();
            ui->graphFourierSpectrum->graph(0)->setData(fourierFreq, fourspec);
            ui->graphFourierSpectrum->xAxis->setRange(0, 250);
            ui->graphFourierSpectrum->yAxis->setRange(analysis.minValue(fourspec), analysis.maxValue(fourspec)+6);
            ui->graphFourierSpectrum->replot();
        // ----------------------------------------
            break;
        }
        case 11:
        // --------------- Task 11 ----------------
            ui->graphCardio->addGraph();
            ui->graphCardio->graph(0)->setData(x, yCardio);
            ui->graphCardio->xAxis->setRange(0, x.length());
            ui->graphCardio->yAxis->setRange(analysis.minValue(yCardio)-1, analysis.maxValue(yCardio)+1);
            ui->graphCardio->replot();
        // ----------------------------------------
            break;
        case 12:{
            ui->graphLPF->addGraph();
            ui->graphLPF->graph(0)->setData(x, lpfWeights);
            ui->graphLPF->xAxis->setRange(0, 2*m);
            ui->graphLPF->yAxis->setRange(analysis.minValue(lpfWeights)-0.01, analysis.maxValue(lpfWeights)+0.01);
            ui->graphLPF->replot();

            for(int i = 0; i<=Loper-1; i++){ i==m ? hpfWeights.append(1 - lpfWeights[i]) : hpfWeights.append(-lpfWeights[i]); }

            ui->graphHPF->addGraph();
            ui->graphHPF->graph(0)->setData(x, hpfWeights);
            ui->graphHPF->xAxis->setRange(0, 2*m);
            ui->graphHPF->yAxis->setRange(analysis.minValue(hpfWeights)-0.01, analysis.maxValue(hpfWeights)+0.1);
            ui->graphHPF->replot();

            for(int i = 0; i<=Loper-1; i++){ pfWeights.append(lpfWeights2[i] - lpfWeights[i]); }

            ui->graphPF->addGraph();
            ui->graphPF->graph(0)->setData(x, pfWeights);
            ui->graphPF->xAxis->setRange(0, 2*m);
            ui->graphPF->yAxis->setRange(analysis.minValue(pfWeights)-0.01, analysis.maxValue(pfWeights)+0.01);
            ui->graphPF->replot();

            for(int i = 0; i<=Loper-1; i++){ i==m ? rfWeights.append(1 + lpfWeights[i] - lpfWeights2[i]) : rfWeights.append(lpfWeights[i] - lpfWeights2[i]); }

            ui->graphRF->addGraph();
            ui->graphRF->graph(0)->setData(x, rfWeights);
            ui->graphRF->xAxis->setRange(0, 2*m);
            ui->graphRF->yAxis->setRange(analysis.minValue(rfWeights)-0.01, analysis.maxValue(rfWeights)+0.01);
            ui->graphRF->replot();
            break;
        }
        case 13:{
            QVector<double> achxLP = analysis.fourierAmplitude(lpfWeights);
            QVector<double> achxHP = analysis.fourierAmplitude(hpfWeights);
            QVector<double> achxPF = analysis.fourierAmplitude(pfWeights);
            QVector<double> achxRF = analysis.fourierAmplitude(rfWeights);

            for(int i = 0; i<m; i++){
                achxLP[i] *= 2*m;
                achxHP[i] *= 2*m;
                achxPF[i] *= 2*m;
                achxRF[i] *= 2*m;
            }

            ui->graphAfcLp->addGraph();
            ui->graphAfcLp->graph(0)->setData(x, achxLP);
            ui->graphAfcLp->xAxis->setRange(0, m);
            ui->graphAfcLp->yAxis->setRange(0, analysis.maxValue(achxLP)+analysis.maxValue(achxLP)/10);
            ui->graphAfcLp->replot();

            ui->graphAfcHp->addGraph();
            ui->graphAfcHp->graph(0)->setData(x, achxHP);
            ui->graphAfcHp->xAxis->setRange(0, m);
            ui->graphAfcHp->yAxis->setRange(0, analysis.maxValue(achxHP)+analysis.maxValue(achxHP)/10);
            ui->graphAfcHp->replot();

            ui->graphAfcPf->addGraph();
            ui->graphAfcPf->graph(0)->setData(x, achxPF);
            ui->graphAfcPf->xAxis->setRange(0, m);
            ui->graphAfcPf->yAxis->setRange(0, analysis.maxValue(achxPF)+analysis.maxValue(achxPF)/10);
            ui->graphAfcPf->replot();

            ui->graphAfcRf->addGraph();
            ui->graphAfcRf->graph(0)->setData(x, achxRF);
            ui->graphAfcRf->xAxis->setRange(0, m);
            ui->graphAfcRf->yAxis->setRange(0, analysis.maxValue(achxRF)+analysis.maxValue(achxRF)/10);
            ui->graphAfcRf->replot();
            break;
        }
        case 14:{
            QVector<double> dataFile = inou().load("C:/data.dat");
            int N = dataFile.length();
            int M = lpfWeights.length();
            QVector<double> filteredLp(N+M);
            int tempN = N + M-1;

            for(auto i(0); i<tempN; ++i){
                int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
                int const jmx = (i < N - 1)  ? i           : N-1;
                for(auto j(jmn); j<=jmx; ++j){
                    filteredLp[i] += dataFile[j] * lpfWeights[i-j];
                }
            }

            for(int i(0); i<=M/2; i++){
                filteredLp.pop_front();
                filteredLp.pop_back();
            }

            ui->graphFilteredLp->addGraph();
            ui->graphFilteredLp->graph(0)->setData(x, filteredLp);
            ui->graphFilteredLp->xAxis->setRange(0, 1000);
            ui->graphFilteredLp->yAxis->setRange(analysis.minValue(filteredLp), analysis.maxValue(filteredLp)+analysis.maxValue(filteredLp)/10);
            ui->graphFilteredLp->replot();

            ui->graphFilteredLpAmpl->addGraph();
            ui->graphFilteredLpAmpl->graph(0)->setData(x, analysis.fourierAmplitude(filteredLp));
            ui->graphFilteredLpAmpl->xAxis->setRange(0, 400);
            ui->graphFilteredLpAmpl->yAxis->setRange(analysis.minValue(analysis.fourierAmplitude(filteredLp)), analysis.maxValue(analysis.fourierAmplitude(filteredLp))+analysis.maxValue(analysis.fourierAmplitude(filteredLp))/10);
            ui->graphFilteredLpAmpl->replot();

            M = hpfWeights.length();
            QVector<double> filteredHp(N+M);
            tempN = N + M-1;

            for(auto i(0); i<tempN; ++i){
                int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
                int const jmx = (i < N - 1)  ? i           : N-1;
                for(auto j(jmn); j<=jmx; ++j){
                    filteredHp[i] += dataFile[j] * hpfWeights[i-j];
                }
            }

            for(int i(0); i<=M/2; i++){
                filteredHp.pop_front();
                filteredHp.pop_back();
            }

            ui->graphFilteredHp->addGraph();
            ui->graphFilteredHp->graph(0)->setData(x, filteredHp);
            ui->graphFilteredHp->xAxis->setRange(0, 1000);
            ui->graphFilteredHp->yAxis->setRange(analysis.minValue(filteredHp), analysis.maxValue(filteredHp)+analysis.maxValue(filteredHp)/10);
            ui->graphFilteredHp->replot();

            ui->graphFilteredHpAmpl->addGraph();
            ui->graphFilteredHpAmpl->graph(0)->setData(x, analysis.fourierAmplitude(filteredHp));
            ui->graphFilteredHpAmpl->xAxis->setRange(0, 400);
            ui->graphFilteredHpAmpl->yAxis->setRange(analysis.minValue(analysis.fourierAmplitude(filteredHp)), analysis.maxValue(analysis.fourierAmplitude(filteredHp))+analysis.maxValue(analysis.fourierAmplitude(filteredHp))/10);
            ui->graphFilteredHpAmpl->replot();

            M = rfWeights.length();
            QVector<double> filteredRf(N+M);
            tempN = N + M-1;

            for(auto i(0); i<tempN; ++i){
                int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
                int const jmx = (i < N - 1)  ? i           : N-1;
                for(auto j(jmn); j<=jmx; ++j){
                    filteredRf[i] += dataFile[j] * rfWeights[i-j];
                }
            }

            for(int i(0); i<=M/2; i++){
                filteredRf.pop_front();
                filteredRf.pop_back();
            }

            ui->graphFilteredPf->addGraph();
            ui->graphFilteredPf->graph(0)->setData(x, filteredRf);
            ui->graphFilteredPf->xAxis->setRange(0, 1000);
            ui->graphFilteredPf->yAxis->setRange(analysis.minValue(filteredRf), analysis.maxValue(filteredRf)+analysis.maxValue(filteredRf)/10);
            ui->graphFilteredPf->replot();

            ui->graphFilteredPfAmpl->addGraph();
            ui->graphFilteredPfAmpl->graph(0)->setData(x, analysis.fourierAmplitude(filteredRf));
            ui->graphFilteredPfAmpl->xAxis->setRange(0, 400);
            ui->graphFilteredPfAmpl->yAxis->setRange(analysis.minValue(analysis.fourierAmplitude(filteredRf)), analysis.maxValue(analysis.fourierAmplitude(filteredRf))+analysis.maxValue(analysis.fourierAmplitude(filteredRf))/10);
            ui->graphFilteredPfAmpl->replot();
            break;
        }
        case 15:{
        // --------------- Task 12 ----------------
            QVector<double> wav = inou().loadWave("../dojdik.wav");
            QVector<double> wavFourier = analysis.fourierAmplitude(wav);
            QVector<double> xWav(wav.length());
            auto min = analysis.minValue(wav);
            auto max = analysis.maxValue(wav);
            for (int i(0); i<wav.length(); i++) { xWav[i] = i; }

            ui->graphWavFile->addGraph();
            ui->graphWavFile->graph(0)->setData(xWav, wav);
            ui->graphWavFile->xAxis->setRange(0, xWav.length());
            ui->graphWavFile->yAxis->setRange(min, max+max/10);
            ui->graphWavFile->replot();

            min = analysis.minValue(wavFourier);
            max = analysis.maxValue(wavFourier);

            ui->graphWavFileFourier->addGraph();
            ui->graphWavFileFourier->graph(0)->setData(xWav, wavFourier);
            ui->graphWavFileFourier->xAxis->setRange(0, wavFourier.length());
            ui->graphWavFileFourier->yAxis->setRange(min, max+max/10);
            ui->graphWavFileFourier->replot();

            //for (int i(0); i<wav.length(); i++) { wav[i] *= 4; }

            min = analysis.minValue(wav);
            max = analysis.maxValue(wav);

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(min, max+max/10);
            ui->graphWavFileLoud->replot();
        // ----------------------------------------
            break;
        }
        case 16:{
            //ПРОДОЛЖИТЬ
            QVector<double> wav = inou().loadWave("../pole.wav");
            QVector<double> Fonem;

            for(int i(4000); i<16000; i++){
                Fonem.append(wav[i]);
            }
            QVector<double> xWav;
            for (int i(0); i<wav.length(); i++) {
                xWav.append(i);
            }


            QVector<double> FonemFourier = analysis.fourierAmplitude(Fonem);

            ui->graph08->addGraph();
            ui->graph08->graph(0)->setData(xWav, Fonem);
            ui->graph08->xAxis->setRange(0, Fonem.length());
            ui->graph08->yAxis->setRange(analysis.minValue(Fonem)-analysis.maxValue(Fonem)/10, analysis.maxValue(Fonem)+analysis.maxValue(Fonem)/10);
            ui->graph08->replot();

            int N = Fonem.length();
            //double dt = 1 / FonemFourier.length(); //
            //if(dt = 0){ throw std::exception(); }
            double f_gr = N/2;
            double df = (2 * f_gr) / N;
            //float df = 22050/FonemFourier.length();

            QVector<double> xWav1;
            for (int i= 0; i<FonemFourier.length(); i++) {
                xWav1.append((double)(i * df));
            }
            inou().exportWave(Fonem, Fonem.length(), "../FonemTest.wav", 1);

            ui->graph08_2->addGraph();
            ui->graph08_2->graph(0)->setData(xWav1, FonemFourier);
            ui->graph08_2->xAxis->setRange(0, FonemFourier.length());
            ui->graph08_2->setInteraction(QCP::iRangeZoom,true);
            ui->graph08_2->setInteraction(QCP::iRangeDrag,true);
            ui->graph08_2->yAxis->setRange(analysis.minValue(FonemFourier)-analysis.maxValue(FonemFourier)/10, analysis.maxValue(FonemFourier)+analysis.maxValue(FonemFourier)/10);
            ui->graph08_2->replot();

            /*ui->graph08_3->addGraph();
            ui->graph08_3->graph(0)->setData(xWav, foo);
            ui->graph08_3->xAxis->setRange(0, foo.length());
            ui->graph08_3->yAxis->setRange(analysis.minValue(foo)-analysis.maxValue(foo)/10, analysis.maxValue(foo)+analysis.maxValue(foo)/10);
            ui->graph08_3->replot();*/

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(analysis.minValue(wav)-analysis.maxValue(wav)/10, analysis.maxValue(wav)+analysis.maxValue(wav)/10);
            ui->graphWavFileLoud->replot();

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(analysis.minValue(wav)-analysis.maxValue(wav)/10, analysis.maxValue(wav)+analysis.maxValue(wav)/10);
            ui->graphWavFileLoud->replot();

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(analysis.minValue(wav)-analysis.maxValue(wav)/10, analysis.maxValue(wav)+analysis.maxValue(wav)/10);
            ui->graphWavFileLoud->replot();
        }


        default:
            return;
    }

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
