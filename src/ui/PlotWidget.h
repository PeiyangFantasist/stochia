#pragma once

#include <QColor>
#include <QPointF>
#include <QWidget>

#include <vector>

class PlotWidget final : public QWidget {
public:
    enum class Style { Curve, Steps, Bars };

    struct Series {
        QString name;
        QColor color;
        Style style = Style::Curve;
        std::vector<QPointF> points;
        bool filled = false;
    };

    explicit PlotWidget(QWidget* parent = nullptr);

    void setSeries(std::vector<Series> series);
    void setAxisLabels(QString xLabel, QString yLabel);
    void setEmptyMessage(QString message);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void leaveEvent(QEvent* event) override;

private:
    QRectF plotRect() const;
    QPointF mapPoint(const QPointF& point, const QRectF& rect) const;
    void updateBounds();
    QString formatNumber(double value) const;

    std::vector<Series> series_;
    QString xLabel_ = QStringLiteral("x");
    QString yLabel_ = QStringLiteral("f(x)");
    QString emptyMessage_ = QStringLiteral("选择一个随机变量开始探索");
    double xMin_ = -1.0;
    double xMax_ = 1.0;
    double yMin_ = 0.0;
    double yMax_ = 1.0;
    QPoint mousePosition_;
    bool mouseInside_ = false;
};
